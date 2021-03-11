"""
Correct single count barcodes by finding matching multiple count barcodes within hamming distance of one.

This is the method used in Chen et al. 2019 (doi: 10.1101/gr.260380.119) to correct TELL-seq barcodes which are
partialy degenerate.
"""
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import logging
import sys

import dnaio

from blr.utils import tqdm, Summary


logger = logging.getLogger(__name__)


def main(args):
    summary = Summary()

    barcodes = count_barcodes(args.uncorrected_barcodes)
    summary["Total barcodes "] += len(barcodes)
    singles, multiples = split_by_count(barcodes, count=1)

    summary["Barcodes with count > 1"] = len(multiples)
    summary["Barcodes with count = 1"] = len(singles)

    correct_singles(singles, multiples, summary)

    summary["Corrected singles (%)"] = 100*summary["Corrected singles"] / summary["Barcodes with count = 1"]

    with open_output(args.output) as output:
        for barcode, cluster in sorted(multiples.items(), key=lambda x: x[1].count, reverse=True):
            print(barcode, cluster.count, ",".join(cluster.barcodes), sep="\t", file=output)

        for barcode in singles:
            print(barcode, 1, barcode, sep="\t", file=output)

    summary.print_stats(__name__)


@contextmanager
def open_output(output_file):
    if output_file == "-":
        yield sys.stdout
    else:
        yield open(output_file, "w")


def count_barcodes(file):
    with dnaio.open(file) as barcodes:
        return Counter([barcode.sequence for barcode in tqdm(barcodes, desc="Count barcodes")])


@dataclass
class Cluster:
    count: int
    barcodes: list

    def add(self, barcode):
        self.count += 1
        self.barcodes.append(barcode)


def split_by_count(barcodes, count=1):
    singles = {barcode for barcode, count in barcodes.items() if count == 1}
    multiples = {barcode: Cluster(count, [barcode]) for barcode, count in barcodes.items() if count > 1}
    return singles, multiples


def correct_singles(singles, multiples, summary):
    nucleotides = {"A", "T", "C", "G"}
    singles_corrected = set()
    for sequence in tqdm(singles, desc="Correcting singles"):
        seq_list = list(sequence)
        matched = None
        found = 0
        for i, base in enumerate(seq_list):
            alt_seq_list = seq_list.copy()
            subsitutions = nucleotides - set(base)
            for sub_base in subsitutions:
                alt_seq_list[i] = sub_base
                alt_sequence = ''.join(alt_seq_list)
                if alt_sequence in multiples:
                    matched = alt_sequence
                    found += 1

        if found == 1:
            multiples[matched].add(sequence)
            singles_corrected.add(sequence)
            summary["Corrected singles"] += 1
        elif found > 1:
            summary["Singles matching multiple"] += 1

    singles -= singles_corrected


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes.")
    parser.add_argument(
        "-o", "--output", default="-",
        help="Output tab separated file for error corrected barcodes. Columns are: corrected barcode sequence, count, "
             "comma-separate barcodes that were corrected to the sequence Default: write to stdout.")
