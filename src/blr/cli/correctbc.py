"""
Correct single count barcodes by finding matching multiple count barcodes within hamming distance of one.

This is the method used in Chen et al. 2019 (doi: 10.1101/gr.260380.119) to correct TELL-seq barcodes which are
partialy degenerate.
"""
from collections import Counter
from dataclasses import dataclass
import logging
from typing import List

import dnaio

from blr.utils import tqdm, Summary, smart_open


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

    with smart_open(args.output) as output:
        for barcode, cluster in sorted(multiples.items(), key=lambda x: x[1], reverse=True):
            print(cluster, file=output)

        for barcode in singles:
            print(f"{barcode}\t1\t{barcode}", file=output)

    summary.print_stats(__name__)


def count_barcodes(file):
    with dnaio.open(file) as barcodes:
        return Counter([barcode.sequence for barcode in tqdm(barcodes, desc="Count barcodes")])


@dataclass(order=False)
class Cluster:
    count: int
    barcodes: List[str]

    def add(self, barcode: str):
        self.count += 1
        self.barcodes.append(barcode)

    def __lt__(self, other):
        return self.count < other.count

    def __str__(self):
        return f"{self.barcodes[0]}\t{self.count}\t{','.join(self.barcodes)}"


def split_by_count(barcodes, count=1):
    singles = {barcode for barcode, count in barcodes.items() if count == 1}
    multiples = {barcode: Cluster(count, [barcode]) for barcode, count in barcodes.items() if count > 1}
    return singles, multiples


def correct_singles(singles, multiples, summary):
    nucleotides = {"A", "T", "C", "G"}
    singles_corrected = set()
    for sequence in tqdm(singles, desc="Correcting singles"):
        matched = [mut_seq for mut_seq in mutate(sequence, nucleotides) if mut_seq in multiples]

        if len(matched) == 1:
            multiples[matched[0]].add(sequence)
            singles_corrected.add(sequence)
            summary["Corrected singles"] += 1
        elif len(matched) > 1:
            summary["Singles matching multiple"] += 1

    singles -= singles_corrected


def mutate(sequence, nucleotides):
    """Generate all strings within Hamming distance of 1"""
    seq_list = list(sequence)
    for i, base in enumerate(seq_list):
        alt_seq_list = seq_list.copy()
        subsitutions = nucleotides - set(base)
        for sub_base in subsitutions:
            alt_seq_list[i] = sub_base
            yield ''.join(alt_seq_list)


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes."
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Output tab separated file for error corrected barcodes. Columns are: corrected barcode sequence, count, "
             "comma-separate barcodes that were corrected to the sequence. Default: write to stdout."
    )
