"""
Parse molecule information from BAM and output stats
"""
import logging
import sys
from contextlib import ExitStack

import pysam
import numpy as np

from blr.cli.buildmolecules import Molecule
from blr.utils import get_bamtag, Summary, tqdm, ACCEPTED_LIBRARY_TYPES, LastUpdatedOrderedDict, calculate_N50

logger = logging.getLogger(__name__)

# This constant should be larger than the maximum possible molecules one could expect. The largest count I have
# seen so far is about 60k so this takes some additional height to this. It could be expanded further if needed.
MAX_MOLECULE_COUNT = 1_000_000_000


def main(args):
    run_readmolecules(
        input=args.input,
        output_tsv=args.output_tsv,
        bed_file=args.bed,
        threshold=args.threshold,
        barcode_tag=args.barcode_tag,
        molecule_tag=args.molecule_tag,
        min_mapq=args.min_mapq,
        library_type=args.library_type
    )


def run_readmolecules(
    input: str,
    output_tsv: str,
    bed_file: str,
    threshold: int,
    barcode_tag: str,
    molecule_tag: str,
    min_mapq: int,
    library_type: str
):

    summary = Summary()
    stats = np.empty((MAX_MOLECULE_COUNT, 2))
    i = 0
    # Read molecules from BAM
    save = pysam.set_verbosity(0)  # Fix for https://github.com/pysam-developers/pysam/issues/939
    with ExitStack() as stack:
        infile = stack.enter_context(pysam.AlignmentFile(input, "rb"))

        # Setup TSV
        if output_tsv is None:
            tsv = sys.stdout
        else:
            tsv = stack.enter_context(open(output_tsv, "w"))
        print("MoleculeID", "Barcode", "Reads", "Length", "BpCovered", sep="\t", file=tsv)

        # Setup BED
        if bed_file:
            bed = stack.enter_context(open(bed_file, "w"))

        for molecule in parse_molecules(openbam=infile, barcode_tag=barcode_tag, molecule_tag=molecule_tag,
                                        library_type=library_type, min_mapq=min_mapq, summary=summary):
            summary["Molecules candidate"] += 1
            if molecule.nr_reads >= threshold:
                summary["Molecules called"] += 1

                stats[i, 0] = molecule.length()
                stats[i, 1] = molecule.bp_covered
                i += 1

                print(molecule.to_tsv(), file=tsv)
                if bed_file:
                    print(molecule.to_bed(), file=bed)

    pysam.set_verbosity(save)

    if i > 0:
        stats.resize((i, 2))
        coverage = 100 * stats[:, 1] / stats[:, 0]
        summary["Fragment N50 (bp)"] = calculate_N50(stats[:, 0])
        summary["Mean fragment size (bp)"] = np.mean(stats[:, 0])
        summary["Median fragment size (bp)"] = np.median(stats[:, 0])
        summary["Longest fragment (bp)"] = np.max(stats[:, 0])
        summary["Mean fragment read coverage (%)"] = np.mean(coverage)
        summary["Median fragment read coverage (%)"] = np.median(coverage)

    summary.print_stats(name=__name__)


def parse_reads(openbam, barcode_tag, molecule_tag, min_mapq, summary):
    for read in tqdm(openbam.fetch(until_eof=True)):
        summary["Total reads"] += 1
        if read.is_duplicate or read.is_unmapped or read.mapping_quality < min_mapq:
            summary["Non analyced reads"] += 1
            continue

        barcode = get_bamtag(pysam_read=read, tag=barcode_tag)
        if not barcode:
            summary["Non analyced reads"] += 1
            continue

        summary[f"Reads with {barcode_tag} tag"] += 1

        molecule_id = get_bamtag(pysam_read=read, tag=molecule_tag)
        if not molecule_id:
            summary["Non analyced reads"] += 1
            continue

        summary[f"Reads with {molecule_tag} tag"] += 1

        yield barcode, molecule_id, read


def parse_molecules(openbam, barcode_tag, molecule_tag, library_type, min_mapq, summary):
    molecules = LastUpdatedOrderedDict()
    prev_chrom = openbam.references[0]
    MAX_DIST = 300_000
    BUFFER_STEP = 1000
    buffer_pos = MAX_DIST + BUFFER_STEP
    for barcode, molecule_id, read in parse_reads(openbam, barcode_tag, molecule_tag, min_mapq, summary):
        if not prev_chrom == read.reference_name:
            yield from molecules.values()
            molecules.clear()
            prev_chrom = read.reference_name

        index = (barcode, molecule_id)
        if index not in molecules:
            molecules[index] = Molecule(read, barcode, index=molecule_id)
        elif molecules[index].has_acceptable_overlap(read, library_type, summary):
            molecules[index].add_read(read)
            molecules.move_to_end(index)

        if read.reference_start > buffer_pos:
            # Leap to current position if ahead of buffer
            buffer_pos = max(read.reference_start, buffer_pos + BUFFER_STEP)

            # Loop over molecules in order of last updated. Get indexes of molecules who are outside MAX_DIST
            # from current position and stop looping if within this distance.
            molcules_to_report = []
            for i, molecule in molecules.items():
                if abs(molecule.stop - read.reference_start) > MAX_DIST:
                    molcules_to_report.append(i)
                else:
                    break

            yield from (molecules.pop(i) for i in molcules_to_report)

    yield from molecules.values()


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="Sorted BAM file tagged with barcode and molecule index."
    )
    parser.add_argument(
        "-o", "--output-tsv",
        help="Write output stats to TSV file. Default: write to stdout."
    )
    parser.add_argument(
        "--bed",
        help="Write molecule bounds to unsorted BED file"
    )
    parser.add_argument(
        "-t", "--threshold", type=int, default=4,
        help="Threshold for how many reads are required for including given molecule in statistics."
             "Default: %(default)s."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
    parser.add_argument(
        "-m", "--molecule-tag", default="MI",
        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
             "Default: %(default)s."
    )
    parser.add_argument(
        "--min-mapq", type=int, default=0,
        help="Minimum mapping-quality to include reads in analysis Default: %(default)s."
    )
    parser.add_argument(
        "-l", "--library-type", default="blr", choices=ACCEPTED_LIBRARY_TYPES,
        help="Select library type from currently available technologies: %(choices)s. Default: %(default)s."
    )
