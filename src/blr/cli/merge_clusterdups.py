"""
Merge barcodes by re-tagging reads in BAM with new barcodes.
"""

import logging

from blr.utils import PySAMIO, get_bamtag, Summary, tqdm

logger = logging.getLogger(__name__)


def main(args):
    run_mergeclusters(
        input=args.input,
        input_merges=args.input_merges,
        output=args.output,
        barcode_tag=args.barcode_tag,
    )


def run_mergeclusters(
    input: str,
    input_merges: str,
    output: str,
    barcode_tag: str,
):
    summary = Summary()

    # Create mapping from old to new barcode from CSV with old-new pairs.
    with open(input_merges) as file:
        old_to_new_barcode = dict(tuple(line.strip().split(",")) for line in file)

    with PySAMIO(input, output, __name__) as (infile, out):
        for read in tqdm(infile, desc="Writing output", total=summary["Total reads"]):
            old_barcode = get_bamtag(pysam_read=read, tag=barcode_tag)
            if old_barcode in old_to_new_barcode:
                summary["Reads with new barcode"] += 1
                new_barcode = old_to_new_barcode[old_barcode]
                read.set_tag(barcode_tag, new_barcode, value_type="Z")

            out.write(read)

    logger.info("Finished")
    summary.print_stats(name=__name__)


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="Coordinate-sorted SAM/BAM file tagged with barcodes."
    )
    parser.add_argument(
        "input_merges",
        help="CSV log file containing all merges to be done. File is in format: {old barcode},{new barcode}."
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Write output BAM to file rather then stdout."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
