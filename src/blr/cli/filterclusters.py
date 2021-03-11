"""
Removes barcode and molecule tags from reads which barcode matches the input set of barcodes.
"""

import logging

from blr.utils import get_bamtag, PySAMIO, Summary, tqdm

logger = logging.getLogger(__name__)


def main(args):
    run_filterclusters(
        input=args.input,
        barcodes=args.barcodes,
        output=args.output,
        barcode_tag=args.barcode_tag,
        molecule_tag=args.molecule_tag
    )


def run_filterclusters(
    input: str,
    barcodes: str,
    output: str,
    barcode_tag: str,
    molecule_tag: str
):
    tags_to_remove = [barcode_tag, molecule_tag]
    removed_tags = {tag: set() for tag in tags_to_remove}
    summary = Summary()
    logger.info("Starting")

    logger.info("Read list of barcodes to filter")
    with open(barcodes, "r") as file:
        barcodes_to_filter = set(file.read().split())

    summary["Barcodes to filter"] = len(barcodes_to_filter)

    logger.info("Filtering BAM")
    with PySAMIO(input, output, __name__) as (openin, openout):
        for read in tqdm(openin.fetch(until_eof=True), desc="Filtering input", unit="reads"):
            summary["Total reads"] += 1

            barcode = get_bamtag(pysam_read=read, tag=barcode_tag)

            if barcode in barcodes_to_filter:
                # Stats
                summary["Removed tags"] += len(tags_to_remove)
                summary["Reads with removed tags"] += 1

                strip_barcode(pysam_read=read, tags_to_be_removed=tags_to_remove, removed_tags=removed_tags)

            summary["Reads written"] += 1
            openout.write(read)

    summary.update({f"Unique {tag} tags removed": len(removed_tags[tag]) for tag in tags_to_remove})

    logger.info("Finished")

    summary.print_stats(name=__name__)


def strip_barcode(pysam_read, tags_to_be_removed, removed_tags):
    """
    Strips an alignment from its barcode sequence. Keeps information in header but adds FILTERED prior to bc info.
    """

    # Modify header
    pysam_read.query_name = f"{pysam_read.query_name}_FILTERED"

    # Remove tags
    for bam_tag in tags_to_be_removed:
        try:
            removed_tags[bam_tag].add(pysam_read.get_tag(bam_tag))
        except KeyError:
            continue

        # Strip read from tag
        pysam_read.set_tag(bam_tag, None, value_type="Z")


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="SAM/BAM file tagged with barcodes information under the tag specified at "
             "-b/--barcode-tag. The file needs to be indexed, sorted & have duplicates marked. "
             "To read from stdin use '-'."
    )
    parser.add_argument(
        "barcodes",
        help="TXT with barcodes to filter out on separate lines."
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Write output BAM to file rather then stdout."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s"
    )
    parser.add_argument(
        "-m", "--molecule-tag", default="MI",
        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
             "Default: %(default)s"
    )
