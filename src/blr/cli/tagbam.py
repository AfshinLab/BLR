"""
Strips headers from tags and depending on mode, set the appropriate SAM tag.
"""

import logging
from itertools import chain
import re

from blr.utils import Summary, PySAMIO, get_bamtag, tqdm

logger = logging.getLogger(__name__)

DNA_BASES = {"A", "T", "C", "G"}


def main(args):
    run_tagbam(
        input=args.input,
        output=args.output,
        sample_number=args.sample_nr,
        barcode_tag=args.barcode_tag,
    )


def get_mode(parser, barcode_tag: str):
    samtags_underline_pattern = re.compile(r"_[A-Z]{2}:[AifZHB]:.*(_|$)")
    barcodes_end_pattern = re.compile(r":[ATGC]{17,}$")
    processing_function = mode_void

    first = []
    # Check first 1000 reads for format
    for i in range(1000):
        try:
            read = next(parser)
        except StopIteration:
            logger.info("Using mode_void")
            break

        first.append(read)
        if samtags_underline_pattern.search(read.query_name):
            logger.info("Using mode_samtags_underline_separation")
            processing_function = mode_samtags_underline_separation
            break
        elif barcodes_end_pattern.search(read.query_name) and read.has_tag(barcode_tag):
            logger.info("Using mode_ema")
            processing_function = mode_ema
            break
    else:  # found nothing
        logger.info("Using mode_void")

    return processing_function, chain(first, parser)


def run_tagbam(
        input: str,
        output: str,
        sample_number: int,
        barcode_tag: str,
):
    logger.info("Starting analysis")

    summary = Summary()

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with PySAMIO(input, output, __name__) as (infile, outfile):
        parser = infile.fetch(until_eof=True)
        processing_function, parser = get_mode(parser, barcode_tag=barcode_tag)

        for read in tqdm(parser, desc="Processing reads", unit=" reads"):
            # Strips header from tag and depending on script mode, possibly sets SAM tag
            summary["Total reads"] += 1
            processing_function(read, sample_number, barcode_tag, summary)
            outfile.write(read)

    summary.print_stats(name=__name__)
    logger.info("Finished")


def mode_samtags_underline_separation(read, sample_nr, barcode_tag, summary):
    """
    Trims tag strings from header and sets them as SAM tags according to values found in header.
    Header should have suffix of SAM tag(s) separeted by underline ('_') e.g.

        header:name:with:stuff_<tag1>:<type1>:<seq1>_<tag2>:<type2>:<seq2>....
    """
    # Strip header
    header = read.query_name.split("_")
    read.query_name = header[0]

    # Set SAM tags
    for tag in header[1:]:
        tag, tag_type, val = tag.split(":")

        # Input from 10x has "-1" attached to the barcode which needs to be removed.
        val = val.split("-")[0]
        assert is_sequence(val)

        if tag == barcode_tag:
            val = f"{val}-{sample_nr}"

        read.set_tag(tag, val, value_type=tag_type)
        summary[f"Reads with tag {tag}"] += 1


def mode_ema(read, sample_nr, barcode_tag, _):  # summary is passed to this function but is not used
    """
    Extract barcode from read header and replace trunkated barcode in SAM tag placed by EMA. Non-barcoded reads
    left intact. Assumes header format as below

        header:and:more...:header:<seq>
    """
    # Check if read is barcoded before doing correction
    tag_barcode = get_bamtag(read, barcode_tag)
    if tag_barcode is not None:
        # Split header into original read name and barcode and check that the header barcode is valid
        read.query_name, header_barcode = read.query_name.rsplit(":", 1)
        assert is_sequence(header_barcode)

        # Remove '-<sample_nr>' added at end by ema e.g 'TTTGTTCATGAGTACG-1' --> 'TTTGTTCATGAGTACG'
        tag_barcode = tag_barcode[:16]

        # Ema also trims the barcode to 16bp (10x Barcode length) so it need to be exchanged for the one in the header.
        # Make sure that the SAM tag barcode is a substring of the header barcode
        assert header_barcode.startswith(tag_barcode)
        read.set_tag(barcode_tag, f"{header_barcode}-{sample_nr}", value_type="Z")


def mode_void(*args, **kwargs):
    pass


def is_sequence(string: str) -> bool:
    """Check if string is DNA sequence"""
    return set(string).issubset(DNA_BASES)


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="BAM file with SAM tag info in header. To read from stdin use '-'."
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Write output BAM to file rather then stdout."
    )
    parser.add_argument(
        "-s", "--sample-nr", default=1, type=int,
        help="Add sample number to each barcode. Default: %(default)s."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
