"""
Strips headers from tags and depending on mode, set the appropriate SAM tag.
"""

import logging
from collections import Counter
from pathlib import Path

from blr.utils import print_stats, PySAMIO, get_bamtag, tqdm
from blr.cli.config import load_yaml

logger = logging.getLogger(__name__)

CONFIGFILE = Path("blr.yaml")


def main(args):
    logger.info("Starting analysis")

    # Load configs to find mapper
    configs, _ = load_yaml(CONFIGFILE)
    mapper = configs["read_mapper"]

    if mapper == "ema":
        processing_function = mode_ema
    else:
        processing_function = mode_samtags_underline_separation

    summary = Counter()

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with PySAMIO(args.input, args.output, __name__) as (infile, outfile):
        for read in tqdm(infile.fetch(until_eof=True), desc="Processing reads", unit=" reads"):
            # Strips header from tag and depending on script mode, possibly sets SAM tag
            summary["Total reads"] += 1
            processing_function(read, args.sample_nr, summary)
            outfile.write(read)

    print_stats(summary, name=__name__)
    logger.info("Finished")


def mode_samtags_underline_separation(read, sample_nr, summary):
    """
    Trims header from tags and sets SAM tags according to values found in header.
    Assumes format: @header_<tag>:<type>:<seq> (can be numerous tags). Constrictions are: Header includes SAM tags
    separated by "_".
    :param read: pysam read alignment
    :param sample_nr: barcodes samples tag.
    :param summary: Collections's Counter object
    :return:
    """

    # Strip header
    header = read.query_name.split("_")
    read.query_name = header[0]

    # Set SAM tags
    for tag in header[1:]:
        tag, tag_type, val = tag.split(":")
        if tag == "BX":
            val += "-" + str(sample_nr)

        read.set_tag(tag, val, value_type=tag_type)
        summary[f"Reads with tag {tag}"] += 1


def mode_ema(read, sample_nr, _):  # summary is passed to this function but is not used
    """
    Trims header from barcode sequences.
    Assumes format @header:and:more...:header:<seq>. Constrictions: There must be exactly 9 elements separated by ":"
    :param read: pysam read alignment
    :param sample_nr:
    :return:
    """

    # Split header into original read name and barcode
    read.query_name, header_barcode = read.query_name.rsplit(":", 1)

    # Modify tag barcode to remove '-1' added at end by ema e.g 'BX:Z:TTTGTTCATGAGTACG-1' --> 'BX:Z:TTTGTTCATGAGTACG'
    # Ema also trims the barcode to 16bp (10x Barcode length) so it need to be exchanged for the one in the header.
    current_barcode = get_bamtag(read, "BX")
    if current_barcode and current_barcode.endswith("-1"):
        modified_barcode = current_barcode[:-2]
        # Make sure that the SAM tag barcode is a substring of the header barcode
        assert header_barcode.startswith(modified_barcode)
        read.set_tag("BX", header_barcode + "-" + str(sample_nr), value_type="Z")


def add_arguments(parser):
    parser.add_argument("input",
                        help="BAM file with SAM tag info in header. To read from stdin use '-'.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("-s", "--sample-nr", default=1, type=int,
                        help="Add sample number to each barcode ")
