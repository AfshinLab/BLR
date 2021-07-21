"""
Strips headers from tags and depending on mode, set the appropriate SAM tag.
"""

import logging

from blr.utils import Summary, PySAMIO, get_bamtag, tqdm, ACCEPTED_READ_MAPPERS

logger = logging.getLogger(__name__)

DNA_BASES = {"A", "T", "C", "G"}


def main(args):
    run_tagbam(
        input=args.input,
        output=args.output,
        mapper=args.mapper,
        sample_number=args.sample_nr,
        barcode_tag=args.barcode_tag,
    )


def run_tagbam(
        input: str,
        output: str,
        mapper: str,
        sample_number: int,
        barcode_tag: str,
):
    logger.info("Starting analysis")

    if mapper == "ema":
        processing_function = mode_ema
    elif mapper == "lariat":
        processing_function = lambda *args: None
    else:
        processing_function = mode_samtags_underline_separation

    summary = Summary()

    # Read SAM/BAM files and transfer barcode information from alignment name to SAM tag
    with PySAMIO(input, output, __name__) as (infile, outfile):
        for read in tqdm(infile.fetch(until_eof=True), desc="Processing reads", unit=" reads"):
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
        "-m", "--mapper", default="bowtie2", choices=ACCEPTED_READ_MAPPERS,
        help="Mapper used for aligning reads. Default: %(default)s."
    )
    parser.add_argument(
        "-s", "--sample-nr", default=1, type=int,
        help="Add sample number to each barcode. Default: %(default)s."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
