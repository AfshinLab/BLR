"""
Process stLFR reads with existing barcodes in header. Barcodes of type '21_325_341' where numbers correspond to
barcode sequences are translated to their original sequences(see Note).

E.g.

Header in input:
    @V100002302L1C001R017000001#53_1482_471/1	1	1

Header out output (mapper is ema):
    @V100002302L1C001R017000001:TATTACTCTC-TAGGATAGTT-GATATAGCGG BX:Z:TATTACTCTC-TAGGATAGTT-GATATAGCGG

Header out output (mapper is other):
    @V100002302L1C001R017000001_BX:Z:TATTACTCTC-TAGGATAGTT-GATATAGCGG 1	1

Note: picard MarkDuplicates in barcode-aware mode has to have the barcode match the regex: ^[ATCGNatcgn-]*$."
"""

import logging
import sys
import dnaio
from tqdm import tqdm
from collections import Counter

from blr.utils import print_stats

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting")

    summary = Counter()

    index_to_barcode = parse_barcodes(args.barcodes)

    in_interleaved = not args.input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    # If no output1 is given output is sent to stdout
    if not args.output1:
        logger.info("Writing output to stdout.")
        args.output1 = sys.stdout.buffer
        args.output2 = None

    out_interleaved = not args.output2
    logger.info(f"Output detected as {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with dnaio.open(args.input1, file2=args.input2, interleaved=in_interleaved, mode="r",
                    fileformat="fastq") as reader, \
            dnaio.open(args.output1, file2=args.output2, interleaved=out_interleaved, mode="w",
                       fileformat="fastq") as writer:
        for read1, read2 in tqdm(reader, desc="Read pairs processed"):
            summary["Reads read"] += 1
            # Header parsing
            name, r1_rest = read1.name.split(maxsplit=1)
            _, r2_rest = read2.name.split(maxsplit=1)

            name, barcode_indeces = name.strip("/1").split("#")

            barcode = translate_indeces(barcode_indeces, index_to_barcode, summary)

            if barcode:
                barcode_id = f"{args.barcode_tag}:Z:{barcode}"

                if args.mapper == "ema":
                    # The EMA aligner requires reads in 10x format e.g.
                    # @READNAME:AAAAAAAATATCTACGCTCA BX:Z:AAAAAAAATATCTACGCTCA
                    new_name = ":".join([name, barcode])
                    new_name = " ".join((new_name, barcode_id))
                    read1.name, read2.name = new_name, new_name
                else:
                    new_name = "_".join([name, barcode_id])
                    read1.name = " ".join([new_name, r1_rest])
                    read2.name = " ".join([new_name, r2_rest])
            else:
                summary["Read pairs missing barcode"] += 1

                # Skip read missing barcode for ema mapping.
                if args.mapper == "ema":
                    continue

            # Write to out
            writer.write(read1, read2)
            summary["Reads written"] += 1
    print(summary)
    print_stats(summary)
    logger.info("Finished")


def translate_indeces(indeces, index_to_barcode, summary):
    if indeces != "0_0_0":  # stLFR reads are tagged with #0_0_0 if the barcode could not be identified.
        barcodes = [index_to_barcode[int(i)] for i in indeces.split("_") if i != ""]
        summary[f"Barcodes of length {len(barcodes)}"] += 1
        if barcodes:
            return "-".join(barcodes)
    else:
        summary["Skipped barcode type 0_0_0"] += 1
    return None


def parse_barcodes(file):
    index_to_barcode = dict()
    with open(file) as f:
        for line in f:
            barcode, index = line.strip().split(maxsplit=1)
            index_to_barcode[int(index)] = barcode

    return index_to_barcode


def add_arguments(parser):
    parser.add_argument(
        "barcodes",
        help="stLFR barocode list for tab separated barcode sequences and indexes."
    )
    parser.add_argument(
        "input1",
        help="Input FASTQ/FASTA file. Assumes to contain read1 if given with second input file. "
             "If only input1 is given, input is assumed to be an interleaved. If reading from stdin"
             "is requested use '-' as a placeholder.")
    parser.add_argument(
        "input2", nargs='?',
        help="Input  FASTQ/FASTA for read2 for paired-end read. Leave empty if using interleaved.")
    parser.add_argument(
        "--output1", "--o1",
        help="Output FASTQ/FASTA file name for read1. If not specified the result is written to "
             "stdout as interleaved. If output1 given but not output2, output will be written as "
             "interleaved to output1.")
    parser.add_argument(
        "--output2", "--o2",
        help="Output FASTQ/FASTA name for read2. If not specified but --o1/--output1 given the "
             "result is written as interleaved .")
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument(
        "-m", "--mapper",
        help="Specify read mapper for labeling reads with barcodes. "
    )
