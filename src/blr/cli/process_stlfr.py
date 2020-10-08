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
from contextlib import ExitStack

from blr.utils import print_stats
from blr.cli.tagfastq import Output, ChunkHandler

logger = logging.getLogger(__name__)


def main(args):
    run_process_stlfr(
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
        barcodes=args.barcodes,
        barcode_tag=args.barcode_tag,
        mapper=args.mapper,
    )


def run_process_stlfr(
        input1: str,
        input2: str,
        output1: str,
        output2: str,
        barcodes: str,
        barcode_tag: str,
        mapper: str,
):
    logger.info("Starting")

    summary = Counter()

    # TODO This translation might not acctually be true since its unsure how the barcode file relates to the tagged
    #  indeces on the FASTQs. Instead one could generate unique barcodes in for each index combo. This would aslo solve
    #  https://github.com/FrickTobias/BLR/issues/219 for stLFR reads
    index_to_barcode = {index: barcode for barcode, index in parse_barcodes(barcodes)}

    in_interleaved = not input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    # If no output1 is given output is sent to stdout
    if not output1:
        logger.info("Writing output to stdout.")
        output1 = sys.stdout.buffer
        output2 = None

    out_interleaved = not output2
    logger.info(f"Output detected as {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with ExitStack() as stack:
        reader = stack.enter_context(dnaio.open(input1, file2=input2, interleaved=in_interleaved, mode="r"))
        writer = stack.enter_context(Output(output1, output2, interleaved=out_interleaved, mapper=mapper))
        if mapper in ["ema", "lariat"]:
            chunks = stack.enter_context(ChunkHandler(chunk_size=1_000))
            heaps = BarcodeHeap()

        for read1, read2 in tqdm(reader, desc="Read pairs processed"):
            summary["Read pairs read"] += 1

            name = read1.name.split("\t")[0]

            # Remove '/1' from read name and split to get barcode_indeces
            name, barcode_indeces = name.strip("/1").split("#")

            barcode = translate_indeces(barcode_indeces, index_to_barcode, summary)

            if barcode:
                barcode_id = f"{barcode_tag}:Z:{barcode}"
                if mapper == "ema":
                    # The EMA aligner requires reads in 10x format e.g.
                    # @READNAME:AAAAAAAATATCTACGCTCA BX:Z:AAAAAAAATATCTACGCTCA
                    name = f"{name}:{barcode} {barcode_id}"
                elif mapper != "lariat":
                    name = f"{name}_{barcode_id}"
            else:
                summary["Read pairs missing barcode"] += 1
                if mapper in ["ema", "lariat"]:
                    continue

            read1.name = name
            read2.name = name

            # Write to out
            if mapper == "ema":
                chunks.build_chunk(
                    heaps.get_heap(barcode),
                    read1.name,
                    read1.sequence,
                    read1.qualities,
                    read2.name,
                    read2.sequence,
                    read2.qualities,
                )
            elif mapper == "lariat":
                corrected_barcode_qual = "K" * len(barcode)
                chunks.build_chunk(
                    heaps.get_heap(barcode),
                    f"@{name}",
                    read1.sequence,
                    read1.qualities,
                    read2.sequence,
                    read2.qualities,
                    f"{barcode}-1",
                    corrected_barcode_qual,
                    "AAAAAA",
                    "KKKKKK"
                )
            else:
                summary["Read pairs written"] += 1
                writer.write(read1, read2)

        # Empty final cache
        if mapper == "ema":
            chunks.write_chunk()

            for entry in chunks.parse_chunks():
                r1 = dnaio.Sequence(*entry[1:4])
                r2 = dnaio.Sequence(*entry[4:7])
                writer.write(r1, r2)
                summary["Read pairs written"] += 1

        elif mapper == "lariat":
            chunks.write_chunk()

            for entry in chunks.parse_chunks():
                print(*entry[1:10], sep="\n", file=writer)
                summary["Read pairs written"] += 1

    print_stats(summary, __name__)
    logger.info("Finished")


def translate_indeces(indeces, index_to_barcode, summary):
    if indeces == "0_0_0":  # stLFR reads are tagged with #0_0_0 if the barcode could not be identified.
        summary["Skipped barcode type 0_0_0"] += 1
        return None
    else:
        # Translate index to barcode sequence, index could be missing e.g. "12_213_" which leades to a empty string
        barcodes = [index_to_barcode[int(i)] for i in indeces.split("_") if i != ""]
        summary[f"Barcodes of length {len(barcodes)}"] += 1
        # TODO Currently partial barcodes are kept, possibly only full 3-part barcodes should be included.
        return "".join(barcodes) if barcodes else None


def parse_barcodes(file):
    with open(file) as f:
        for line in f:
            barcode, index = line.strip().split(maxsplit=1)
            yield barcode, int(index)


class BarcodeHeap:
    def __init__(self):
        self._heap = 0
        self._barcodes = {}

    def get_heap(self, barcode: str) -> str:
        heap = self._barcodes.setdefault(barcode, self._heap + 1)
        if heap > self._heap:
            self._heap += 1
        return str(heap)


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
