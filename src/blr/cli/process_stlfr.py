"""
Process stLFR reads with existing barcodes in header.

Barcodes of type '21_325_341' where numbers correspond to barcode sequences are translated to unique DNA strings-type
barcode of length 16 nt (similar to 10x barcode format)
"""

from contextlib import ExitStack
from itertools import product
import logging
import sys
from pathlib import Path

import dnaio

from blr.utils import Summary, ACCEPTED_READ_MAPPERS, tqdm
from blr.cli.tagfastq import Output, ChunkHandler, write_ema_output, write_lariat_output

logger = logging.getLogger(__name__)


def main(args):
    run_process_stlfr(
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
        output_nobc1=args.output_nobc1,
        output_nobc2=args.output_nobc2,
        output_bins=args.output_bins,
        nr_bins=args.nr_bins,
        output_translations=args.output_translations,
        barcode_tag=args.barcode_tag,
        mapper=args.mapper,
        sample_number=args.sample_nr,
    )


def run_process_stlfr(
        input1: str,
        input2: str,
        output1: str,
        output2: str,
        output_nobc1: str,
        output_nobc2: str,
        output_bins: str,
        nr_bins: int,
        output_translations: str,
        barcode_tag: str,
        mapper: str,
        sample_number: int,
):
    logger.info("Starting")

    summary = Summary()

    barcodes = BarcodeGenerator()

    in_interleaved = not input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    special_fmt = False
    out_interleaved = not output2 or not output1
    if output_bins is not None:
        logger.info(f"Writing output as binned interleaved FASTQ to {output_bins}.")
        output_bins = Path(output_bins)
        output_bins.mkdir(exist_ok=True)
        output1, output2 = None, None
        special_fmt = mapper == "ema"
    else:
        if not output1:
            logger.info("Writing output to stdout.")
            output1 = sys.stdout.buffer
            output2 = None
        logger.info(f"Output detected as {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    if output_nobc1 is not None and mapper != "ema":
        logger.warning(f"Writing non barcoded reads to {output_nobc1} is only available with option '--mapper ema'.")
        output_nobc1, output_nobc2 = None, None

    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with ExitStack() as stack:
        reader = stack.enter_context(dnaio.open(input1, file2=input2, interleaved=in_interleaved, mode="r"))
        writer = stack.enter_context(Output(file1=output1, file2=output2, interleaved=out_interleaved, mapper=mapper,
                                            file_nobc1=output_nobc1, file_nobc2=output_nobc2, bins_dir=output_bins,
                                            nr_bins=nr_bins))
        chunks = None
        if mapper == "lariat" or (mapper == "ema" and output_bins is None):
            chunks = stack.enter_context(ChunkHandler(chunk_size=1_000_000))
            heaps = BarcodeHeap()

        for read1, read2, barcode in parse_stlfr_reads(reader, barcodes, barcode_tag, mapper, summary, special_fmt):
            if barcode is None:
                summary["Read pairs missing barcode"] += 1

                # Write non barcoded reads to separate file if exists for ema.
                if mapper == "ema" and output_nobc1 is not None:
                    summary["Read pairs written"] += 1
                    writer.write_nobc(read1, read2)
                    continue

                if mapper in ["ema", "lariat"]:
                    continue

            # Write to out
            if special_fmt:
                writer.write_ema_special(read1, read2, barcode)
                summary["Read pairs written"] += 1
            elif mapper == "ema":
                chunks.build_chunk(
                    f"{heaps.get_heap(barcode)}\t"
                    f"{read1.name}\t"
                    f"{read1.sequence}\t"
                    f"{read1.qualities}\t"
                    f"{read2.sequence}\t"
                    f"{read2.qualities}\n"
                )
            elif mapper == "lariat":
                chunks.build_chunk(
                    f"{heaps.get_heap(barcode)}\t"
                    f"@{read1.name}\t"
                    f"{read1.sequence}\t"
                    f"{read1.qualities}\t"
                    f"{read2.sequence}\t"
                    f"{read2.qualities}\t"
                    f"{barcode}-{sample_number}\t"
                    f"KKKKKKKKKKKKKKKK\t"
                    "AAAAAA\t"
                    "KKKKKK\n"
                )
            else:
                summary["Read pairs written"] += 1
                writer.write(read1, read2)

        # Empty final cache
        if chunks is not None:
            chunks.write_chunk()

        if mapper == "lariat" and output_bins is not None:
            remaining_reads = summary["Read pairs read"] - summary["Reads missing barcode"] - summary["Read pairs written"]  # noqa: E501
            bin_size = remaining_reads // nr_bins
            logger.info(f"Using bin of size {bin_size}.")
            writer.set_bin_size(bin_size)

        if mapper == "ema" and output_bins is None:
            write_ema_output(chunks, writer, summary)
        elif mapper == "lariat":
            write_lariat_output(chunks, writer, summary)

    if output_translations is not None:
        with open(output_translations, "w") as f:
            for index, barcode in barcodes.translate_barcode.items():
                print(f"{index},{barcode}", file=f)

    summary.print_stats(__name__)
    logger.info("Finished")


def parse_stlfr_reads(reader, barcodes, barcode_tag, mapper, summary, special_fmt):
    for read1, read2 in tqdm(reader, desc="Read pairs processed"):
        summary["Read pairs read"] += 1

        name = read1.name.split("\t")[0]

        # Remove '/1' from read name and split to get barcode_indices
        name, barcode_indices = name.strip("/1").split("#")

        barcode = translate_indeces(barcode_indices, barcodes, summary)

        if barcode is not None:
            barcode_id = f"{barcode_tag}:Z:{barcode}"
            if special_fmt:
                pass
            elif mapper == "ema":
                # The EMA aligner requires reads in 10x FASTQ format e.g.
                # @READNAME:AAAAAAAATATCTACGCTCA BX:Z:AAAAAAAATATCTACGCTCA
                name = f"{name}:{barcode} {barcode_id}"
            elif mapper != "lariat":
                name = f"{name}_{barcode_id}"

        read1.name = name
        read2.name = name
        yield read1, read2, barcode


def translate_indeces(index_string, barcodes, summary):
    if index_string == "0_0_0":  # stLFR reads are tagged with 0_0_0 if the barcode could not be identified.
        summary["Skipped barcode type 0_0_0"] += 1
        return None
    elif len(index_string.split("_")) < 3:
        summary["Skipped barcode of incorrect length"] += 1
        return None
    return barcodes.get(index_string)


class BarcodeGenerator:
    """
    Generate 16 nt unique barcode barcodes for particular stLFR index string, e.g. '1_2_4'
     """
    def __init__(self,):
        self._barcode_generator = self._generate_barcodes()
        # Skip first entry which is AAAAAAAAAAAAAAAA. Reads tagged with this barcode causes
        # error `Segmentation fault: 11` in ema.
        _ = next(self._barcode_generator)
        self.translate_barcode = {}

    def get(self, index_string):
        if index_string in self.translate_barcode:
            return self.translate_barcode[index_string]

        barcode = next(self._barcode_generator)
        self.translate_barcode[index_string] = barcode
        return barcode

    @staticmethod
    def _generate_barcodes():
        """
        Iterator to generate unique 16 nt barcode
        """
        # Note: Theortical barcode numbers.
        #   Possible stLFR barcodes: 1536*1536*1536 = 3,623,878,656
        #   Possible barcodes of len 16: 4^16 =       4,294,967,296
        for barcode in product("ATCG", repeat=16):
            yield "".join(barcode)


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
        "input1",
        help="Input FASTQ/FASTA file. Assumes to contain read1 if given with second input file. "
             "If only input1 is given, input is assumed to be an interleaved. If reading from stdin"
             "is requested use '-' as a placeholder."
    )
    parser.add_argument(
        "input2", nargs='?',
        help="Input  FASTQ/FASTA for read2 for paired-end read. Leave empty if using interleaved."
    )
    output = parser.add_mutually_exclusive_group(required=False)
    output.add_argument(
        "--output1", "--o1",
        help="Output FASTQ/FASTA file name for read1. If not specified the result is written to "
             "stdout as interleaved. If output1 given but not output2, output will be written as "
             "interleaved to output1."
    )
    parser.add_argument(
        "--output2", "--o2",
        help="Output FASTQ/FASTA name for read2. If not specified but --o1/--output1 given the "
             "result is written as interleaved ."
    )
    parser.add_argument(
        "--output-nobc1", "--n1",
        help="Only for ema! Output FASTQ/FASTA file name to write non-barcoded read1 reads. "
             "If output_nobc1 given but not output_nobc2, output will be written as interleaved to "
             "output_nobc1."
    )
    parser.add_argument(
        "--output-nobc2", "--n2",
        help="Only for ema! Output FASTQ/FASTA file name to write non-barcoded read2 reads."
    )
    output.add_argument(
        "--output-bins",
        help=f"Output barcoded reads to bins named '{Output.BIN_FASTQ_TEMPLATE}' in the provided directory. Only "
             f"used for ema mapping and uses ema special format."
    )
    parser.add_argument(
        "--nr-bins", type=int, default=100,
        help="Number of bins to split reads into when using the '--output-bins' alternative. Default: %(default)s."
    )
    parser.add_argument(
        "-t", "--output-translations",
        help="Output CSV with combinatorial barcode index to translated barcode sequence."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
    parser.add_argument(
        "-m", "--mapper", default="bowtie2", choices=ACCEPTED_READ_MAPPERS,
        help="Specify read mapper for labeling reads with barcodes. Default: %(default)s."
    )
    parser.add_argument(
        "--sample-nr", type=int, default=1,
        help="Sample number to append to barcode string. Default: %(default)s."
    )
