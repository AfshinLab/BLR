"""
Tag FASTQ/FASTA headers with corrected barcodes from separate raw barcode FASTQ with matching
headers and starcode CLSTR output file with corrected barcodes.

ABOUT:

First the raw barcodes FASTQ is parser to get a dictionary:

    raw_barcodes[<HEADER>] = <RAW_BARCODE>

Then the corrected barcodes CLSTR file from starcode are parsed to get a dictionary:

    corrected_barcodes[<RAW_BARCODE>] = <CORRECTED_BARCODE>

The for each read-pair in the input FASTQ(s) the corrected barcode is recovered and used to tag
the read by including it in the header.

    <HEADER> ==> <RAW_BARCODE> ==> <CORRECTED_BARCODE>
"""

import logging
import sys
import dnaio
from xopen import xopen
from itertools import islice
from collections import Counter
import heapq
from pathlib import Path
import tempfile
from contextlib import ExitStack

from blr.utils import tqdm, print_stats

logger = logging.getLogger(__name__)


def main(args):
    run_tagfastq(
        uncorrected_barcodes=args.uncorrected_barcodes,
        corrected_barcodes=args.corrected_barcodes,
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
        barcode_tag=args.barcode_tag,
        sequence_tag=args.sequence_tag,
        mapper=args.mapper,
        skip_singles=args.skip_singles
    )


def run_tagfastq(
        uncorrected_barcodes: str,
        corrected_barcodes: str,
        input1: str,
        input2: str,
        output1: str,
        output2: str,
        barcode_tag: str,
        sequence_tag: str,
        mapper: str,
        skip_singles: bool
):
    logger.info("Starting")
    summary = Counter()
    # Get the corrected barcodes and create a dictionary pointing each raw barcode to its
    # canonical sequence.
    with open(corrected_barcodes, "r") as reader:
        corrected_barcodes, heap = parse_corrected_barcodes(reader, summary, mapper, skip_singles=skip_singles)

    in_interleaved = not input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    # If no output1 is given output is sent to stdout
    if not output1:
        logger.info("Writing output to stdout.")
        output1 = sys.stdout.buffer
        output2 = None

    out_interleaved = not output2
    logger.info(f"Output detected as {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    if mapper == "ema":
        # Required for ema sorted output
        # Inspired by: https://stackoverflow.com/questions/56948292/python-sort-a-large-list-that-doesnt-fit-in-memory
        output_chunk = list()
        chunk_size = 1_000_000
        logger.info(f"Using chunks of size {chunk_size} for sorting.")
        chunk_id = 0
        tmpdir = Path(tempfile.mkdtemp(prefix="tagfastq_sort"))
        chunk_file_template = "chunk_*.tsv"
        chunk_sep = "\t"
        tmp_writer = open(tmpdir / chunk_file_template.replace("*", str(chunk_id)), "w")

    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with dnaio.open(input1, file2=input2, interleaved=in_interleaved, mode="r",
                    fileformat="fastq") as reader, \
            Output(output1, output2, interleaved=out_interleaved, mapper=mapper) as writer, \
            BarcodeReader(uncorrected_barcodes) as uncorrected_barcode_reader:

        for read1, read2 in tqdm(reader, desc="Read pairs processed", disable=False):
            # Header parsing
            summary["Read pairs read"] += 1
            name_and_pos, nr_and_index1 = read1.name.split(maxsplit=1)
            _, nr_and_index2 = read2.name.split(maxsplit=1)

            sample_index = nr_and_index1.split(':')[-1]
            uncorrected_barcode_seq = uncorrected_barcode_reader.get_barcode(name_and_pos)
            corrected_barcode_seq = None

            # Check if barcode was found and update header with barcode info.
            if uncorrected_barcode_seq and uncorrected_barcode_seq in corrected_barcodes:
                corrected_barcode_seq = corrected_barcodes[uncorrected_barcode_seq]

                raw_barcode_id = f"{sequence_tag}:Z:{uncorrected_barcode_seq}"
                corr_barcode_id = f"{barcode_tag}:Z:{corrected_barcode_seq}"

                # Create new name with barcode information.
                if mapper == "ema":
                    # The EMA aligner requires reads in 10x format e.g.
                    # @READNAME:AAAAAAAATATCTACGCTCA BX:Z:AAAAAAAATATCTACGCTCA
                    new_name = ":".join([name_and_pos, corrected_barcode_seq])
                    new_name = " ".join((new_name, corr_barcode_id))
                    read1.name, read2.name = new_name, new_name
                elif args.mapper != "lariat":
                    new_name = "_".join([name_and_pos, raw_barcode_id, corr_barcode_id])
                    read1.name = " ".join([new_name, nr_and_index1])
                    read2.name = " ".join([new_name, nr_and_index2])
            else:
                summary["Reads missing barcode"] += 1

                # EMA and lairat aligner cannot handle reads without barcodes so these are skipped.
                if mapper in ["ema", "lariat"]:
                    continue

            # Write to out
            if mapper == "ema":
                output_chunk.append(
                    chunk_sep.join((
                        str(heap[corrected_barcode_seq]),
                        read1.name,
                        read1.sequence,
                        read1.qualities,
                        read2.name,
                        read2.sequence,
                        read2.qualities,
                        "\n"
                    ))
                )
                if not summary["Read pairs read"] % chunk_size:
                    output_chunk = sorted(output_chunk, key=lambda x: int(x.split(chunk_sep)[0]))
                    tmp_writer.writelines(output_chunk)
                    output_chunk.clear()
                    tmp_writer.close()
                    chunk_id += 1
                    tmp_writer = open(tmpdir / chunk_file_template.replace("*", str(chunk_id)), "w")

            elif mapper == "lariat":
                corrected_barcode_qual = "K" * len(corrected_barcode_seq)
                sample_index_qual = "K" * len(sample_index)
                print(
                    f"@{name_and_pos}",
                    read1.sequence,
                    read1.qualities,
                    read2.sequence,
                    read2.qualities,
                    f"{corrected_barcode_seq}-1",
                    corrected_barcode_qual,
                    sample_index,
                    sample_index_qual,
                    sep="\n", file=writer
                )
            else:
                summary["Read pairs written"] += 1
                writer.write(read1, read2)

        # Empty final cache
        if mapper == "ema":
            if output_chunk:
                logger.info(f"Writing chunk {chunk_id}, final")
                output_chunk = sorted(output_chunk, key=lambda x: int(x.split(chunk_sep)[0]))
                tmp_writer.writelines(output_chunk)

            tmp_writer.close()

            with ExitStack() as chunkstack:
                logger.info("Opening chunks for merge")
                chunks = []
                for chunk in tmpdir.iterdir():
                    chunks.append(chunkstack.enter_context(open(chunk)))

                logger.info("Merging chunks")
                for entry in heapq.merge(*chunks, key=lambda x: int(x.split(chunk_sep)[0])):
                    entry = entry.split(chunk_sep)
                    r1 = dnaio.Sequence(*entry[1:4])
                    r2 = dnaio.Sequence(*entry[4:7])
                    writer.write(r1, r2)
                    summary["Read pairs written"] += 1

    print_stats(summary, __name__)

    logger.info("Finished")


def parse_corrected_barcodes(open_file, summary, mapper, skip_singles=False):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing to a
    corrected canonical sequence
    :param open_file: starcode tabular output file.
    :param summary: collections.Counter object
    :param skip_singles: bool. Skip clusters of size 1
    :return: dict: raw sequences pointing to a corrected canonical sequence.
    """
    corrected_barcodes = dict()
    canonical_seqs = list()
    heap_index = {}
    for index, cluster in tqdm(enumerate(open_file), desc="Clusters processed"):
        canonical_seq, size, cluster_seqs = cluster.strip().split("\t", maxsplit=3)
        summary["Corrected barcodes"] += 1
        summary["Reads with corrected barcodes"] += int(size)
        summary["Uncorrected barcodes"] += len(cluster_seqs.split(","))

        if skip_singles and int(size) == 1:
            summary["Clusters size 1 skipped"] += 1
            continue

        corrected_barcodes.update({raw_seq: canonical_seq for raw_seq in cluster_seqs.split(",")})
        canonical_seqs.append(canonical_seq)

    if mapper == "ema":
        logger.info("Creating heap index for sorting barcodes for ema mapping.")
        # Scramble seqs to ensure no seqs sharing 16-bp prefix are neighbours.
        scramble(canonical_seqs)
        heap_index = {seq: nr for nr, seq in enumerate(canonical_seqs)}

    return corrected_barcodes, heap_index


def scramble(seqs, maxiter=10):
    """Scramble sequences by moving pairs with similar prefix appart. Loosely based on bubble sort"""
    swapped = True
    iteration = 1
    while swapped and iteration < maxiter:
        iteration += 1
        swapped = False
        for i in range(len(seqs) - 2):
            # Compare 16 first bases of both seqs
            if seqs[i][:16] == seqs[i + 1][:16]:
                # Swap the next-comming elements
                seqs[i + 2], seqs[i + 1] = seqs[i + 1], seqs[i + 2]
                # Set the flag to True so we'll loop again
                swapped = True
        logger.info("Iteration " + str(iteration))


class BarcodeReader:
    def __init__(self, filename):
        self._cache = dict()
        self._file = dnaio.open(filename, mode="r")
        self.barcodes = self.parse()

    def parse(self):
        for barcode in tqdm(self._file, desc="Uncorrected barcodes processed", disable=True):
            read_name, _ = barcode.name.split(maxsplit=1)
            yield read_name, barcode.sequence

    def get_barcode(self, read_name, maxiter=10):
        if read_name in self._cache:
            return self._cache.pop(read_name)

        for barcode_read_name, barcode_sequence in islice(self.barcodes, maxiter):
            # If read_name in next pair then parser lines are synced --> drop cache.
            if read_name == barcode_read_name:
                self._cache.clear()
                return barcode_sequence

            self._cache[barcode_read_name] = barcode_sequence
        return None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()


class Output:
    def __init__(self, file1, file2=None, interleaved=False, mapper=None):
        if mapper == "lariat":
            self._open_file = xopen(file1, mode='w')
        else:
            if interleaved:
                self._open_file = dnaio.open(file1, interleaved=True, mode="w", fileformat="fastq")
            else:
                self._open_file = dnaio.open(file1, file2=file2, mode="w", fileformat="fastq")

    def __enter__(self):
        return self._open_file

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._open_file.close()

    def write(self, *args):
        return self._open_file.write(*args)


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes.")
    parser.add_argument(
        "corrected_barcodes",
        help="FASTQ/FASTA for error corrected barcodes. Currently accepts output from starcode "
             "clustering with '--print-clusters' enabled.")
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
        "-s", "--sequence-tag", default="RX",
        help="SAM tag for storing the uncorrected barcode sequence. Default: %(default)s")
    parser.add_argument(
        "-m", "--mapper",
        help="Specify read mapper for labeling reads with barcodes. "
    )
    parser.add_argument(
        "--skip-singles", default=False, action="store_true",
        help="Skip adding barcode information for corrected barcodes with only one supporting read pair"
    )
