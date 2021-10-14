"""
Tag FASTQ headers with barcodes.

Takes in separate raw barcode FASTQ with matching headers and CLSTR file with corrected barcodes. Header formatting
depends on which mapper is requested.

ABOUT:

First the raw barcodes FASTQ is parser to get a dictionary:

    raw_barcodes[<HEADER>] = <RAW_BARCODE>

Then the corrected barcodes CLSTR file from starcode are parsed to get a dictionary:

    corrected_barcodes[<RAW_BARCODE>] = <CORRECTED_BARCODE>

The for each read-pair in the input FASTQ(s) the corrected barcode is recovered and used to tag
the read by including it in the header.

    <HEADER> ==> <RAW_BARCODE> ==> <CORRECTED_BARCODE>
"""

from contextlib import ExitStack
from heapq import merge
from itertools import islice
import logging
from pathlib import Path
import sys
import tempfile

import dnaio
from xopen import xopen

from blr.utils import tqdm, Summary, ACCEPTED_READ_MAPPERS

logger = logging.getLogger(__name__)

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT"
}


def main(args):
    run_tagfastq(
        uncorrected_barcodes=args.uncorrected_barcodes,
        corrected_barcodes=args.corrected_barcodes,
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
        output_nobc1=args.output_nobc1,
        output_nobc2=args.output_nobc2,
        output_bins=args.output_bins,
        nr_bins=args.nr_bins,
        barcode_tag=args.barcode_tag,
        sequence_tag=args.sequence_tag,
        mapper=args.mapper,
        min_count=args.min_count,
        pattern_match=args.pattern_match,
        sample_number=args.sample_nr,
    )


def run_tagfastq(
        uncorrected_barcodes: str,
        corrected_barcodes: str,
        input1: str,
        input2: str,
        output1: str,
        output2: str,
        output_nobc1: str,
        output_nobc2: str,
        output_bins: str,
        nr_bins: int,
        barcode_tag: str,
        sequence_tag: str,
        mapper: str,
        min_count: int,
        pattern_match: str,
        sample_number: int,
):
    logger.info("Starting")
    summary = Summary()
    # Get the corrected barcodes and create a dictionary pointing each raw barcode to its
    # canonical sequence.
    template = [set(IUPAC[base]) for base in pattern_match] if pattern_match else []
    with xopen(corrected_barcodes) as reader:
        corrected_barcodes, heap = parse_corrected_barcodes(reader, summary, mapper, template,
                                                            min_count=min_count)

    in_interleaved = not input2
    logger.info(f"Input detected as {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    out_interleaved = not output2 or not output1
    if output_bins is not None:
        logger.info(f"Writing output as binned interleaved FASTQ to {output_bins}.")
        output_bins = Path(output_bins)
        output_bins.mkdir(exist_ok=True)
        output1, output2 = None, None
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
        writer = stack.enter_context(Output(file1=output1, file2=output2, interleaved=out_interleaved,
                                            file_nobc1=output_nobc1, file_nobc2=output_nobc2, mapper=mapper,
                                            bins_dir=output_bins))
        uncorrected_barcode_reader = stack.enter_context(BarcodeReader(uncorrected_barcodes))
        chunks = None
        if mapper in ["ema", "lariat"]:
            chunks = stack.enter_context(ChunkHandler(chunk_size=1_000_000))

        for read1, read2, corrected_barcode_seq in parse_reads(reader, corrected_barcodes, uncorrected_barcode_reader,
                                                               barcode_tag, sequence_tag, mapper):
            summary["Read pairs read"] += 1
            if corrected_barcode_seq is None:
                summary["Reads missing barcode"] += 1

                # Write non barcoded reads to separate file if exists for ema.
                if mapper == "ema" and output_nobc1 is not None:
                    summary["Read pairs written"] += 1
                    writer.write_nobc(read1, read2)
                    continue

                # EMA and lairat aligner cannot handle reads without barcodes so these are skipped.
                if mapper in ["ema", "lariat"]:
                    continue

            # Write to out
            if mapper == "ema":
                chunks.build_chunk(
                    f"{str(heap[corrected_barcode_seq])}\t"
                    f"{read1.name}\t"
                    f"{read1.sequence}\t"
                    f"{read1.qualities}\t"
                    f"{read2.sequence}\t"
                    f"{read2.qualities}\n"
                )
            elif mapper == "lariat":
                corrected_barcode_qual = "K" * len(corrected_barcode_seq)
                chunks.build_chunk(
                    f"{str(heap[corrected_barcode_seq])}\t"
                    f"@{read1.name}\t"
                    f"{read1.sequence}\t"
                    f"{read1.qualities}\t"
                    f"{read2.sequence}\t"
                    f"{read2.qualities}\t"
                    f"{corrected_barcode_seq}-{sample_number}\t"
                    f"{corrected_barcode_qual}\t"
                    "AAAAAA\t"
                    "KKKKKK\n"
                )
            else:
                summary["Read pairs written"] += 1
                writer.write(read1, read2)

        # Empty final cache
        if chunks is not None:
            chunks.write_chunk()

        if output_bins is not None:
            bin_size = (summary["Read pairs read"] - summary["Reads missing barcode"]) // nr_bins
            logger.info(f"Using bin of size {bin_size}.")
            writer.set_bin_size(bin_size)

        if mapper == "ema":
            write_ema_output(chunks, writer, summary)
        elif mapper == "lariat":
            write_lariat_output(chunks, writer, summary)

    summary.print_stats(__name__)

    logger.info("Finished")


def write_ema_output(chunks, writer, summary):
    for entry in chunks.parse_chunks():
        r1 = dnaio.Sequence(*entry[1:4])
        r2 = dnaio.Sequence(entry[1], *entry[4:6])
        writer.write(r1, r2, heap=entry[0])
        summary["Read pairs written"] += 1


def write_lariat_output(chunks, writer, summary):
    for entry in chunks.parse_chunks():
        lines = "\n".join(entry[1:]) + "\n"
        writer.write(lines, heap=entry[0])
        summary["Read pairs written"] += 1


def parse_reads(reader, corrected_barcodes, uncorrected_barcode_reader, barcode_tag, sequence_tag, mapper):
    for read1, read2 in tqdm(reader, desc="Read pairs processed"):
        # Header parsing
        # TODO Handle reads with single header
        name_and_pos, nr_and_index1 = read1.name.split(maxsplit=1)

        uncorrected_barcode_seq = uncorrected_barcode_reader.get_barcode(name_and_pos)
        corrected_barcode_seq = corrected_barcodes.get(uncorrected_barcode_seq, None)

        # Check if barcode was found and update header with barcode info.
        if corrected_barcode_seq is not None:
            corr_barcode_id = f"{barcode_tag}:Z:{corrected_barcode_seq}"

            # Create new name with barcode information.
            if mapper == "ema":
                # The EMA aligner requires reads in 10x format e.g.
                # @READNAME:AAAAAAAATATCTACGCTCA BX:Z:AAAAAAAATATCTACGCTCA
                read1.name = f"{name_and_pos}:{corrected_barcode_seq} {corr_barcode_id}"
                read2.name = read1.name
            elif mapper != "lariat":
                _, nr_and_index2 = read2.name.split(maxsplit=1)
                raw_barcode_id = f"{sequence_tag}:Z:{uncorrected_barcode_seq}"
                read1.name = f"{name_and_pos}_{raw_barcode_id}_{corr_barcode_id} {nr_and_index1}"
                read2.name = f"{name_and_pos}_{raw_barcode_id}_{corr_barcode_id} {nr_and_index2}"

        yield read1, read2, corrected_barcode_seq


def parse_corrected_barcodes(open_file, summary, mapper, template, min_count=0):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing to a
    corrected canonical sequence
    :param open_file: starcode tabular output file.
    :param summary: collections.Counter object
    :param min_count: int. Skip clusters with fewer than min_count reads.
    :return: dict: raw sequences pointing to a corrected canonical sequence.
    """
    corrected_barcodes = {}
    canonical_seqs = []
    heap_index = {}
    for canonical_seq, size, cluster_seqs in tqdm(parse_clstr(open_file), desc="Clusters processed"):
        summary["Corrected barcodes"] += 1
        summary["Reads with corrected barcodes"] += size
        summary["Uncorrected barcodes"] += len(cluster_seqs)

        if size < min_count:
            summary["Barcodes with too few reads"] += 1
            summary["Reads with barcodes with too few reads"] += size
            continue

        if template and not match_template(canonical_seq, template):
            summary["Barcodes not matching pattern"] += 1
            summary["Reads with barcodes not matching pattern"] += size
            continue

        corrected_barcodes.update({raw_seq: canonical_seq for raw_seq in cluster_seqs})
        canonical_seqs.append(canonical_seq)

    if mapper in ["ema", "lariat"]:
        logger.info("Creating heap index for sorting barcodes for ema mapping.")

        # Scramble seqs to ensure no seqs sharing 16-bp prefix are neighbours for ema.
        if mapper == "ema":
            scramble(canonical_seqs, maxiter=100)

        heap_index = {seq: nr for nr, seq in enumerate(canonical_seqs)}

    return corrected_barcodes, heap_index


def parse_clstr(clstr_file):
    """"Generator to parse open .clstr files from starcode."""
    for cluster in clstr_file:
        canonical_seq, size, cluster_seqs_list = cluster.strip().split("\t")
        yield canonical_seq, int(size), cluster_seqs_list.split(",")


def scramble(seqs, maxiter=10):
    """Scramble sequences by moving pairs with similar prefix appart. Loosely based on bubble sort"""
    swapped = True
    iteration = 1
    start_from = 0
    while swapped and iteration < maxiter:
        iteration += 1
        swapped = False
        swap_pos = []
        for i in range(start_from, len(seqs) - 2):
            # Compare 16 first bases of both seqs
            if seqs[i][:16] == seqs[i + 1][:16]:
                swap_pos.append(i)
                # Swap the next-comming elements
                seqs[i + 2], seqs[i + 1] = seqs[i + 1], seqs[i + 2]
                # Set the flag to True so we'll loop again
                swapped = True

        start_from = min(swap_pos) if swap_pos else 0

    if not swapped:
        logger.info("Scrambling done!")
    else:
        logger.warning("Scrambling reached maxiter")


def match_template(sequence: str, template) -> bool:
    if len(sequence) != len(template):
        return False

    for base, accepted_bases in zip(sequence, template):
        if base not in accepted_bases:
            return False
    return True


class BarcodeReader:
    def __init__(self, filename):
        self._cache = {}
        self._file = dnaio.open(filename, mode="r")
        self.barcodes = self.parse()

    def parse(self):
        for barcode_entry in self._file:
            name, *_ = barcode_entry.name.split(" ")
            yield name, barcode_entry.sequence

    def get_barcode(self, read_name, maxiter=128):
        if read_name in self._cache:
            return self._cache.pop(read_name)

        for barcode_read_name, barcode_sequence in islice(self.barcodes, maxiter):
            # If read_name in next pair then parser lines are synced --> drop cache.
            if read_name == barcode_read_name:
                self._cache.clear()
                return barcode_sequence

            self._cache[barcode_read_name] = barcode_sequence

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()


class Output:
    """
    Output handler for different output formats required by different mappers.
    """
    BIN_FASTQ_TEMPLATE = "ema-bin-*"  # Same name as in `ema preproc`.

    def __init__(self, file1=None, file2=None, interleaved=False, file_nobc1=None, file_nobc2=None, mapper=None,
                 bins_dir=None):
        self._mapper = mapper

        self._bin_nr = 0
        self._reads_written = 0
        self._bin_size = None
        self._bins_dir = bins_dir
        self._prev_heap = None
        self._bin_filled = False
        self._pre_write = lambda *args: None
        self._post_write = lambda *args: None

        self._open_file = None
        if file1 is not None:
            self._open_file = self._setup_single_output(file1, file2, interleaved)
        elif bins_dir is not None:
            self._open_new_bin()
            self._pre_write = self._open_new_bin_if_full
            self._post_write = self._check_bin_full
        else:
            sys.exit("Either file1 or bins_dir need to be provided.")

        self._open_file_nobc = None
        if file_nobc1 is not None:
            self._open_file_nobc = self._setup_single_output(file_nobc1, file_nobc2,
                                                             interleaved=file_nobc2 is None)
        # Setup writers based on mapper.
        if self._mapper == "lariat":
            self._write = self._write_lariat
            self._write_nobc = self._write_nobc_lariat
        else:
            self._write = self._write_default
            self._write_nobc = self._write_nobc_default

    def set_bin_size(self, value):
        self._bin_size = value

    def _setup_single_output(self, file1, file2, interleaved):
        if self._mapper == "lariat":
            if file2 is not None:
                Path(file2).touch()
            return xopen(file1, mode='w')
        else:
            return dnaio.open(file1, file2=file2, interleaved=interleaved, mode="w", fileformat="fastq")

    def _open_new_bin(self):
        if self._open_file is not None:
            self._open_file.close()
        bin_nr_str = str(self._bin_nr).zfill(3)
        file_name = self._bins_dir / Output.BIN_FASTQ_TEMPLATE.replace("*", bin_nr_str)
        self._open_file = dnaio.open(file_name, interleaved=True, mode="w", fileformat="fastq")
        self._bin_nr += 1

    def _open_new_bin_if_full(self, heap):
        # Start a new bin if the current is full while not splitting heaps over separate bins
        if self._bin_filled and heap != self._prev_heap:
            self._bin_filled = False
            self._open_new_bin()
            logger.debug(f"Bin overflow = {self._reads_written} ")

        self._prev_heap = heap

    def _check_bin_full(self):
        if self._reads_written > self._bin_size:
            self._bin_filled = True
            self._reads_written = 0

    def _write_default(self, read1, read2):
        self._open_file.write(read1, read2)

    def _write_lariat(self, read1, _):
        self._open_file.write(read1)

    def write(self, read1, read2=None, heap=None):
        self._pre_write(heap)

        # Write reads to output file
        self._write(read1, read2)
        self._reads_written += 1

        self._post_write()

    def _write_nobc_default(self, read1, read2):
        self._open_file_nobc.write(read1, read2)

    def _write_nobc_lariat(self, read1, _):
        self._open_file_nobc.write(read1)

    def write_nobc(self, read1, read2=None):
        if self._open_file_nobc is not None:
            self._write_nobc(read1, read2)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._open_file.close()
        if self._open_file_nobc is not None:
            self._open_file_nobc.close()


class ChunkHandler:
    def __init__(self, chunk_size: int = 100_000):
        # Required for ema sorted output
        # Inspired by: https://stackoverflow.com/questions/56948292/python-sort-a-large-list-that-doesnt-fit-in-memory
        self._output_chunk = []
        self._chunk_size = chunk_size
        self._chunk_id = 0
        self._tmpdir = Path(tempfile.mkdtemp(prefix="tagfastq_sort"))
        self._chunk_file_template = "chunk_*.tsv"
        self._chunk_files = []
        self._chunk_sep = "\t"
        self._tmp_writer = self.create_writer()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._tmp_writer.close()

    def create_writer(self):
        if hasattr(self, "_tmp_writer"):
            self._tmp_writer.close()
        tmpfile = self._tmpdir / self._chunk_file_template.replace("*", str(self._chunk_id))
        self._chunk_files.append(tmpfile)
        self._chunk_id += 1
        return open(tmpfile, "w")

    def build_chunk(self, line: str):
        """Add entry to write to temporary file, first argument should be the heap index"""
        self._output_chunk.append(line)

        if len(self._output_chunk) > self._chunk_size:
            self.write_chunk()
            self._tmp_writer = self.create_writer()

    def write_chunk(self):
        self._output_chunk.sort(key=self._get_heap)
        self._tmp_writer.writelines(self._output_chunk)
        self._output_chunk *= 0  # Clear list faster than list.clear(). See https://stackoverflow.com/a/44349418

    def parse_chunks(self):
        if not self._tmp_writer.closed:
            self._tmp_writer.close()

        with ExitStack() as chunkstack:
            logger.info("Opening chunks for merge")
            chunks = [chunkstack.enter_context(chunk.open(mode="r")) for chunk in self._chunk_files]
            logger.info("Merging chunks")
            for entry in merge(*chunks, key=self._get_heap):
                yield entry.strip().split(self._chunk_sep)

    def _get_heap(self, x):
        return int(x.split(self._chunk_sep)[0])


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes."
    )
    parser.add_argument(
        "corrected_barcodes",
        help="FASTQ/FASTA for error corrected barcodes. Currently accepts output from starcode "
             "clustering with '--print-clusters' enabled."
    )
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
             "result is written as interleaved."
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
        help=f"Output interleaved FASTQ split into bins named '{Output.BIN_FASTQ_TEMPLATE}' in the provided "
             f"directory. Entries will be grouped based on barcode. Only used for ema mapping."
    )
    parser.add_argument(
        "--nr-bins", type=int, default=100,
        help="Number of bins to split reads into when using the '--output-bins' alternative. Default: %(default)s."
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s."
    )
    parser.add_argument(
        "-s", "--sequence-tag", default="RX",
        help="SAM tag for storing the uncorrected barcode sequence. Default: %(default)s."
    )
    parser.add_argument(
        "-m", "--mapper", default="bowtie2", choices=ACCEPTED_READ_MAPPERS,
        help="Specify read mapper for labeling reads with barcodes. Selecting 'ema' or 'lariat' produces output "
             "required for these particular mappers. Default: %(default)s."
    )
    parser.add_argument(
        "-c", "--min-count", default=0, type=int,
        help="Minimum number of reads per barcode to tag read name. Default: %(default)s."
    )
    parser.add_argument(
        "-p", "--pattern-match",
        help="IUPAC barcode string to match against corrected barcodes e.g. for DBS it is usualy "
             "BDHVBDHVBDHVBDHVBDHV. Non-matched barcodes will be removed."
    )
    parser.add_argument(
        "--sample-nr", type=int, default=1,
        help="Sample number to append to barcode string. Default: %(default)s."
    )
