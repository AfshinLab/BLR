import re
from pathlib import Path
import sys
import pysam
from dataclasses import dataclass
import numpy as np
import os
from collections import namedtuple, Counter, defaultdict, OrderedDict
import contextlib

from blr import __version__

if sys.stderr.isatty():
    from tqdm import tqdm
else:
    def tqdm(iterable, **kwargs):
        return iterable


ACCEPTED_LIBRARY_TYPES = ["dbs", "blr", "10x", "stlfr", "tellseq"]  # TODO Remove blr
ACCEPTED_READ_MAPPERS = ["ema", "lariat", "bwa", "bowtie2", "minimap2"]


def is_1_2(s, t):
    """
    Determine whether s and t are identical except for a single character of
    which one of them is '1' and the other is '2'.
    """
    differences = 0
    one_two = {"1", "2"}
    for c1, c2 in zip(s, t):
        if c1 != c2:
            differences += 1
            if differences == 2:
                return False
            if {c1, c2} != one_two:
                return False
    return differences == 1


def guess_paired_path(path: Path):
    """
    Given the path to a file that contains the sequences for the first read in a
    pair, return the file that contains the sequences for the second read in a
    pair. Both files must have identical names, except that the first must have
    a '1' in its name, and the second must have a '2' at the same position.

    Return None if no second file was found or if there are too many candidates.

    >>> guess_paired_path(Path('file.1.fastq.gz'))  # doctest: +SKIP
    'file.2.fastq.gz'  # if that file exists
    """
    name = path.name
    # All lone 1 digits replaced with '?'
    name_with_globs = re.sub(r"(?<![0-9])1(?![0-9])", "?", name)
    paths = [p for p in path.parent.glob(name_with_globs) if is_1_2(str(p), str(path))]
    if len(paths) == 1:
        return paths[0]
    return None


def get_bamtag(pysam_read: pysam.AlignedSegment, tag: str, default=None):
    """
    Fetches tags from bam files. Return default value of the same type if not found.
    """
    try:
        return pysam_read.get_tag(tag)
    except KeyError:
        return default


class Summary(Counter):

    def print_stats(self, name=None, value_width=15, print_to=sys.stderr):
        """
        Prints stats in nice table with two column for the key and value pairs in summary
        :param name: name of script for header e.g. '__name__'
        :param value_width: width for values column in table
        :param print_to: Where to direct output. Default: stderr
        """
        # Get widths for formatting
        max_name_width = max(map(len, self.keys()), default=10)
        width = value_width + max_name_width + 1

        # Header
        print("="*width, file=print_to)
        print(f"STATS SUMMARY - {name}", file=print_to)
        print("-"*width, file=print_to)

        # Print stats in columns
        for name, value in self.items():
            value_str = str(value)
            if isinstance(value, (int, np.integer)):
                value_str = f"{value:>{value_width},}"
            elif isinstance(value, (float, np.float)):
                value_str = f"{value:>{value_width+4},.3f}"

            print(f"{name:<{max_name_width}} {value_str}", file=print_to)
        print("="*width, file=print_to)


class PySAMIO:
    """ Reader and writer for BAM/SAM files that automatically attaches processing step information to header """

    def __init__(self, inname: str, outname: str, name: str, inmode: str = "rb", outmode: str = "wb"):
        """
        :param inname: Path to input SAM/BAM file.
        :param outname: Path to output SAM/BAM file.
        :param name: __name__ variable from script.
        :param inmode: Reading mode for input file. 'r' for SAM and 'rb' for BAM.
        :param outmode: Reading mode for output file. 'r' for SAM and 'rb' for BAM.
        """
        self._save = pysam.set_verbosity(0)  # Fix for https://github.com/pysam-developers/pysam/issues/939
        self.infile = pysam.AlignmentFile(inname, inmode)
        self.header = self._make_header(name)
        self.outfile = pysam.AlignmentFile(outname, outmode, header=self.header)

    def __enter__(self):
        return self.infile, self.outfile

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.infile.close(), self.outfile.close()
        pysam.set_verbosity(self._save)
        return isinstance(exc_val, OSError)

    def _make_header(self, name):
        """
        Create SAM header dict with new tool and command line argument information based on template file. Appends
        new PG entry with tool name (ID), software name (PN), command line arguments (CL) to track the tools applied
        to the file. Use in output SAM/BAM file as 'header' attribute.

        :param name: string. Pass '__name__' variable to be used to get program and tool name.
        :return: pysam.AlignmentHeader object
        """

        def make_unique(identifier, prev_entries):
            """
            Compare idenifier to other identifiers in header and make unique if neccessary by adding suffix.
            :param identifier: string. Program record identifier.
            :param prev_entries: dict. Dictionary of "PG" header entries from pysam.AlignmentFile header.
            :return: string. Updated identifier
            """
            nr = 0
            updated_identifier = identifier
            while any(updated_identifier == e["ID"] for e in prev_entries):
                nr += 1
                updated_identifier = ".".join([identifier, str(nr)])
            return updated_identifier

        # Process info strings for header
        id_name = name.split(".")[-1]
        program_name = name.split(".")[0]
        cmd_line = f"\"{' '.join(sys.argv)}\""

        header = self.infile.header.to_dict()
        pg_entries = header["PG"]
        # Make sure id_name is unique by adding numbers at end if needed.
        id_name = make_unique(id_name, pg_entries)

        pg_entries.append({
            "ID": id_name,       # Program record identifier. Must be unique
            "PN": program_name,  # Program name
            "CL": cmd_line,      # Command line arguments string.
            "VN": __version__    # Version information
        })
        header["PG"] = pg_entries

        return pysam.AlignmentHeader.from_dict(header)


def calculate_N50(lengths):
    """
    Calculate N50 metric for list of integers.
    Based on https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c.
    :param lengths: list containing integers.
    :return int. N50 metric
    """
    lengths = np.array(lengths)
    lengths[::-1].sort()

    csum = np.cumsum(lengths)
    n2 = int(sum(lengths) / 2)

    # get index for cumsum >= N/2
    csumn2 = min(csum[csum >= n2])
    ind = np.where(csum == csumn2)
    return lengths[ind[0][0]]


class ReadGroup:
    """
    Read group information for read tagging.
    See SAM format for details: https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    def __init__(self, identifier, library, sample, platfrom_unit, platform):
        self.identifier = identifier
        self.library = library
        self.sample = sample
        self.platform_unit = platfrom_unit
        self.platform = platform

        self.ID = f"ID:{self.identifier}"
        self.LB = f"LB:{self.library}"
        self.SM = f"SM:{self.sample}"
        self.PU = f"PU:{self.platform_unit}"
        self.PL = f"PL:{self.platform}"

    def __repr__(self):
        return r"\t".join(["@RG", self.ID, self.LB, self.SM, self.PU, self.PL])


@dataclass
class FastaIndexRecord:
    name: str
    length: int


def parse_fai(file):
    """Parse a FASTA index file (.fai) and return a list of FastaIndexRecords"""
    chromosomes = []
    for line in file:
        if not line.strip():
            continue
        fields = line.split("\t")
        name = fields[0]
        length = int(fields[1])
        chromosomes.append(FastaIndexRecord(name, length))
    return chromosomes


def chromosome_chunks(index_records, size=20_000_000):
    """
    Given a list of chromosomes (as FastaIndexRecords), split into chunks such that each
    chunk has at least one chromosome and then as many additional chromosomes as possible
    without exceeding the chunksize.
    """
    chunk = []
    chunk_size = 0
    for index_record in index_records:
        if chunk and chunk_size + index_record.length > size:
            yield chunk
            chunk = []
            chunk_size = 0
        chunk.append(index_record)
        chunk_size += index_record.length
    if chunk:
        yield chunk


def generate_chunks(reference, size=200_000_000, phasing_contigs_string=None, contigs_skipped=None):
    """
    Group contigs in reference into chunks
    """
    chunks = defaultdict(list)
    if reference is not None:
        try:
            with open(reference + ".fai") as f:
                reference_contigs = list(parse_fai(f))
        except FileNotFoundError:
            sys.exit(
                f"The genome index file {reference}.fai is missing. Please create it with 'samtools faidx'"
            )

        # Define primary contigs which go through every processing step (possibly excluding the phasing steps).
        primary_contigs = {c.name for c in reference_contigs}
        if contigs_skipped is not None:
            skip_pattern = re.compile(contigs_skipped)
            primary_contigs = {c.name for c in reference_contigs if not re.match(skip_pattern, c.name)}

        # Define phasing contigs as a subset of the primary contigs that also are phased.
        phasing_contigs = primary_contigs.copy()
        if phasing_contigs_string is not None:
            phasing_contigs = set(phasing_contigs_string.split(","))
            if not all(name in {c.name for c in reference_contigs} for name in phasing_contigs):
                sys.exit(f"The contigs in phasing_contigs ({', '.join(phasing_contigs)}) does not match those found in"
                         f"the reference {reference} ({', '.join(c.name for c in reference_contigs)}).")

            if len(phasing_contigs-primary_contigs) > 0:
                sys.exit("The phasing contigs should be a subset of the primary contigs. Currently the contig(s) "
                         f"{','.join(phasing_contigs-primary_contigs)} are missing from the primary contigs. Adjust "
                         "the 'contigs_skipped' pattern or contigs defined in 'phasing_contigs'.")

        primary_not_phased = primary_contigs - phasing_contigs

        # We want to make sure that contigs destined for different rules are not combined into the same chunks.
        # Therefore we first chunk the phasing set then chunk the remaining contigs that are part of the primary set.
        # Finally any remaining contigs are chunked.
        chunks["phased"] = list(chromosome_chunks(filter(lambda x: x.name in phasing_contigs, reference_contigs),
                                                  size=size))
        chunks["primary"] = chunks["phased"].copy()
        chunks["primary"] += list(chromosome_chunks(filter(lambda x: x.name in primary_not_phased, reference_contigs),
                                                    size=size))
        chunks["all"] = chunks["primary"].copy()
        chunks["all"] += list(chromosome_chunks(filter(lambda x: x.name not in primary_contigs, reference_contigs),
                                                size=size))

        chunks["not_phased"] = [chunk for chunk in chunks["all"] if chunk[0].name not in phasing_contigs]
        chunks["not_primary"] = [chunk for chunk in chunks["all"] if chunk[0].name not in primary_contigs]
        chunks["primary_not_phased"] = [chunk for chunk in chunks["all"] if chunk[0].name in primary_not_phased]

    return chunks


def symlink_relpath(source, target):
    """
    Generate a symlink to source that is relative to the target location. Corresponds roughtly to 'ln -rs'.
    """
    commonpath = os.path.commonpath([source, target])
    relpath_source = os.path.relpath(source, commonpath)
    os.symlink(relpath_source, target)


def parse_phaseblocks(file):
    """
    Format  description from https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

    Example entry.
    ```
    BLOCK: offset: 6 len: 4 phased: 2 SPAN: 23514 fragments 1
    ```

    offset: <SNV offset>
    len: <SNV span of block>
    phased: <# SNVs phased>
    SPAN: <base pair span of block>
    fragments <# of fragments in block>

    """
    Phaseblock = namedtuple("Phaseblock", ["snv_span", "phased_snvs", "length", "fragments"])
    for line in file:
        if line.startswith("BLOCK:"):
            contents = line.split()
            yield Phaseblock(int(contents[4]), int(contents[6]), int(contents[8]), int(contents[10]))
        else:
            continue


@contextlib.contextmanager
def smart_open(filename=None):
    # From https://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    if filename and filename is not None and filename != "-":
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


class LastUpdatedOrderedDict(OrderedDict):
    """Store items in the order the keys were last added"""
    # Taken from https://docs.python.org/3/library/collections.html#ordereddict-examples-and-recipes
    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.move_to_end(key)


def parse_filters(string: str):
    """Parse filter string with semicolon separate pairs of commad separated filter expressions and names for
    GATK VariantFiltration commands."""
    filters = ""
    for pair in string.split(";"):
        filter_expression, filter_name = pair.split(",")
        filters += f" -filter {filter_expression} --filter-name {filter_name}"
    return filters
