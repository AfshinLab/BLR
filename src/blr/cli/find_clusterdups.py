"""
Find barcode sequences originating from the same droplet/compartment.

Condition to call cluster duplicate:

    Two barcodes are considered cluster duplicates if they share two identical read pair positions
    (start, stop, direction) within a defined window (-w/--window). Positions are also required to not overlap in any
    way that is not allowed by the library technology. Specifically for BLR and stLFR libraries that use the Tn5
    tagmentase overlaps of 9 +/- 1 bp are allowed. For TELL-seq libraries that use the MuA tagmentase overlaps of
    5 +/- 1 bp are allowed. For 10x Chromium libraries any overlaps are allowed.
"""

from pysam import AlignmentFile, AlignedSegment, set_verbosity
from argparse import ArgumentError
import logging
from collections import Counter, deque, OrderedDict, defaultdict
import pickle
import numpy as np
from math import ceil

from blr.utils import get_bamtag, print_stats, tqdm, ACCEPTED_LIBRARY_TYPES

logger = logging.getLogger(__name__)


def main(args):
    if not (args.output_pickle or args.output_merges):
        raise ArgumentError(None, "Arguments --output-merges and/or --output-pickle required")

    run_find_clusterdups(
        input=args.input,
        output_pickle=args.output_pickle,
        output_merges=args.output_merges,
        barcode_tag=args.barcode_tag,
        buffer_size=args.buffer_size,
        window=args.window,
        min_mapq=args.min_mapq,
        quantile_threshold=args.quantile_threshold,
        library_type=args.library_type
    )


def run_find_clusterdups(
    input: str,
    output_pickle: str,
    output_merges: str,
    barcode_tag: str,
    buffer_size: int,
    window: int,
    min_mapq: int,
    quantile_threshold: float,
    library_type: str
):
    non_acceptable_overlap = get_non_acceptable_overlap_func(library_type)
    logger.info("Starting Analysis")
    summary = Counter()
    positions = OrderedDict()
    dup_positions = OrderedDict()
    chrom_prev = None
    pos_prev = 0
    uf = UnionFind()
    for read, mate in tqdm(paired_reads(input, min_mapq, summary), desc="Reading pairs"):
        barcode = get_bamtag(read, barcode_tag)
        if not barcode:
            summary["Non tagged reads"] += 2
            continue

        orientation = "F" if mate.is_read1 and read.is_read2 else "R"

        summary["Reads analyzed"] += 2

        chrom_new = read.reference_name
        pos_new = read.reference_start

        # Store position (5'-ends of R1 and R2) and orientation ('F' or 'R') which is used to group
        # duplicates. Based on picard MarkDuplicates definition, see
        # https://sourceforge.net/p/samtools/mailman/message/25062576/
        current_position = (mate.reference_start, read.reference_end, orientation)

        if abs(pos_new - pos_prev) > buffer_size or chrom_new != chrom_prev:
            find_duplicate_positions(positions, dup_positions)

            if chrom_new != chrom_prev:
                if chrom_prev is not None:
                    threshold = get_barcode_threshold(dup_positions, quantile=quantile_threshold)

                    summary["Barcode duplicate positions"] += len(dup_positions)

                    logger.info(f"Removing duplicate positions with greater than or equal to {threshold} barcodes "
                                f"for {chrom_prev}")
                    query_barcode_duplicates(dup_positions, uf, threshold, window,
                                             non_acceptable_overlap, summary)
                    positions.clear()
                    dup_positions.clear()
                chrom_prev = chrom_new

            pos_prev = pos_new

        positions.setdefault(current_position, PositionTracker(current_position)).add_barcode(barcode)

    # Process last chunk
    find_duplicate_positions(positions, dup_positions)
    threshold = get_barcode_threshold(dup_positions, quantile=quantile_threshold)

    summary["Barcode duplicate positions"] += len(dup_positions)

    logger.info(f"Using threshold {threshold} to filter duplicate position on contig {chrom_prev}")
    query_barcode_duplicates(dup_positions, uf, threshold, window,
                             non_acceptable_overlap, summary)

    summary["Duplicate compartments"] = len({v for _, v in uf.items()})
    summary["Barcodes removed"] = sum(1 for _ in uf.items()) - summary["Duplicate compartments"]

    # Write outputs
    if output_pickle:
        logger.info(f"Writing pickle object to {output_pickle}")
        with open(output_pickle, 'wb') as file:
            pickle.dump(uf.parents, file, pickle.HIGHEST_PROTOCOL)

    if output_merges:
        logger.info(f"Writing merges to {output_merges}")
        with open(output_merges, 'w') as file:
            for old_barcode, new_barcode in uf.items():
                if old_barcode != new_barcode:
                    print(old_barcode, new_barcode, sep=",", file=file)

    logger.info("Finished")
    print_stats(summary, name=__name__)


def paired_reads(path: str, min_mapq: int, summary):
    """
    Yield (forward_read, reverse_read) pairs for all properly paired read pairs in the input file.

    :param path: str, path to SAM file
    :param min_mapq: int
    :param summary: dict
    :return: read, mate: both as pysam AlignedSegment objects.
    """
    cache = dict()
    save = set_verbosity(0)  # Fix for https://github.com/pysam-developers/pysam/issues/939
    with AlignmentFile(path) as openin:
        for read in openin:
            summary["Total reads"] += 1
            # Requirements: read mapped, mate mapped and read has barcode tag
            # Cache read if matches requirements, continue with pair.
            if read.query_name in cache:
                mate = cache.pop(read.query_name)
            else:
                if pair_is_mapped_and_proper(read, min_mapq, summary):
                    cache[read.query_name] = read
                continue
            if pair_orientation_is_fr(read, mate, summary):
                yield read, mate
    set_verbosity(save)


def pair_is_mapped_and_proper(read: AlignedSegment, min_mapq: int, summary) -> bool:
    """
    Checks so read pair meets requirements before being used in analysis.
    :param read: pysam read
    :param min_mapq: int
    :param summary: dict
    :return: bool
    """
    if read.is_unmapped:
        summary["Unmapped reads"] += 1
        return False

    if read.mate_is_unmapped:
        summary["Unmapped reads"] += 1
        return False

    if read.mapping_quality < min_mapq:
        summary["Reads low MAPQ"] += 1
        return False

    if not read.is_proper_pair:
        summary["Reads not proper pair"] += 1
        return False
    return True


def pair_orientation_is_fr(read: AlignedSegment, mate: AlignedSegment, summary) -> bool:
    # Proper layout of read pair.
    # PAIR      |       mate            read
    # ALIGNMENTS|    ---------->      <--------
    # CHROMOSOME| ==================================>
    if not mate.is_reverse and read.is_reverse:
        return True
    summary["Reads with wrong orientation"] += 2
    return False


def find_duplicate_positions(positions, dup_positions):
    """
    Parse positions to find duplicate positions, i.e. containing more than one barcode, and add these to the
    dup_positions list. Non duplicate positions are removed.
    """
    positions_to_remove = list()
    for position, tracked_position in positions.items():
        if tracked_position.has_updated_barcodes:
            if tracked_position.is_duplicate():
                dup_positions[position] = tracked_position
            tracked_position.has_updated_barcodes = False
        else:
            positions_to_remove.append(position)

    for position in positions_to_remove:
        del positions[position]


def get_barcode_threshold(dup_positions, quantile: float = 0.99, min_threshold=6):
    """
    Calculate upper threshold for number of barcodes per position to include in duplicate query.
    """
    barcode_coverage = np.array([len(position.barcodes) for position in dup_positions.values()])
    if barcode_coverage.size > 0:
        return max(min_threshold, ceil(np.quantile(barcode_coverage, quantile)))
    else:
        return min_threshold


def query_barcode_duplicates(dup_positions, uf, threshold: float, window: int, non_acceptable_overlap,
                             summary):
    """
    Query barcode duplicates from list of duplicate positions. Position are filtered using the set threshold.
    """
    buffer_dup_pos = deque()
    for position, tracked_position in tqdm(dup_positions.items(), desc="Seeding dups"):
        if len(tracked_position.barcodes) < threshold:
            summary["Filtered barcode duplicate positions"] += 1
            seed_duplicates(
                uf=uf,
                buffer_dup_pos=buffer_dup_pos,
                position=tracked_position.position,
                position_barcodes=tracked_position.barcodes,
                window=window,
                non_acceptable_overlap=non_acceptable_overlap,
                summary=summary
            )


def get_non_acceptable_overlap_func(library_type: str):
    if library_type in {"blr", "stlfr"}:  # Tn5-type tagmentation
        return lambda x: x < 0 and x not in {-8, -9, -10}
    elif library_type in {"tellseq"}:  # MuA-type tagmentation
        return lambda x: x < 0 and x not in {-4, -5, -6}
    else:
        return lambda x: False


class PositionTracker:
    """
    Stores barcodes related to a position. The position is considered duplicate if more than one barcode is present.
    """
    def __init__(self, position):
        self.position = position
        self.reads = 0
        self.barcodes = set()
        self.has_updated_barcodes = False

    def add_barcode(self, barcode: str):
        self.has_updated_barcodes = barcode not in self.barcodes
        self.reads += 2
        self.barcodes.add(barcode)

    def is_duplicate(self) -> bool:
        return len(self.barcodes) > 1


def seed_duplicates(uf, buffer_dup_pos, position, position_barcodes, window: int, non_acceptable_overlap,
                    summary):
    """
    Identifies connected barcodes i.e. barcodes sharing two duplicate positions within the current window which is used
    to construct a graph. For Tn5 type libraries, overlapping positions are not compared unless they are allowed by Tn5
    tagmentation (i.e overlap 9±1 bp).
    :param uf: dict: Tracks which barcodes should be merged.
    :param buffer_dup_pos: list: Tracks previous duplicate positions and their barcode sets.
    :param position: tuple: Positions (start, stop) to be analyzed and subsequently saved to buffer.
    :param position_barcodes: seq: Barcodes at analyzed position
    :param window: int: Max distance allowed between postions to call barcode duplicate.
    """
    # Loop over list to get the positions closest to the analyzed position first. When positions
    # are out of the window size of the remaining buffer is removed.
    if not uf.same_component(*position_barcodes):
        for index, (compared_position, compared_barcodes) in enumerate(buffer_dup_pos):
            distance = position[0] - compared_position[1]

            if non_acceptable_overlap(distance):
                continue

            if distance <= window:
                barcode_intersect = position_barcodes & compared_barcodes

                # If two or more unique barcodes are found, update merge dict
                if len(barcode_intersect) >= 2:
                    uf.union(*barcode_intersect)
            else:
                # Remove positions outside of window (at end of list) since positions are sorted.
                for i in range(len(buffer_dup_pos) - index):
                    buffer_dup_pos.pop()
                break

        # Add new position at the start of the list.
        buffer_dup_pos.appendleft((position, position_barcodes))


class UnionFind:
    """Union-find data structure.
    Each UnionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:
    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.
    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.

    Based on Josiah Carlson's code,
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
    with significant additional changes by D. Eppstein.
    https://www.ics.uci.edu/~eppstein/PADS/UnionFind.py
    """

    def __init__(self, mapping=None):
        """Create a new  union-find structure."""
        self.parents = mapping if isinstance(mapping, dict) else {}

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __contains__(self, item):
        return item in self.parents

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def items(self):
        """Iterate over tuples of items and their root"""
        for x in self:
            yield x, self[x]

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]

        # Use lexicographical ordering to set main root in set
        heaviest = sorted(roots)[0]
        for r in roots:
            if r != heaviest:
                self.parents[r] = heaviest

    def connected_components(self):
        """Iterator for sets"""
        components = defaultdict(list)
        for item, root in self.items():
            components[root].append(item)

        for component in components.values():
            yield component

    def same_component(self, *objects) -> bool:
        """Returns true if all objects are present in the same set"""
        if all(x in self for x in objects):
            return len({self[x] for x in objects}) == 1
        return False

    def update(self, other):
        """Update sets based on other UnionFind instance"""
        for x, root in other.items():
            self.union(x, root)

    @classmethod
    def from_dict(cls, mapping):
        return cls(mapping.copy())


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="Coordinate-sorted SAM/BAM file tagged with barcodes.")
    parser.add_argument(
        "--output-pickle",
        help="Output python dict of barcodes to merge as pickle object.")
    parser.add_argument(
        "--output-merges",
        help="Output a CSV log file containing all merges done. File is in format: {old barcode id},{new barcode id}")
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument(
        "-w", "--window", type=int, default=30000,
        help="Window size. Duplicate positions within this distance will be used to find cluster "
        "duplicates. Default: %(default)s")
    parser.add_argument(
        "--buffer-size", type=int, default=200,
        help="Buffer size for collecting duplicates. Should be around read length. "
        "Default: %(default)s")
    parser.add_argument(
        "--min-mapq", type=int, default=0,
        help="Minimum mapping-quality to include reads in analysis Default: %(default)s")
    parser.add_argument(
        "-q", "--quantile-threshold", type=float, default=0.99,
        help="Quantile to filter out positions with to high barcode coverage. Default: %(default)s"
    )
    parser.add_argument(
        "-l", "--library-type", default="blr", choices=ACCEPTED_LIBRARY_TYPES,
        help="Library type of data. Default: %(default)s"
    )
