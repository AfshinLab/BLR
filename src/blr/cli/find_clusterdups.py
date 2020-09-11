"""
Find barcode cluster duplicates (two different barcode sequences origin to the same droplet, tagging the
same tagmented long molecule).

Condition to call barcode duplicate:

Two positions (positions defined as a unique set of read_start, read_stop, mate_start, mate_stop))
at a maximum of W (--window, default 100kbp, between = max(downstream_pos)-max(downstream_pos)) bp
apart sharing more than one barcode (share = union(bc_set_pos1, bc_set_pos2)).
"""

from pysam import AlignmentFile, AlignedSegment
from argparse import ArgumentError
import logging
from collections import Counter, deque, OrderedDict, defaultdict
import pickle
from copy import deepcopy
import numpy as np
from math import ceil

from blr.utils import get_bamtag, print_stats, tqdm

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
    tn5 = library_type in {'blr', 'stlfr'}
    logger.info("Starting Analysis")
    summary = Counter()
    positions = OrderedDict()
    dup_positions = OrderedDict()
    chrom_prev = None
    pos_prev = 0
    barcode_graph = BarcodeGraph()
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

                    logger.info(f"Removing duplicate positions with greater than or equal to {threshold} barcodes for {chrom_prev}")
                    query_barcode_duplicates(dup_positions, barcode_graph, threshold, window, tn5, summary)
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
    query_barcode_duplicates(dup_positions, barcode_graph, threshold, window, tn5, summary)

    # Write outputs
    if output_pickle:
        logger.info(f"Writing pickle object to {output_pickle}")
        with open(output_pickle, 'wb') as file:
            pickle.dump(barcode_graph.graph, file, pickle.HIGHEST_PROTOCOL)

    if output_merges:
        logger.info(f"Writing merges to {output_merges}")
        with open(output_merges, 'w') as file:
            for old_barcode, new_barcode in barcode_graph.iter_merges():
                summary["Barcodes removed"] += 1
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


def get_barcode_threshold(dup_positions, quantile: float = 0.999, min_threshold=6):
    """
    Calculate upper threshold for number of barcodes per position to include in duplicate query.
    """
    barcode_coverage = np.array([len(position.barcodes) for position in dup_positions.values()])
    return max(min_threshold, ceil(np.quantile(barcode_coverage, quantile)))


def query_barcode_duplicates(dup_positions, barcode_graph, threshold: float, window: int, tn5: bool, summary):
    """
    Query barcode duplicates from list of duplicate positions. Position are filtered using the set threshold.
    """
    buffer_dup_pos = deque()
    for position, tracked_position in tqdm(dup_positions.items(), desc="Seeding dups"):
        if len(tracked_position.barcodes) < threshold:
            summary["Filtered barcode duplicate positions"] += 1
            seed_duplicates(
                barcode_graph=barcode_graph,
                buffer_dup_pos=buffer_dup_pos,
                position=tracked_position.position,
                position_barcodes=tracked_position.barcodes,
                window=window,
                tn5=tn5,
                summary=summary
            )


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


def seed_duplicates(barcode_graph, buffer_dup_pos, position, position_barcodes, window: int, tn5: bool, summary):
    """
    Identifies connected barcodes i.e. barcodes sharing two duplicate positions within the current window which is used
    to construct a graph. For Tn5 type libraries, overlapping positions are not compared unless they are allowed by Tn5
    tagmentation (i.e overlap 9±1 bp).
    :param barcode_graph: dict: Tracks which barcodes should be merged.
    :param buffer_dup_pos: list: Tracks previous duplicate positions and their barcode sets.
    :param position: tuple: Positions (start, stop) to be analyzed and subsequently saved to buffer.
    :param position_barcodes: seq: Barcodes at analyzed position
    :param window: int: Max distance allowed between postions to call barcode duplicate.
    :param tn5: bool: Libray constructed using Tn5 tagmentation
    """
    # Loop over list to get the positions closest to the analyzed position first. When positions
    # are out of the window size of the remaining buffer is removed.
    for index, (compared_position, compared_barcodes) in enumerate(buffer_dup_pos):
        distance = position[0] - compared_position[1]

        # Skip comparison against overlapping reads unless they are from Tn5 tagmentation for Tn5-type libraries.
        if tn5 and distance < 0 and distance not in {-8, -9, -10}:
            continue

        if distance <= window:
            barcode_intersect = position_barcodes & compared_barcodes

            # If two or more unique barcodes are found, update merge dict
            if len(barcode_intersect) >= 2:
                barcode_graph.add_connected_barcodes(barcode_intersect)
        else:
            # Remove positions outside of window (at end of list) since positions are sorted.
            for i in range(len(buffer_dup_pos) - index):
                buffer_dup_pos.pop()
            break

    # Add new position at the start of the list.
    buffer_dup_pos.appendleft((position, position_barcodes))


class BarcodeGraph:
    def __init__(self, graph=None):
        self.graph = graph if graph else defaultdict(set)
        self._seen = set()

    def add_connected_barcodes(self, barcodes):
        for barcode in barcodes:
            self.graph[barcode].update(barcodes - set(barcode))

    def _update_component(self, nodes, component):
        """Breadth-first search of nodes"""
        new_nodes = set()
        for n in nodes:
            if n not in self._seen:
                self._seen.add(n)
                component.add(n)
                new_nodes |= self.graph[n]

        new_nodes = new_nodes.difference(component)

        if new_nodes:
            self._update_component(new_nodes, component)

    def components(self):
        """Generate all connected components of graph"""
        self._seen.clear()
        for node, neigbours in self.graph.items():
            if node not in self._seen:
                self._seen.add(node)
                component = {node}
                self._update_component(neigbours, component)
                yield component

    def get_merges(self):
        """Get dict of barcodes to merge in style of current_barcode -> new_barcode"""
        return {current: new for current, new in self.iter_merges()}

    def iter_merges(self):
        """Iterate over merges. Yields tuple with current_barcode and new_barcode"""
        for component in self.components():
            barcodes_sorted = sorted(component)
            barcode_min = barcodes_sorted[0]
            for barcode in barcodes_sorted[1:]:
                yield barcode, barcode_min

    def merge(self, other):
        """Merge BarcodeGraph objects"""
        for barcode, connected_barcodes in other.graph.items():
            self.graph[barcode] |= connected_barcodes

    @classmethod
    def from_graph(cls, graph):
        """Construct BarcodeGraph instance from existing graph"""
        return cls(deepcopy(graph))

    @classmethod
    def from_merges(cls, merges):
        """Construct BarcodeGraph instance from existing dict of merges"""
        graph = defaultdict(set)
        for current_barcode, new_barcode in merges.items():
            graph[current_barcode].add(new_barcode)
            graph[new_barcode].add(current_barcode)
        return cls(graph)


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
        "-q", "--quantile-threshold", type=float, default=0.999,
        help="Quantile to filter out positions with to high barcode coverage. Default: %(default)s"
    )
    parser.add_argument(
        "-l", "--library-type", default="blr", choices=["blr", "10x", "stlfr"],
        help="Library type of data. Default: %(default)s"
    )
