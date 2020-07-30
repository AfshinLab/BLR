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
from collections import Counter, OrderedDict, defaultdict
import pickle
from copy import deepcopy
from scipy import sparse, stats
import numpy as np

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
        mapq_threshold=args.mapping_quality,
    )


def run_find_clusterdups(
    input: str,
    output_pickle: str,
    output_merges: str,
    barcode_tag: str,
    buffer_size: int,
    window: int,
    mapq_threshold: int,
):
    logger.info("Starting Analysis")
    summary = Counter()
    positions = OrderedDict()
    chrom_prev = None
    pos_prev = 0
    barcode_graph = BarcodeGraph()
    dup_pos = defaultdict(set)
    for read, mate in tqdm(paired_reads(input, mapq_threshold, summary), desc="Reading pairs"):
        barcode = get_bamtag(read, barcode_tag)
        if not barcode:
            summary["Non tagged reads"] += 2
            continue

        if mate.is_read1 and read.is_read2:
            orientation = "F"
        else:
            orientation = "R"

        summary["Reads analyzed"] += 2

        chrom_new = read.reference_name
        pos_new = read.reference_start

        # Store position (5'-ends of R1 and R2) and orientation ('F' or 'R') with is used to group
        # duplicates. Based on picard MarkDuplicates definition, see
        # https://sourceforge.net/p/samtools/mailman/message/25062576/
        current_position = (read.reference_name, mate.reference_start, read.reference_end, orientation)

        if abs(pos_new - pos_prev) > buffer_size or chrom_new != chrom_prev:
            find_barcode_duplicates(positions, dup_pos)

            if chrom_new != chrom_prev:
                positions.clear()
                chrom_prev = chrom_new

            pos_prev = pos_new

        if current_position not in positions:
            positions[current_position] = PositionTracker(current_position)
        positions[current_position].add_barcode(barcode)

    # Process last chunk
    find_barcode_duplicates(positions, dup_pos)

    # Create index for each position
    pos_index = {pos: index for index, pos in enumerate(dup_pos)}

    # Assign each barcode to a set of position indexes.
    barcode_dup_pos = defaultdict(set)
    for pos, barcodes in dup_pos.items():
        for barcode in barcodes:
            barcode_dup_pos[barcode].add(pos_index[pos])

    # Generate coordinates for creating a matrix with barcodes as row and positions as columns
    # TODO Should we skip barcodes with to few positions here?
    rows = []
    cols = []
    index_barcode = dict()
    for row, (barcode, positions) in enumerate(barcode_dup_pos.items()):
        index_barcode[row] = barcode
        for col in positions:
            rows.append(row)
            cols.append(col)

    # -------------------------------------------------------------------------------------------------
    # The code below is inspired by the script coalescence.py from the 10x Genomics longranger software
    # https://github.com/10XGenomics/longranger
    # Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
    # -------------------------------------------------------------------------------------------------

    # Create sparse matrix
    nbarcodes = len(barcode_dup_pos)
    npostions = len(pos_index)
    matrix = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nbarcodes, npostions), dtype=int)

    # Get the position occupancy
    position_counts = np.array((matrix > 0).sum(axis=0)).flatten().astype('int')
    logger.info(f"Nr of duplicate positions: {len(position_counts)}")

    # switch off high-coverage positions
    percentile = 99.95
    high_cov_threshold = np.percentile(position_counts, percentile)
    logger.info(f"High coverage threshold set to: {high_cov_threshold} using percentile {percentile}")

    high_cov_positions = np.where(position_counts > high_cov_threshold)[0]
    logger.info(f"Nr of positions removed: {len(high_cov_positions)}")

    (r, c) = matrix.nonzero()
    for pos in high_cov_positions:
        matrix.data[c == pos] = 0
        position_counts[pos] = 0

    # TODO How does this work?
    # Adjust for 'effective genome size' based on the distribution over positions
    # i.e. more skewed distribution -> fewer effective position
    effective_positions_factor = ((position_counts / position_counts.sum()) ** 2).sum()
    effective_positions = 1/effective_positions_factor
    logger.info(f"Nr of effective positions: {effective_positions}")

    # number of positions on each bc
    num_positions = np.array((matrix > 0.0).sum(axis=1)).flatten()
    logger.info(f"mean num positions: {num_positions.mean()},  median: {np.median(num_positions)}")

    # TODO How does this work?
    # compute an initial threshold for coalescence
    # first, compute the expected number of overlaps between two BCs having the mean number of fragments
    expected_overlaps = (num_positions.mean() ** 2) / effective_positions

    # Barcode need to share more duplicate positions than this for overlaps to be considered
    min_threshold = 2

    # now use a ~5 sigma threshold to represent an initial cutoff for significance
    # note: this is more informative for GemCode than Chromium, because there are many more fragments per BC
    # and thus the expected number of overlaps due to chance is much higher. for Chromium, the threshold will usually
    # be 2.
    overlap_threshold = np.float32(max(min_threshold, round(expected_overlaps + 5 * np.sqrt(expected_overlaps))))
    logger.info("expected overlaps: %f -- using overlap threshold: %f" % (expected_overlaps, overlap_threshold))

    # Chunk out matrix in x and y and find significant overlaps in each pair of chunks
    bc_chunk_size = 1000

    # Choose a p-value that accounts for many comparisons
    fpr = 1
    n_comparisons = nbarcodes ** 2
    pvalue_cut = fpr / n_comparisons  # bernoulli corrected p-value threshold
    logger.info("using pvalue cutoff: %e" % pvalue_cut)

    # Compare the barcodes in each chunk towards themselfs and all other chunks.
    for x in tqdm(range(0, nbarcodes, bc_chunk_size), desc="Comparing chunks"):
        for y in range(x, nbarcodes, bc_chunk_size):
            window_intersection_slices(matrix, x, y, bc_chunk_size, pvalue_cut, overlap_threshold, num_positions,
                                       effective_positions, index_barcode, barcode_graph)

    # Delete the fragment matrix to save memory
    del matrix

    # -------------------------------------------------------------------------------------------------
    # End of code from coalescence.py
    # -------------------------------------------------------------------------------------------------

    # Write outputs
    if output_pickle:
        logger.info(f"Writing pickle object to {output_pickle}")
        with open(output_pickle, 'wb') as file:
            pickle.dump(barcode_graph.graph, file, pickle.HIGHEST_PROTOCOL)

    merge_dict = barcode_graph.get_merges()
    summary["Barcodes removed"] = len(merge_dict)

    if output_merges:
        logger.info(f"Writing merges to {output_merges}")
        with open(output_merges, 'w') as file:
            for old_barcode, new_barcode in merge_dict.items():
                print(old_barcode, new_barcode, sep=",", file=file)

    logger.info("Finished")
    print_stats(summary, name=__name__)


def window_intersection_slices(matrix, bc_offset1, bc_offset2, bc_chunk_size, pvalue_threshold, overlap_threshold,
                               num_positions, effective_positions, index_barcode, barcode_graph):
    """Count the window intersection for bc slices of the full matrix"""

    nbarcodes, _ = matrix.shape
    slice1 = matrix[bc_offset1:min(nbarcodes, bc_offset1 + bc_chunk_size)]
    slice2 = matrix[bc_offset2:min(nbarcodes, bc_offset2 + bc_chunk_size)]

    # Multiply matrix  transpose, non-zero position will indicate overlaps
    shared_positions = slice1 * slice2.T

    # overlap matrix is normalized so that the average value for uncorrelated BCs
    # is 1.  Use that to pick a threshold of 'interesting' BC pairs

    # TODO How do this work?
    # Take a very conservative bite to reduce false positives
    overlapping_bcs = shared_positions > overlap_threshold
    _i1, _i2 = overlapping_bcs.nonzero()
    first = _i1 < _i2

    # scipy.sparse has bug for empty selection -- bail out here to avoid it
    if first.sum() == 0:
        return

    shared_positions_pass = np.array(shared_positions[_i1[first], _i2[first]]).flatten()

    i1s = _i1[first] + bc_offset1
    i2s = _i2[first] + bc_offset2

    for (elem_idx, (i1, i2)) in enumerate(zip(i1s, i2s)):
        # Compute the p-value of the overlaps under a null model with uncorrelated positions
        lmbda = float(num_positions[i1] * num_positions[i2]) / effective_positions
        overlaps = shared_positions_pass[elem_idx]

        # Do a quick check that we have a significant hit before doing detailed p-value calculation
        if overlaps < lmbda + 6 * np.sqrt(lmbda):
            continue

        pvalue = stats.poisson.sf(overlaps, lmbda)
        if pvalue < pvalue_threshold:
            barcode_graph.add_connected_barcodes({index_barcode[i1], index_barcode[i2]})


def paired_reads(path: str, mapq_threshold: int, summary):
    """
    Yield (forward_read, reverse_read) pairs for all properly paired read pairs in the input file.

    :param path: str, path to SAM file
    :param mapq_threshold: int
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
                if pair_is_mapped_and_proper(read, mapq_threshold, summary):
                    cache[read.query_name] = read
                continue
            if pair_orientation_is_fr(read, mate, summary):
                yield read, mate


def pair_is_mapped_and_proper(read: AlignedSegment, mapq_threshold: int, summary) -> bool:
    """
    Checks so read pair meets requirements before being used in analysis.
    :param read: pysam read
    :param mapq_threshold: int
    :param summary: dict
    :return: bool
    """
    if read.is_unmapped:
        summary["Unmapped reads"] += 1
        return False

    if read.mate_is_unmapped:
        summary["Unmapped reads"] += 1
        return False

    if not read.is_proper_pair:
        summary["Reads not proper pair"] += 1
        return False

    if read.mapping_quality < mapq_threshold:
        summary[f"Reads with MAPQ<{mapq_threshold}"] += 1
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


def find_barcode_duplicates(positions, dup_pos):
    """
    Parse positions to check for valid duplicate positions that can be quired to find barcodes to merge
    :param positions: list: Position to check for duplicates
    :param dup_pos: list: Tracks previous duplicate positions and their barcode sets.
    """
    positions_to_remove = list()
    for position, tracked_position in positions.items():
        if tracked_position.has_updated_barcodes:
            if tracked_position.is_duplicate():
                dup_pos[position] |= tracked_position.barcodes
            tracked_position.has_updated_barcodes = False
        else:
            positions_to_remove.append(position)

    for position in positions_to_remove:
        del positions[position]


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
        self.has_updated_barcodes = True
        self.reads += 2
        self.barcodes.add(barcode)

    def is_duplicate(self) -> bool:
        return len(self.barcodes) > 1


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
        merges_dict = dict()
        for component in self.components():
            barcodes_sorted = sorted(component)
            barcode_min = barcodes_sorted[0]
            merges_dict.update({barcode: barcode_min for barcode in barcodes_sorted[1:]})
        return merges_dict

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
        "-w", "--window", type=int, default=100000,
        help="Window size. Duplicate positions within this distance will be used to find cluster "
        "duplicates. Default: %(default)s")
    parser.add_argument(
        "--buffer-size", type=int, default=200,
        help="Buffer size for collecting duplicates. Should be around read length. "
        "Default: %(default)s")
    parser.add_argument(
        "-m", "--mapping-quality", type=int, default=0,
        help="Mapping quality threshold. Default: %(default)s")

