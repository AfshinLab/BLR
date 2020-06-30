"""
Remove barcode duplicates (two different barcode sequences origin to the same droplet, tagging the
same tagmented long molecule) by merging barcode sequences for reads sharing duplicates.

Condition to call barcode duplicate:

Two positions (positions defined as a unique set of read_start, read_stop, mate_start, mate_stop))
at a maximum of W (--window, default 100kbp, between = max(downstream_pos)-max(downstream_pos)) bp
apart sharing more than one barcode (share = union(bc_set_pos1, bc_set_pos2)).
"""

from pysam import AlignmentFile, AlignedSegment
import logging
from collections import Counter, deque, OrderedDict

from blr.utils import PySAMIO, get_bamtag, print_stats, tqdm

logger = logging.getLogger(__name__)


def main(args):
    run_clusterrmdup(
        input=args.input,
        output=args.output,
        merge_log=args.merge_log,
        barcode_tag=args.barcode_tag,
        buffer_size=args.buffer_size,
        window=args.window,
    )


def run_clusterrmdup(
    input: str,
    output: str,
    merge_log: str,
    barcode_tag: str,
    buffer_size: int,
    window: int,
):
    logger.info("Starting Analysis")
    summary = Counter()
    positions = OrderedDict()
    chrom_prev = None
    pos_prev = 0
    merge_dict = dict()
    buffer_dup_pos = deque()
    for read, mate in tqdm(paired_reads(input, summary), desc="Reading pairs"):
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
        current_position = (mate.reference_start, read.reference_end, orientation)

        if abs(pos_new - pos_prev) > buffer_size or chrom_new != chrom_prev:
            find_barcode_duplicates(positions, buffer_dup_pos, merge_dict, window, summary)

            if chrom_new != chrom_prev:
                positions.clear()
                buffer_dup_pos.clear()
                chrom_prev = chrom_new

            pos_prev = pos_new

        if current_position not in positions:
            positions[current_position] = PositionTracker(current_position)
        positions[current_position].add_barcode(barcode)

    # Process last chunk
    find_barcode_duplicates(positions, buffer_dup_pos, merge_dict, window, summary)

    # Remove several step redundancy (5 -> 3, 3 -> 1) => (5 -> 1, 3 -> 1)
    reduce_several_step_redundancy(merge_dict)
    summary["Barcodes removed"] = len(merge_dict)

    # Write outputs
    barcodes_written = set()
    with PySAMIO(input, output, __name__) as (infile, out), \
            open(merge_log, "w") as bc_merge_file:
        print("Previous_barcode,New_barcode", file=bc_merge_file)
        for read in tqdm(infile, desc="Writing output", total=summary["Total reads"]):

            # If read barcode in merge dict, change tag and header to compensate.
            previous_barcode = get_bamtag(pysam_read=read, tag=barcode_tag)
            if previous_barcode in merge_dict:
                summary["Reads with new barcode"] += 1
                new_barcode = str(merge_dict[previous_barcode])
                read.set_tag(barcode_tag, new_barcode, value_type="Z")

                # Merge file writing
                if previous_barcode not in barcodes_written:
                    barcodes_written.add(previous_barcode)
                    print(f"{previous_barcode},{new_barcode}", file=bc_merge_file)

            out.write(read)

    logger.info("Finished")
    print_stats(summary, name=__name__)


def paired_reads(path: str, summary):
    """
    Yield (forward_read, reverse_read) pairs for all properly paired read pairs in the input file.

    :param path: str, path to SAM file
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
                if pair_is_mapped_and_proper(read, summary):
                    cache[read.query_name] = read
                continue
            if pair_orientation_is_fr(read, mate, summary):
                yield read, mate


def pair_is_mapped_and_proper(read: AlignedSegment, summary) -> bool:
    """
    Checks so read pair meets requirements before being used in analysis.
    :param read: pysam read
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


def find_barcode_duplicates(positions, buffer_dup_pos, merge_dict, window: int, summary):
    """
    Parse position to check if they are valid duplicate positions. If so start looking for barcode
    duplicates.
    :param positions: list: Position to check for duplicates
    :param merge_dict: dict: Tracks which barcodes shuold be merged .
    :param buffer_dup_pos: list: Tracks previous duplicate positions and their barcode sets.
    :param window: int: Max distance allowed between positions to call barcode duplicate.
    :param summary: dict
    """
    positions_to_remove = list()
    for position in positions.keys():
        tracked_position = positions[position]
        # Check if position is updated i.e. new reads added, else remove it.
        if tracked_position.is_updated:
            if tracked_position.is_duplicate():
                seed_duplicates(
                    merge_dict=merge_dict,
                    buffer_dup_pos=buffer_dup_pos,
                    position=tracked_position.position,
                    position_barcodes=tracked_position.barcodes,
                    window=window
                )
            tracked_position.is_updated = False
        else:
            positions_to_remove.append(position)

    for position in positions_to_remove:
        del positions[position]


class PositionTracker:
    """
    Stores read pairs information relatec to a position and keeps track if reads/mates are marked
    as duplicate for that set of reads.
    """

    def __init__(self, position):
        self.position = position
        self.reads = 0
        self.barcodes = set()
        self.is_updated = False

    def add_barcode(self, barcode: str):
        self.is_updated = True
        self.reads += 2
        self.barcodes.add(barcode)

    def is_duplicate(self) -> bool:
        return len(self.barcodes) >= 2


def seed_duplicates(merge_dict, buffer_dup_pos, position, position_barcodes, window: int):
    """
    Builds up a merge dictionary for which any keys should be overwritten by their value. Also
    keeps all previous positions saved in a list in which all reads which still are withing the
    window size are saved.
    :param merge_dict: dict: Tracks which barcodes should be merged.
    :param buffer_dup_pos: list: Tracks previous duplicate positions and their barcode sets.
    :param position: tuple: Positions (start, stop) to be analyzed and subsequently saved to buffer.
    :param position_barcodes: seq: Barcodes at analyzed position
    :param window: int: Max distance allowed between postions to call barcode duplicate.
    """

    pos_start_new = position[0]
    # Loop over list to get the positions closest to the analyzed position first. When position
    # are out of the window size of the remaining buffer is removed.
    for index, (compared_position, compared_barcodes) in enumerate(buffer_dup_pos):

        # Skip comparison against self.
        if position == compared_position:
            continue

        compared_position_stop = compared_position[1]
        if compared_position_stop + window >= pos_start_new:
            barcode_intersect = position_barcodes & compared_barcodes

            # If two or more unique barcodes are found, update merge dict
            if len(barcode_intersect) >= 2:
                update_merge_dict(merge_dict, barcode_intersect)
        else:
            # Remove positions outside of window (at end of list) since positions are sorted.
            for i in range(len(buffer_dup_pos) - index):
                buffer_dup_pos.pop()
            break

    # Add new position at the start of the list.
    buffer_dup_pos.appendleft((position, position_barcodes))


def update_merge_dict(merge_dict, barcodes):
    """
    Add new barcodes to merge_dict.
    :param merge_dict: dict: Barcode pairs directing merges
    :param barcodes: set: Barcodes to add
    """
    barcodes_sorted = sorted(barcodes)
    barcodes_real_min = find_min_barcode(barcodes_sorted[0], merge_dict)
    for barcode_to_merge in barcodes_sorted[1:]:
        # Never add give one key multiple values (connect the prev/new values instead)
        if barcode_to_merge in merge_dict:
            previous_min = find_min_barcode(barcode_to_merge, merge_dict)
            if not previous_min == barcodes_real_min:
                if min(barcodes_real_min, previous_min) == barcodes_real_min:
                    merge_dict[previous_min] = barcodes_real_min
                else:
                    merge_dict[barcodes_real_min] = previous_min

        # Normal case, just add high_bc_id => min_bc_id
        else:
            merge_dict[barcode_to_merge] = barcodes_real_min


def find_min_barcode(barcode, merge_dict):
    """
    Goes through merge dict and finds the alphabetically top string for a chain of key-value
    entries. E.g if merge_dict has TAGA => GGAT, GGAT => CTGA, CTGA => ACGA it will return ACGA if
    any of the values CTGA, GGAT, TAGA or ACGA are given.
    :param: bc_minimum: str: Barcode
    :param: merge_dict: dict: Barcode pairs directing merges
    :return: str: alphabetically top barcode string
    """
    while barcode in merge_dict:
        barcode = merge_dict[barcode]
    return barcode


def reduce_several_step_redundancy(merge_dict):
    """
    Takes translation  dict saved in object and makes sure 5->3, 3->1 becomes 5->1, 3->1
    :param merge_dict: "messy" merge_dict
    :return: "clean" merge_dict
    """

    for barcode_to_remove in sorted(merge_dict.keys()):
        merge_dict[barcode_to_remove] = find_min_barcode(barcode_to_remove, merge_dict)
    return merge_dict


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="Coordinate-sorted SAM/BAM file tagged with barcodes.")
    parser.add_argument(
        "merge_log",
        help="CSV log file containing all merges done. File is in format: "
        "{old barcode id},{new barcode id}")
    parser.add_argument(
        "-o", "--output", default="-",
        help="Write output BAM to file rather then stdout.")
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
