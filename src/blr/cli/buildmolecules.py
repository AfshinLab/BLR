"""
Tags SAM/BAM file with molecule information based on barcode sequence and genomic proximity.

A molecule is defined by having 1) minimum --threshold reads and including all reads with the same barcode which are 2)
a maximum distance of --window between any given reads.
"""

from collections import defaultdict
import logging
import statistics

import pandas as pd
import pysam

from blr.utils import PySAMIO, get_bamtag, Summary, calculate_N50, tqdm, ACCEPTED_LIBRARY_TYPES, LastUpdatedOrderedDict

logger = logging.getLogger(__name__)

# For reads not associated to a specific molecule the molecule id is set to -1.
DEFAULT_MOLECULE_ID = -1


def main(args):
    run_buildmolecules(
        input=args.input,
        output=args.output,
        threshold=args.threshold,
        window=args.window,
        barcode_tag=args.barcode_tag,
        stats_tsv=args.stats_tsv,
        molecule_tag=args.molecule_tag,
        min_mapq=args.min_mapq,
        library_type=args.library_type
    )


def run_buildmolecules(
    input: str,
    output: str,
    threshold: int,
    window: int,
    barcode_tag: str,
    stats_tsv: str,
    molecule_tag: str,
    min_mapq: int,
    library_type: str
):
    summary = Summary()

    # Build molecules from BCs and reads
    save = pysam.set_verbosity(0)  # Fix for https://github.com/pysam-developers/pysam/issues/939
    with pysam.AlignmentFile(input, "rb") as infile:
        bc_to_mol_dict, header_to_mol_dict = build_molecules(pysam_openfile=infile,
                                                             barcode_tag=barcode_tag,
                                                             window=window,
                                                             min_reads=threshold,
                                                             library_type=library_type,
                                                             min_mapq=min_mapq,
                                                             summary=summary)
    pysam.set_verbosity(save)

    # Writes filtered out
    with PySAMIO(input, output, __name__) as (openin, openout):
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            name = read.query_name

            molecule_id = header_to_mol_dict.get(name, DEFAULT_MOLECULE_ID)
            read.set_tag(molecule_tag, molecule_id)

            openout.write(read)

    header_to_mol_dict.clear()

    # Write molecule/barcode file stats
    if stats_tsv:
        logger.info(f"Writing {stats_tsv}")
        df = compute_molecule_stats_dataframe(bc_to_mol_dict)
        df.to_csv(stats_tsv, sep="\t", index=False)
        if not df.empty:
            update_summary_from_molecule_stats(df, summary)
    summary.print_stats(name=__name__)


def parse_reads(pysam_openfile, barcode_tag, min_mapq, summary):
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):
        summary["Total reads"] += 1
        if read.is_duplicate or read.is_unmapped or read.mapping_quality < min_mapq:
            summary["Non analyced reads"] += 1
            continue

        barcode = get_bamtag(pysam_read=read, tag=barcode_tag)
        if not barcode:
            summary["Non analyced reads"] += 1
            continue

        yield barcode, read


def build_molecules(pysam_openfile, barcode_tag, window, min_reads, library_type, min_mapq, summary):
    """
    Builds all_molecules.bc_to_mol ([barcode][moleculeID] = molecule) and
    all_molecules.header_to_mol ([read_name]=mol_ID)
    :param pysam_openfile: Pysam open file instance.
    :param barcode_tag: Tag used to store barcode in bam file.
    :param window: Max distance between reads to include in the same molecule.
    :param min_reads: Minimum reads to include molecule in all_molecules.bc_to_mol
    :param library_type: str. Library construction method.
    :param min_mapq: int
    :param summary: dict for stats collection
    :return: dict[barcode][molecule] = moleculeInstance, dict[read_name] = mol_ID
    """

    all_molecules = AllMolecules(min_reads=min_reads, window=window, library_type=library_type)

    prev_chrom = None
    prev_window_stop = window
    logger.info("Dividing barcodes into molecules")
    for barcode, read in parse_reads(pysam_openfile, barcode_tag, min_mapq, summary):
        # Commit molecules between chromosomes
        if prev_chrom != read.reference_name:
            all_molecules.report_and_remove_all()
            prev_chrom = read.reference_name

        all_molecules.assign_read(read, barcode, summary)

        if read.reference_start > prev_window_stop:
            all_molecules.update_cache(read.reference_start)
            prev_window_stop = read.reference_start + window

    all_molecules.report_and_remove_all()

    return all_molecules.bc_to_mol, all_molecules.header_to_mol


class Molecule:
    """
    A Splitting of barcode read groups into several molecules based on mapping proximity. Equivalent to several
    molecules being barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """
    molecule_counter = 0

    def __init__(self, read, barcode, id=None):
        """
        :param read: pysam.AlignedSegment
        :param barcode: barcode ID
        """
        self.barcode = barcode
        self.start = read.reference_start
        self.stop = read.reference_end
        self.read_headers = {read.query_name}
        self.number_of_reads = 1
        self.bp_covered = self.stop - self.start

        Molecule.molecule_counter += 1
        self.id = Molecule.molecule_counter if id is None else id

    def length(self):
        return self.stop - self.start

    def add_read(self, read):
        """
        Updates molecule's stop position, number of reads and header name set()
        """
        self.bp_covered += max(read.reference_end - max(read.reference_start, self.stop), 0)
        self.stop = max(read.reference_end, self.stop)
        self.read_headers.add(read.query_name)

        self.number_of_reads += 1

    def has_acceptable_overlap(self, read, library_type, summary):
        if read.query_name in self.read_headers:  # Within pair
            return True

        if self.stop < read.reference_start:  # No overlap
            return True

        # If Tn5 transposase was used for library construction, overlaps of ~9 bp are accepted.
        if library_type in {'dbs', 'blr', 'stlfr'} and 8 <= self.stop - read.reference_start <= 10:  # TODO Remove blr
            summary["Tn5-overlapping reads"] += 1
            return True

        # If MuA transposase was used for library construction, overlaps of ~5 bp are accepted.
        if library_type == 'tellseq' and 4 <= self.stop - read.reference_start <= 6:
            summary["MuA-overlapping reads"] += 1
            return True

        # Overlapping reads are ok for 10x genomics libraries.
        if library_type == "10x":
            return True

        summary["Overlapping reads in molecule"] += 1
        return False

    def to_dict(self):
        return {
            "MoleculeID": self.id,
            "Barcode": self.barcode,
            "Reads": self.number_of_reads,
            "Length": self.length(),
            "BpCovered": self.bp_covered,
        }


class AllMolecules:
    """
    Tracks all molecule information, with finished molecules in .bc_to_mol, and molecules which still might get more
    reads in .molecule_cache.
    """

    def __init__(self, min_reads, window, library_type):
        """
        :param min_reads: Minimum reads required to add molecule to .bc_to_mol from .cache_dict
        :param window: Current window for detecting molecules.
        :param library_type: str. Library construction method
        """

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Window for calling molecule
        self.window = window

        # Bool for checking Tn5 overlaps
        self.library_type = library_type

        # Molecule tracking system
        self.molecule_cache = LastUpdatedOrderedDict()

        # Dict for finding mols belonging to the same BC
        self.bc_to_mol = defaultdict(list)

        # Dict for finding mol ID when writing out
        self.header_to_mol = dict()

    def assign_read(self, read, barcode, summary):
        """
        Assign read to current molecule while checking for overlaps or start new molecule.
        """
        if barcode in self.molecule_cache:
            if self.read_is_in_window(read, barcode):
                if self.check_overlaps_ok(read, barcode, summary):
                    self.add_read_to_molecule(read, barcode)
            else:
                self.report(barcode)
                self.terminate(barcode)

                self.create_new_molecule(read, barcode)
        else:
            self.create_new_molecule(read, barcode)

    def check_overlaps_ok(self, read, barcode, summary):
        """
        Check if overlap between current read and those in the molecule are acceptable.
        """
        molecule = self.molecule_cache[barcode]
        return molecule.has_acceptable_overlap(read, self.library_type, summary)

    def read_is_in_window(self, read, barcode):
        """
        Check if read is within window of molecule in cache
        """
        return self.molecule_cache[barcode].stop + self.window >= read.reference_start

    def add_read_to_molecule(self, read, barcode):
        """
        Add read to existing molecule
        """
        self.molecule_cache[barcode].add_read(read)
        self.molecule_cache.move_to_end(barcode)

    def create_new_molecule(self, read, barcode):
        """
        Create new molecule and add to cache
        """
        self.molecule_cache[barcode] = Molecule(read=read, barcode=barcode)

    def update_cache(self, current_start):
        """
        Go through cache and report any molecules outside current window
        """
        window_start = current_start - self.window
        barcodes_to_remove = set()
        for barcode, molecule in self.molecule_cache.items():
            if molecule.stop < window_start:
                self.report(barcode)
                barcodes_to_remove.add(barcode)
            else:
                break

        for barcode in barcodes_to_remove:
            del self.molecule_cache[barcode]

    def report(self, barcode):
        """
        Commit molecule to .bc_to_mol, if molecule.reads >= min_reads. If molecule in cache only barcode is required.
        """
        molecule = self.molecule_cache[barcode]
        if molecule.number_of_reads >= self.min_reads:
            self.bc_to_mol[barcode].append(molecule.to_dict())
            self.header_to_mol.update(
                {header: molecule.id for header in molecule.read_headers}
            )

    def terminate(self, barcode):
        """
        Removes a specific molecule from .molecule_cache
        """
        del self.molecule_cache[barcode]

    def report_and_remove_all(self):
        """
        Commit all .molecule_cache molecules to .bc_to_mol and empty .molecule_cache (provided they meet criterias by
        report function).
        """
        for barcode in self.molecule_cache:
            self.report(barcode)
        self.molecule_cache.clear()


def compute_molecule_stats_dataframe(bc_to_mol_dict):
    """
    Writes stats file for molecules and barcode with information like how many reads, barcodes, molecules etc they
    have
    """
    molecule_data = list()
    while bc_to_mol_dict:
        _, molecules = bc_to_mol_dict.popitem()
        molecule_data.extend(molecules)

    return pd.DataFrame(molecule_data)


def update_summary_from_molecule_stats(df, summary):
    summary["Fragment N50 (bp)"] = calculate_N50(df["Length"])
    summary["Mean fragment size (bp)"] = statistics.mean(df["Length"])
    summary["Median fragment size (bp)"] = statistics.median(df["Length"])
    summary["Longest fragment (bp)"] = max(df["Length"])
    summary["Mean fragment read coverage (%)"] = statistics.mean(100 * df["BpCovered"] / df["Length"])
    summary["Median fragment read coverage (%)"] = statistics.median(100 * df["BpCovered"] / df["Length"])


def add_arguments(parser):
    parser.add_argument(
        "input",
        help="Sorted SAM/BAM file tagged with barcode in the same tag as specified in "
             "-b/--barcode-tag."
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Write output BAM to file rather then stdout."
    )
    parser.add_argument(
        "-t", "--threshold", type=int, default=4,
        help="Threshold for how many reads are required for including given molecule in statistics "
             "(except_reads_per_molecule). Default: %(default)s"
    )
    parser.add_argument(
        "-w", "--window", type=int, default=30000,
        help="Window size cutoff for maximum distance in between two reads in one molecule. Default: "
             "%(default)s"
    )
    parser.add_argument(
        "-b", "--barcode-tag", default="BX",
        help="SAM tag for storing the error corrected barcode. Default: %(default)s"
    )
    parser.add_argument(
        "-s", "--stats-tsv", metavar="FILE",
        help="Write molecule stats in TSV format to FILE"
    )
    parser.add_argument(
        "-m", "--molecule-tag", default="MI",
        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
             "Default: %(default)s"
    )
    parser.add_argument(
        "--min-mapq", type=int, default=0,
        help="Minimum mapping-quality to include reads in analysis Default: %(default)s"
    )
    parser.add_argument(
        "-l", "--library-type", default="dbs", choices=ACCEPTED_LIBRARY_TYPES,
        help="Select library type from currently available technologies: %(choices)s. Default: %(default)s"
    )
