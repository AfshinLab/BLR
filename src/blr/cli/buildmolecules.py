"""
Tags SAM/BAM file with molecule information based on barcode sequence and genomic proximity.

A molecule is defined by having 1) minimum --threshold reads and including all reads with the same barcode which are 2)
a maximum distance of --window between any given reads.
"""

import pysam
import logging
from collections import Counter, OrderedDict, defaultdict
import pandas as pd
import statistics

from blr.utils import PySAMIO, get_bamtag, print_stats, calculate_N50, tqdm

logger = logging.getLogger(__name__)


def main(args):
    summary = Counter()

    # Build molecules from BCs and reads
    with pysam.AlignmentFile(args.input, "rb") as infile:
        library_type = infile.header.to_dict()["RG"][0]["LB"]
        bc_to_mol_dict, header_to_mol_dict = build_molecules(pysam_openfile=infile,
                                                             barcode_tag=args.barcode_tag,
                                                             window=args.window,
                                                             min_reads=args.threshold,
                                                             tn5=library_type == "blr",
                                                             summary=summary)
    # Writes filtered out
    with PySAMIO(args.input, args.output, __name__) as (openin, openout):
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            name = read.query_name
            barcode = get_bamtag(pysam_read=read, tag=args.barcode_tag)

            # If barcode is not in bc_to_mol_dict the barcode does not have enough proximal reads to make a single
            # molecule.
            bc_num_molecules = len(bc_to_mol_dict.get(barcode, {}))
            read.set_tag(args.number_tag, bc_num_molecules)

            # If the read name is in header_to_mol_dict then it is associated to a specific molecule.
            # For reads not associated to a specific molecule the molecule id is set to -1.
            molecule_id = header_to_mol_dict.get(name, -1)
            read.set_tag(args.molecule_tag, molecule_id)

            openout.write(read)

    header_to_mol_dict.clear()

    # Write molecule/barcode file stats
    if args.stats_tsv:
        logger.info(f"Writing {args.stats_tsv}")
        df = compute_molecule_stats_dataframe(bc_to_mol_dict)
        df.to_csv(args.stats_tsv, sep="\t", index=False)
        if not df.empty:
            update_summary_from_molecule_stats(df, summary)
    print_stats(summary, name=__name__)


def parse_reads(pysam_openfile, barcode_tag, summary):
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):
        summary["Total reads"] += 1
        if read.is_duplicate or read.is_unmapped:
            summary["Non analyced reads"] += 1
            continue

        barcode = get_bamtag(pysam_read=read, tag=barcode_tag)
        if not barcode:
            summary["Non analyced reads"] += 1
            continue

        yield barcode, read


def build_molecules(pysam_openfile, barcode_tag, window, min_reads, tn5, summary):
    """
    Builds all_molecules.bc_to_mol ([barcode][moleculeID] = molecule) and
    all_molecules.header_to_mol ([read_name]=mol_ID)
    :param pysam_openfile: Pysam open file instance.
    :param barcode_tag: Tag used to store barcode in bam file.
    :param window: Max distance between reads to include in the same molecule.
    :param min_reads: Minimum reads to include molecule in all_molecules.bc_to_mol
    :param tn5: boolean. Library is constructed using Tn5 transposase and has possible 9-bp overlaps.
    :param summary: dict for stats collection
    :return: dict[barcode][molecule] = moleculeInstance, dict[read_name] = mol_ID
    """

    all_molecules = AllMolecules(min_reads=min_reads, window=window, tn5=tn5)

    prev_chrom = pysam_openfile.references[0]
    logger.info("Dividing barcodes into molecules")
    for nr, (barcode, read) in enumerate(parse_reads(pysam_openfile, barcode_tag, summary)):
        # Commit molecules between chromosomes
        if not prev_chrom == read.reference_name:
            all_molecules.report_and_remove_all()
            prev_chrom = read.reference_name

        all_molecules.assign_read(read, barcode, summary)

        if nr % 5000 == 0:
            all_molecules.update_cache(read.reference_start)

    all_molecules.report_and_remove_all()

    return all_molecules.bc_to_mol, all_molecules.header_to_mol


class Molecule:
    """
    A Splitting of barcode read groups into several molecules based on mapping proximity. Equivalent to several
    molecules being barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """

    molecule_counter = int()

    def __init__(self, read, barcode):
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
        self.id = Molecule.molecule_counter

    def length(self):
        return self.stop - self.start

    def add_read(self, read):
        """
        Updates molecule's stop position, number of reads and header name set()
        """
        self.bp_covered += read.reference_end - max(read.reference_start, self.stop)
        self.stop = read.reference_end
        self.read_headers.add(read.query_name)

        self.number_of_reads += 1

    def has_acceptable_overlap(self, read, tn5, summary):
        if read.query_name in self.read_headers:  # Within pair
            return True

        if self.stop < read.reference_start:  # No overlap
            return True

        # If tn5 was used for library construction, overlaps of ~9 bp are accepted.
        if tn5 and 8 <= self.stop - read.reference_start <= 10:
            summary["Tn5-overlapping reads"] += 1
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

    def __init__(self, min_reads, window, tn5):
        """
        :param min_reads: Minimum reads required to add molecule to .bc_to_mol from .cache_dict
        :param window: Current window for detecting molecules.
        :param tn5: bool. If library is of Tn5 type.
        """

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Window for calling molecule
        self.window = window

        # Bool for checking Tn5 overlaps
        self.tn5 = tn5

        # Molecule tracking system
        self.molecule_cache = OrderedDict()

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
        return molecule.has_acceptable_overlap(read, self.tn5, summary)

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
        barcode, molecules = bc_to_mol_dict.popitem()
        nr_molecules = len(molecules)
        for molecule in molecules:
            molecule["NrMolecules"] = nr_molecules
            molecule_data.append(molecule)

    return pd.DataFrame(molecule_data)


def update_summary_from_molecule_stats(df, summary):
    summary["Fragment N50 (bp)"] = calculate_N50(df["Length"])
    summary["Mean fragment size (bp)"] = statistics.mean(df["Length"])
    summary["Median fragment size (bp)"] = statistics.median(df["Length"])
    summary["Longest fragment (bp)"] = max(df["Length"])
    summary["Mean fragment bp covered by reads"] = statistics.mean(df["BpCovered"] / df["Length"])
    summary["Median fragment bp covered by reads"] = statistics.median(df["BpCovered"] / df["Length"])


def add_arguments(parser):
    parser.add_argument("input",
                        help="Sorted SAM/BAM file tagged with barcode in the same tag as specified in "
                             "-b/--barcode-tag.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output BAM to file rather then stdout.")
    parser.add_argument("-t", "--threshold", type=int, default=4,
                        help="Threshold for how many reads are required for including given molecule in statistics "
                             "(except_reads_per_molecule). Default: %(default)s")
    parser.add_argument("-w", "--window", type=int, default=30000,
                        help="Window size cutoff for maximum distance in between two reads in one molecule. Default: "
                             "%(default)s")
    parser.add_argument("-b", "--barcode-tag", default="BX",
                        help="SAM tag for storing the error corrected barcode. Default: %(default)s")
    parser.add_argument("-s", "--stats-tsv", metavar="FILE",
                        help="Write molecule stats in TSV format to FILE")
    parser.add_argument("-m", "--molecule-tag", default="MI",
                        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
                             "Default: %(default)s")
    parser.add_argument("-n", "--number-tag", default="MN",
                        help="SAM tag for storing molecule count for a particular barcode. Default: %(default)s")
