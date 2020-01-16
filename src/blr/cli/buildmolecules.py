"""
Tags SAM/BAM file with molecule information based on barcode sequence and genomic proximity.

A molecule is defined by having 1) minimum --threshold reads and including all reads with the same barcode which are 2)
a maximum distance of --window between any given reads.
"""

import pysam
import logging
from collections import Counter
from tqdm import tqdm

from blr import utils

logger = logging.getLogger(__name__)


def main(args):
    summary = Counter()

    # Build molecules from BCs and reads
    with pysam.AlignmentFile(args.input, "rb") as infile:
        bc_to_mol_dict, header_to_mol_dict = build_molecules(pysam_openfile=infile,
                                                             barcode_tag=args.barcode_tag,
                                                             window=args.window,
                                                             min_reads=args.threshold,
                                                             summary=summary)

    # Writes filtered out
    with pysam.AlignmentFile(args.input, "rb") as openin, \
            pysam.AlignmentFile(args.output, "wb", template=openin) as openout:
        logger.info("Writing filtered bam file")
        for read in tqdm(openin.fetch(until_eof=True)):
            name = read.query_name
            barcode = utils.get_bamtag(pysam_read=read, tag=args.barcode_tag)

            # If barcode is not in bc_to_mol_dict the barcode does not have enough proximal reads to make a single
            # molecule.
            if barcode in bc_to_mol_dict:
                summary[f"Output reads {args.number_tag} tagged"] += 1
                bc_num_molecules = len(bc_to_mol_dict[barcode])
                read.set_tag(args.number_tag, bc_num_molecules)

            # If the read name is in header_to_mol_dict then it is associated to a specific molecule.
            if name in header_to_mol_dict:
                summary[f"Output reads {args.molecule_tag} tagged"] += 1
                molecule_ID = header_to_mol_dict[name]
                read.set_tag(args.molecule_tag, molecule_ID)

            openout.write(read)

    utils.print_stats(summary, name=__name__)

    # Write molecule/barcode file stats
    if args.stats_files:
        logger.info("Writing statistics files")
        write_molecule_stats(bc_to_mol_dict)


def parse_reads(pysam_openfile, barcode_tag, summary):
    for read in tqdm(pysam_openfile.fetch(until_eof=True)):
        is_good = True
        summary["Total reads"] += 1
        if read.is_duplicate:
            summary["Duplicates"] += 1
            is_good = False

        barcode = utils.get_bamtag(pysam_read=read, tag=barcode_tag)
        if not barcode:
            summary["Reads without barcode"] += 1
            is_good = False

        if read.is_unmapped and barcode:
            summary[f"Unmapped {barcode_tag} tagged read"] += 1
            is_good = False

        if is_good:
            yield barcode, read
        else:
            summary["Non analyced reads"] += 1


def build_molecules(pysam_openfile, barcode_tag, window, min_reads, summary):
    """
    Builds all_molecules.bc_to_mol ([barcode][moleculeID] = molecule) and
    all_molecules.header_to_mol ([read_name]=mol_ID)
    :param pysam_openfile: Pysam open file instance.
    :param barcode_tag: Tag used to store barcode in bam file.
    :param window: Max distance between reads to include in the same molecule.
    :param min_reads: Minimum reads to include molecule in all_molecules.bc_to_mol
    :param summary: dict for stats collection
    :return: dict[barcode][molecule] = moleculeInstance, dict[read_name] = mol_ID
    """

    all_molecules = AllMolecules(min_reads=min_reads)

    prev_chrom = pysam_openfile.references[0]
    logger.info("Dividing barcodes into molecules")
    for barcode, read in parse_reads(pysam_openfile, barcode_tag, summary):
        read_start = read.reference_start
        read_stop = read.reference_end

        # Commit molecules between chromosomes
        if not prev_chrom == read.reference_name:
            all_molecules.report_and_remove_all()
            prev_chrom = read.reference_name

        if barcode in all_molecules.cache_dict:
            molecule = all_molecules.cache_dict[barcode]

            # Read is within window => add read to molecule (don't include overlapping reads).
            if (molecule.stop + window) >= read_start:
                if molecule.stop >= read_start and read.query_name not in molecule.read_headers:
                    summary["Overlapping reads in molecule"] += 1
                else:
                    molecule.add_read(stop=read_stop, read_header=read.query_name)
                    all_molecules.cache_dict[barcode] = molecule

            # Read is not within window => report old and initiate new molecule for that barcode.
            else:
                all_molecules.report(molecule=molecule)
                all_molecules.terminate(molecule=molecule)

                molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
                all_molecules.cache_dict[molecule.barcode] = molecule

        else:
            molecule = Molecule(barcode=barcode, start=read_start, stop=read_stop, read_header=read.query_name)
            all_molecules.cache_dict[molecule.barcode] = molecule

    all_molecules.report_and_remove_all()

    return all_molecules.bc_to_mol, all_molecules.header_to_mol


class Molecule:
    """
    A Splitting of barcode read groups into several molecules based on mapping proximity. Equivalent to several
    molecules being barcoded simultaneously in the same emulsion droplet (meaning with the same barcode).
    """

    molecule_counter = int()

    def __init__(self, barcode, start, stop, read_header):
        """
        :param barcode: barcode ID
        :param start: min(read_mapping_positions)
        :param stop: max(read_mapping_positions)
        :param read_header: read ID
        """
        self.barcode = barcode
        self.start = start
        self.stop = stop
        self.read_headers = {read_header}
        self.number_of_reads = 1

        Molecule.molecule_counter += 1
        self.ID = Molecule.molecule_counter

    def length(self):
        return self.stop - self.start

    def add_read(self, stop, read_header):
        """
        Updates molecule's stop position, number of reads and header name set()
        """

        self.stop = stop
        self.read_headers.add(read_header)

        self.number_of_reads += 1


class AllMolecules:
    """
    Tracks all molecule information, with finished molecules in .bc_to_mol, and molecules which still might get more
    reads in .cache_dict.
    """

    def __init__(self, min_reads):
        """
        :param min_reads: Minimum reads required to add molecule to .bc_to_mol from .cache_dict
        """

        # Min required reads for calling proximal reads a molecule
        self.min_reads = min_reads

        # Molecule tracking system
        self.cache_dict = dict()

        # Dict for finding mols belonging to the same BC
        self.bc_to_mol = dict()

        # Dict for finding mol ID when writing out
        self.header_to_mol = dict()

    def report(self, molecule):
        """
        Commit molecule to .bc_to_mol, if molecule.reads >= min_reads
        """

        if molecule.number_of_reads >= self.min_reads:
            if molecule.barcode not in self.bc_to_mol:
                self.bc_to_mol[molecule.barcode] = set()
            self.bc_to_mol[molecule.barcode].add(molecule)
            for header in molecule.read_headers:
                self.header_to_mol[header] = molecule.ID

    def terminate(self, molecule):
        """
        Removes a specific molecule from .cache_dict
        """

        del self.cache_dict[molecule.barcode]

    def report_and_remove_all(self):
        """
        Commit all .cache_dict molecules to .bc_to_mol and empty .cache_dict (provided they meet criterias by report
        function).
        """

        for molecule in self.cache_dict.values():
            self.report(molecule=molecule)
        self.cache_dict = dict()


def write_molecule_stats(molecule_dict):
    """
    Writes stats file for molecules and barcode with information like how many reads, barcodes, molecules etc they
    have
    """

    # Opening all files
    molecules_per_bc = open("molecules_per_bc.tsv", "w")
    molecule_stats = open("molecule_stats.tsv", "w")

    # Write headers
    print(f"Barcode\tNrMolecules", file=molecules_per_bc)
    print(f"Reads\tLength\tBarcode\tNrMolecules",
          file=molecule_stats)

    # Writing molecule-dependant stats
    for barcode in tqdm(molecule_dict):
        number_of_molecules = len(molecule_dict[barcode])
        print(f"{barcode}\t{number_of_molecules}", file=molecules_per_bc)
        for molecule in (molecule_dict[barcode]):
            print(f"{molecule.number_of_reads}\t{molecule.length()}\t{barcode}\t{number_of_molecules}",
                  file=molecule_stats)

    # Close files
    for output_file in (molecules_per_bc, molecule_stats):
        output_file.close()


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
    parser.add_argument("-s", "--stats-files", action="store_true",
                        help="Write barcode/molecule statistics files.")
    parser.add_argument("-m", "--molecule-tag", default="MI",
                        help="SAM tag for storing molecule index specifying a identified molecule for each barcode. "
                             "Default: %(default)s")
    parser.add_argument("-n", "--number-tag", default="MN",
                        help="SAM tag for storing molecule count for a particular barcode. Default: %(default)s")
