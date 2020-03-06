"""
Phases BAM files according to phased SNVs from a HapCUT2 output .phase file.

Requirements:
  - HapCUT2 .phase output
  - BAM file:
    - molecule tags
    - sorted

Outputs:
  - BAM file with haplotype and phase set tag, where that combination corresponds to a continuously phased region.
"""

import pysam
import logging
from tqdm import tqdm
from collections import Counter, OrderedDict

from blr.utils import print_stats, get_bamtag, PySAMIO

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting analysis")
    anomaly_file = open(args.anomaly_file_name, 'w')
    summary = Counter()

    # Save hetSNV info
    phased_snv_dict = get_phased_snvs_hapcut2_format(args.hapcut2_phase_file, args.min_phred_switch_error,
                                                     args.discard_pruning, summary)

    # Phase reads & corresponding molecules at hetSNV sites
    read_phase_dict, molecule_phase_dict = phase_reads_and_molecules(args.input_bam, args.molecule_tag,
                                                                     phased_snv_dict, anomaly_file, summary)

    # Write output (setting phasing)
    with PySAMIO(args.input_bam, args.output, __name__) as (infile, out):
        for read in tqdm(infile.fetch(until_eof=True), desc="Writing output", unit=" reads",
                         total=summary["Total reads"]):
            molecule_phase = None
            read_phase = None

            # Get molecule and read phasing info
            molecule = get_bamtag(read, args.molecule_tag)
            if molecule and molecule in molecule_phase_dict:
                molecule_phase = molecule_phase_dict[molecule]
            if read.query_name in read_phase_dict:
                read_phase = read_phase_dict[read.query_name]

            # Choose phase from the read/molecule phasing information
            phase = decide_haplotype(read_phase, molecule_phase, read, anomaly_file, summary)
            if phase:
                summary["Total reads phased"] += 1
                ps_tag = phase[0]
                hp_tag = int(phase[1])
                read.set_tag("HP", hp_tag)
                read.set_tag("PS", ps_tag)

            out.write(read)

    anomaly_file.close()
    print_stats(summary, name=__name__)
    logger.info("Finished")


def get_phased_snvs_hapcut2_format(hapcut2_file, min_phred_switch_error, discard_pruning, summary):
    """
    Gets phasing information from HapCUT2 output phase file and turns into dictionary.

    HapCUT2 phase file format:
    BLOCK: <block-header>
    <tsv>   <information>   <for>   <phased>    <snv>   <...>

    Example:
                        snv1  snv2
     Reference   ========A=====C=========
     Haplotype1  --------A-----G---------
     Haplotype2  --------T-----C---------

    HapCUT2 entries:
     snv1: _, 0, 1, _, _, A, T, _, _, _, _, _
     snv2: _, 1, 0, _, _, C, G, _, _, _, _, _

    Will be output as:
     dict[chrom][pos_snv1] = (block_number, A, T)
     dict[chrom][pos_snv2] = (block_number. G, C)

    :param hapcut2_file:
    :param min_phred_switch_error:
    :param discard_pruning:
    :param summary
    :return: dict
    """

    phased_snv_dict = dict()

    with PhaseBlockReader(hapcut2_file) as reader:
        for position in reader:
            print(position)
            # Variant filtering
            if skip_variant(position, min_phred_switch_error, discard_pruning, summary):
                continue

            # TODO: Add support for phased SVs (1 bp deletion, insertion and bigger should currently not work)
            for call in position.variants:
                if len(position.ref) != len(call):
                    summary["HapCUT2 file: Phased SVs (not used)"] += 1
                    continue

            # Add phasing info to dict
            if position.chrom not in phased_snv_dict:
                phased_snv_dict[position.chrom] = OrderedDict()

            phased_snv_dict[position.chrom][position.pos] = (position.phaseblock, position.hap1, position.hap2)
            summary["HapCUT2 file: Usable phased sites"] += 1

        summary["HapCUT2 file: Phase blocks"] = position.phaseblock

    return phased_snv_dict


class PhaseBlockReader:
    def __init__(self, filename):
        self.filename = filename
        self._file = open(self.filename, "r")
        self.block_count = 0

    def __enter__(self):
        return self

    def __iter__(self):
        for line in self._file:
            if line.startswith("BLOCK:"):
                self.block_count += 1
                continue

            if line.strip() == "*" * 8:
                continue

            yield(PhasedSNV(line, self.block_count))

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()


class PhasedSNV:
    """
    Phased SNV information from HapCUT2 output file.
    See: https://github.com/vibansal/HapCUT2/blob/master/outputformat.md for specifics.
    """
    def __init__(self, line, phaseblock):
        self.phaseblock = phaseblock

        line_entries = line.strip().split(maxsplit=11)

        # VCF file index (1-based index of the line in the input VCF describing variant)
        self.index = int(line_entries[0])

        # allele on haploid chromosome copy A (0 means reference allele, 1 means variant allele, - for an unphased
        # variant)
        self.call1 = int(line_entries[1]) if line_entries[1] != "-" else None

        # allele on haploid chromosome copy B (0 means reference allele, 1 means variant allele, - for an unphased
        # variant)
        self.call2 = int(line_entries[2]) if line_entries[2] != "-" else None

        self.chrom = line_entries[3]
        self.pos = int(line_entries[4])
        self.ref = line_entries[5]         # reference allele (allele corresponding to 0 in column 2 or 3)
        self.alt = line_entries[6]         # variant allele (allele corresponding to 1 in column 2 or 3)

        self.variants = [self.ref] + self.alt.split(",")

        self.hap1 = self.variants[self.call1] if self.call1 is not None else None
        self.hap2 = self.variants[self.call2] if self.call2 is not None else None

        self.genotype = line_entries[7]    # VCF genotype field (unedited, directly from original VCF)

        self.is_pruned = int(line_entries[8]) == 1
        self.is_phased = self.call1 is not None and self.call2 is not None

        # switch quality: phred-scaled estimated probability that there is a switch error starting at this SNV (0
        # means switch error is likely, 100 means switch is unlikely)
        self.switch_qual = float(line_entries[9]) if line_entries[9] != "." else None

        # mismatch quality: phred-scaled estimated probability that there is a mismatch [single SNV] error at this
        # SNV (0 means SNV is low quality, 100 means SNV is high quality)
        self.missmatch_qual = float(line_entries[10])

    def __str__(self):
        return ", ".join(f"{k} = {v}" for k, v in vars(self).items())


def skip_variant(position, min_phred_switch_error, discard_pruning, summary):
    # Filter: Unphased variant
    if not position.is_phased:
        summary["HapCUT2 file: Non-phased entries"] += 1
        return True

    # Filter: Phred score for switch error here
    if position.switch_qual and position.switch_qual < min_phred_switch_error:
        summary["HapCUT2 file: Phase entries under min quality"] += 1
        return True

    # Filter: Phred score for mismatch at SNV
    if discard_pruning and position.is_pruned:
        summary["HapCUT2 file: Pruned entries"] += 1
        return True

    return False


def phase_reads_and_molecules(bam_file, molecule_tag, phased_snv_dict, anomaly_file, summary):
    """
    Goes through a SAM file and matches reads at variants to one of the two haplotypes. Then sets the reads molecule
    phase accordingly.
    :param bam_file: file, SAM format
    :param molecule_tag: str, SAM tag
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info
    :anomaly_file: open tsv file for writing weird reads to
    :param summary: Counter instance, from collections
    :return: dict[read.query_name] = haplotype, dict[molecule] = haplotype
    """
    molecule_phase_dict = dict()
    read_phase_dict = dict()
    prev_read_pos = (str(), int(), int())
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        for read in tqdm(infile.fetch(until_eof=True), desc="Phasing reads at hetSNV positions", unit=" reads"):
            summary["Total reads"] += 1

            # Init criteria for read
            if skip_read(read, phased_snv_dict, summary):
                continue

            # Get variants at the current read and remove any variants whenever reads have passed their position
            if prev_read_pos != (read.reference_name, read.reference_start, read.reference_end):
                variants_at_read, phased_snv_dict = get_phased_variants(read, phased_snv_dict)
                prev_read_pos = (read.reference_name, read.reference_start, read.reference_end)

            # Link molecule to phasing information
            if len(variants_at_read) >= 1:

                # Phase read according to the haplotype with most concordant SNVs.
                haplotype = phase_read(read, variants_at_read, summary)
                if not haplotype:
                    continue

                # If equal support for both haplotypes, don't phase
                if len(haplotype.values()) == 2 and len(set(haplotype.values())) == 1:
                    line = "\t".join([read.query_name, "h1=h2"])
                    print(line, file=anomaly_file)
                    summary["Reads with equal support for both haplotypes"] += 1
                else:
                    read_haplotype = max(haplotype, key=haplotype.get)
                    read_phase_dict[read.query_name] = read_haplotype

                # Increment support for molecule belonging to haplotype (except whenever molecule is stripped, == -1)
                molecule = get_bamtag(read, molecule_tag)
                if molecule and molecule != -1:
                    summary["Phased reads with molecule info"] += 1
                    molecule_phase_dict = add_phase_info_to_molecules(molecule_phase_dict, molecule, haplotype)
                else:
                    summary["Phased reads without molecule info"] += 1

    # Phase molecules, decided by the most SNV concordance to a given haplotype
    molecule_phase_dict = decide_molecule_phase(molecule_phase_dict, anomaly_file, summary)

    return read_phase_dict, molecule_phase_dict


def skip_read(read, phased_snv_dict, summary):
    if read.is_unmapped:
        return True
    if read.reference_name not in phased_snv_dict:
        summary["Reads without phased SNV in chr"] += 1
        return True
    # if first variant is after read, continue to next read
    if read.reference_end < first_item(phased_snv_dict[read.reference_name], read.reference_end + 1):
        return True


def first_item(iterable, stop_iteration_value=None):
    """
    Returns first item of a given iterable, or stop_iteration_value if last element has been reached.
    :param iterable: Any iterable
    :param stop_iteration_value: What should be returned if the iterable does not have any entries
    :return:
    """
    try:
        return next(iter(iterable))
    except StopIteration:
        return stop_iteration_value


def get_phased_variants(read, phased_snv_dict):
    """
    Finds all called SNVs withing read.reference_start and read.reference_end. Removes variants downstream of
    read.reference_start from phased_snv_dict.

    :param read: pysam read alignment
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info, must be sorted within chromosomes
    :return: dict[pos][nucleotide] = phase_info, updated phased_snv_dict
    """

    variants_at_read = dict()
    removal_list = list()

    for var_pos, haplotype_info in phased_snv_dict[read.reference_name].items():

        # Variant upstream of read => remove var from dict
        if var_pos < read.reference_start:
            removal_list.append(var_pos)
            continue
        # Variant at read => save pos, keep in dict
        elif read.reference_start <= var_pos and var_pos <= read.reference_end:
            variants_at_read[var_pos] = haplotype_info
        # Variant downstream of read => keep in dict and return SNV positions found
        else:
            break

    # Remove variants upstream of reads
    for pos_to_remove in removal_list:
        del phased_snv_dict[read.reference_name][pos_to_remove]

    return variants_at_read, phased_snv_dict


def phase_read(read, variants_at_read, summary):
    """
    Assigns a read a haplotype from the variants at that read
    :param read: pysam read alignment
    :param variants_at_read: dict, dict[ref_pos][phase_info] = allele_seq
    :param summary: Counter instance, from collections
    :return: str, phase_info for the variant with most support.
    """

    # Translate to read.seq positions to reference positions, where read.seq[seq_pos] is equal to <ref_seq>[ref_pos]
    pos_translations = translate_positions(read=read, ref_pos_list=list(variants_at_read.keys()))
    read_haplotype = Counter()

    for ref_pos, seq_pos in pos_translations.items():

        # TODO: Add handling of deletions
        if not seq_pos:
            continue

        # Using read.query_alignment_sequence compensates for insertions (rather than read.seq)
        nt = read.query_alignment_sequence[seq_pos]
        if nt not in variants_at_read[ref_pos]:
            continue
        elif nt == variants_at_read[ref_pos][1]:
            haplotype = "1"
        elif nt == variants_at_read[ref_pos][2]:
            haplotype = "2"
        phase_set = variants_at_read[ref_pos][0]

        read_haplotype[(phase_set, haplotype)] += 1

    # Read did not fit any of the called haplotypes
    if len(read_haplotype) == 0:
        summary["Reads not fitting called alleles"] += 1
        return

    return read_haplotype


def translate_positions(read, ref_pos_list):
    """
    Translates a reference position to the nucleotide position in an read.seq string, where deletions are accounted
    for.

    Discrepancy problem between Ref pos/read.seq pos. ref start =/= read start & deletions.

    Ref        =======================
    Ref pos     1 3     9
    Read        |-- . . ---> (. = Deletion)
    read pos    0 2     3

    :param ref_pos_list: list with ints, position in reference to be translated
    :param read: pysam read alignment
    :return: dict, dict[ref_pos] = seq_pos
    """

    ref_pos_iterator = iter(ref_pos_list)
    ref_pos = next(ref_pos_iterator)
    pos_translations = dict()

    # Loop adjusting for deletion (extra bases in ref / missing bases in read)
    seq_pos = int()
    for block in read.get_blocks():
        block_start, block_stop = block
        while True:

            # Pos upstream of block, add block len to position
            if block_stop < ref_pos:
                seq_pos += block_stop - block_start
                break
            # Pos between blocks
            elif block_start > ref_pos:
                pos_translations[ref_pos] = None  # TODO: Change this to whatever is used to describe deletions

            # Pos found
            else:
                block_pos = ref_pos - block_start
                pos_translations[ref_pos] = seq_pos + block_pos - 1

            # When last variant, return (even if not last alignment block)
            try:
                ref_pos = next(ref_pos_iterator)
            except StopIteration:
                return pos_translations


def add_phase_info_to_molecules(molecule_phase_dict, molecule, haplotype):
    """

    :param molecule_phase_dict:
    :param haplotype:
    :param summary:
    :return:
    """
    for hap, sup in haplotype.items():
        if molecule not in molecule_phase_dict:
            molecule_phase_dict[molecule] = Counter()
        molecule_phase_dict[molecule][hap] += sup

    return molecule_phase_dict


def decide_molecule_phase(molecule_phase_dict, anomaly_file, summary):
    """

    :param molecule_phase_dict:
    :param anomaly_file:
    :return:
    """

    for molecule, phase_info in molecule_phase_dict.copy().items():
        # Check if there are two max values (equal support for top contending haplotypes)
        if len(phase_info.values()) >= 2 and len(set(sorted(phase_info.values(), reverse=True)[0:2])) == 1:
            line = "\t".join([f"MI:{molecule}", "h1=h2"])
            print(line, file=anomaly_file)
            summary["Moleules with equal support for haplotypes"] += 1
            continue

        molecule_haplotype = max(phase_info, key=phase_info.get)
        molecule_phase_dict[molecule] = molecule_haplotype

    return molecule_phase_dict


def decide_haplotype(read_phase, molecule_phase, read, anomaly_file, summary):
    """

    :param read_phase:
    :param molecule_phase:
    :param anomaly_file:
    :param summary:
    :return:
    """

    phase = None
    if molecule_phase and read_phase:
        if molecule_phase == read_phase:
            phase = read_phase
            summary["Concordant read/mol phasing"] += 1
        else:
            line = "\t".join([read.query_name, "mol_hap!=read_hap"])
            print(line, file=anomaly_file)
            summary["Discordant read/mol phasing"] += 1
    elif read_phase:
        summary["Read phased from read (no mol info)"] += 1
        phase = read_phase
    elif molecule_phase:
        summary["Read phased from only mol info"] += 1
        phase = molecule_phase
    else:
        summary["Reads with no phasing info"] += 1

    return phase


def add_arguments(parser):
    parser.add_argument("input_bam",
                        help="BAM file. To read from stdin use '-'. Must be sorted.")
    parser.add_argument("hapcut2_phase_file",
                        help="Phaseblock file from HapCUT2. Must be sorted.")

    parser.add_argument("-o", "--output", default="-",
                        help="Write output phased BAM to file rather then stdout.")
    parser.add_argument("--anomaly-file-name", default="haplotag_anomalies.tsv",
                        help="File to output information with anomalous reads, such as those with conflicting read/"
                             "molecule phasing information. These will also be written to output but will not have "
                             "any phasing information added to them. Default: &(default)s")
    parser.add_argument("--min-phred-switch-error", default=30, type=float,
                        help="Minimum phred score for switch error at any given variant. Require HapCUT2 to have been"
                             " run using the '--error_analysis_mode 1' option. Default: %(default)s.")
    parser.add_argument("--discard-pruning", default=True,
                        help="Discard phasing events marked as pruned by HapCUT2. Default: %(default)s.")
    parser.add_argument("--molecule-tag", default="MI", help="Molecule SAM tag. Default: %(default)s.")
    parser.add_argument("--haplotype-tag", default="HP", help="Haplotype SAM tag. Default: %(default)s")
    parser.add_argument("--phase-set-tag", default="PS", help="Phase set SAM tag. Default: %(default)s")
