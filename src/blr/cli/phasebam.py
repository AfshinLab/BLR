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

from collections import Counter, OrderedDict
import logging

import pysam

from blr.utils import Summary, get_bamtag, PySAMIO, tqdm

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting analysis")
    anomaly_file = open(args.anomaly_file_name, 'w')
    summary = Summary()

    # Save hetSNV info
    phased_snv_dict = get_phased_snvs_hapcut2_format(args.hapcut2_phase_file, args.min_switch_error,
                                                     args.min_missmatch_error, args.include_pruned, summary)

    # Phase reads & corresponding molecules at hetSNV sites
    read_phase_dict, molecule_phase_dict = phase_reads_and_molecules(args.input_bam, args.molecule_tag,
                                                                     phased_snv_dict, anomaly_file, args.min_mapq,
                                                                     summary)

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
            phase = decide_haplotype(read_phase, molecule_phase, read, args.tag_priority, anomaly_file, summary)
            if phase and not read.is_unmapped:  # Don't assign phase to unmapped reads.
                summary["Total reads phased"] += 1
                ps_tag = phase[0]
                hp_tag = phase[1]
                read.set_tag("HP", hp_tag)

                # Tag reads with phase inferred from read phase info with quality.
                if len(phase) == 3:
                    read.set_tag("PC", phase[2])
                read.set_tag("PS", ps_tag)
            elif args.discard_untagged:
                continue

            out.write(read)
    summary["Total reads phased (%)"] = 100 * summary["Total reads phased"]/summary["Total reads"]

    anomaly_file.close()
    summary.print_stats(name=__name__)
    logger.info("Finished")


def get_phased_snvs_hapcut2_format(hapcut2_file, min_switch_error, min_missmatch_error, include_pruned, summary):
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

    :param hapcut2_file: str. Path to hapcut2 phaseblock file
    :param min_switch_error: int. Threshold for switch error rate for variant
    :param min_missmatch_error: int. Threshold for missmatch error rate for variant
    :param include_pruned: bool. Whether to include pruned variants.
    :param summary. dict. Summary instance.
    :return: dict
    """

    phased_snv_dict = dict()

    with PhaseBlockReader(hapcut2_file) as reader:
        for variant in tqdm(reader, desc="Parsing variants"):
            summary["HapCUT2 file: Total variants"] += 1

            # Variant filtering
            if skip_variant(variant, min_switch_error, min_missmatch_error, include_pruned, summary):
                summary["HapCUT2 file: Skipped sites"] += 1
                continue

            # Add phasing info to dict
            if variant.chrom not in phased_snv_dict:
                phased_snv_dict[variant.chrom] = OrderedDict()

            phased_snv_dict[variant.chrom][variant.pos] = (variant.phaseset, variant.hap1, variant.hap2)
            summary["HapCUT2 file: Usable phased sites"] += 1

        summary["HapCUT2 file: Phase blocks"] = reader.block_count

    # Variants are ordered within phasesets but different phasesets might overlap. Thus sorting is required.
    sorted_snvs = {ref: OrderedDict(sorted(snvs.items())) for ref, snvs in phased_snv_dict.items()}

    return sorted_snvs


class PhaseBlockReader:
    """Parser for HapCUT2 phaseblock files.
    See https://github.com/vibansal/HapCUT2/blob/master/outputformat.md for format details.
    File has format:

    ```
    BLOCK: offset: 2 len: 2963 phased: 2223 SPAN: 273194 fragments 13333
    <phaseblock 1 pos 1>
    <phaseblock 1 pos 2>
    ...
    <phaseblock 1 pos n>
    *******
    BLOCK: offset: 27490 len: 33 phased: 19 SPAN: 19698 fragments 32
    <phaseblock 2 pos 1>
    <phaseblock 2 pos 2>
    ...
    <phaseblock 2 pos n>
    *******
    ```
    """
    def __init__(self, filename):
        self.filename = filename
        self._file = open(self.filename, "r")
        self.block_count = 0

    def __enter__(self):
        return self

    def __iter__(self):
        phaseset = None
        for line in self._file:
            if line.startswith("BLOCK:"):
                self.block_count += 1
                phaseset = None
                continue

            if line.strip() == "*" * 8:
                continue

            # Use first block entry to define phaseset.
            if not phaseset:
                first_enty = PhasedVariant(line)
                phaseset = first_enty.phaseset
                yield first_enty
            else:
                yield PhasedVariant(line, phaseset=phaseset)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()


class PhasedVariant:
    """
    Phased variant information from HapCUT2 output file.
    See: https://github.com/vibansal/HapCUT2/blob/master/outputformat.md for specifics.
    """
    def __init__(self, line, phaseset=None, max_phred=100.0):
        line_entries = line.strip().split(maxsplit=11)

        # VCF file index (1-based index of the line in the input VCF describing variant)
        self.index = int(line_entries[0])

        # allele on haploid chromosome copy A (0 means reference allele, 1 means variant allele, - for an unphased
        # variant). In some cases the call can also be >1 if there are multiple ALT on a position.
        self.call1 = int(line_entries[1]) if line_entries[1] != "-" else None

        # allele on haploid chromosome copy B (0 means reference allele, 1 means variant allele, - for an unphased
        # variant)
        self.call2 = int(line_entries[2]) if line_entries[2] != "-" else None

        self.chrom = line_entries[3]
        self.pos = int(line_entries[4]) - 1  # Translate positions to 0-based start (same as pysam).
        self.ref = line_entries[5]         # reference allele (allele corresponding to 0 in column 2 or 3)
        self.alt = line_entries[6]         # variant allele(s) (allele corresponding to 1 in column 2 or 3)

        self.variants = [self.ref] + self.alt.split(",")

        self.hap1 = self.variants[self.call1] if self.call1 is not None else None
        self.hap2 = self.variants[self.call2] if self.call2 is not None else None

        self.genotype = line_entries[7]    # VCF genotype field (unedited, directly from original VCF)

        self.is_pruned = int(line_entries[8]) == 1
        self.is_phased = self.call1 is not None and self.call2 is not None

        # switch quality: phred-scaled estimated probability that there is a switch error starting at this SNV (0
        # means switch error is likely, 100 means switch is unlikely). If '.' (i.e. switch error not calculated) the
        # max_phred value is used.
        self.switch_qual = float(line_entries[9]) if line_entries[9] != "." else max_phred

        # mismatch quality: phred-scaled estimated probability that there is a mismatch [single SNV] error at this
        # SNV (0 means SNV is low quality, 100 means SNV is high quality)
        self.missmatch_qual = float(line_entries[10])

        # phaseset is defined as the first position (1-based) of the current phase block. This is the same notation
        # as used for the PS tag in phased VCF files.
        self.phaseset = phaseset if phaseset else (self.pos + 1)

        # Attempt to normalize variants that are not SNVs.
        if not self.is_snv():
            self.normalize()

    def is_snv(self):
        """Check if varinats are all single-nucleotide variants (SNVs)"""
        return (self.hap1 != self.hap2) and (len(self.hap1) == len(self.hap2) == 1)

    def __str__(self):
        return f"chrom={self.chrom}, position={self.pos}, phaseset={self.phaseset}, haps={self.hap1}|{self.hap2}, " \
            f"phased={self.is_phased}, pruned={self.is_pruned}, switch={self.switch_qual}, miss={self.missmatch_qual}"

    def normalize(self):
        """
        Return a normalized version of this variant.

        Common prefixes and/or suffixes between the reference and alternative allele are removed,
        and the position is adjusted as necessary.

        Ex.
            hap1 = GGGATG --> A
            hap2 = GGGTTG --> T

        Based on function from Whatshap: https://bitbucket.org/whatshap/whatshap/src/v0.18/whatshap/vcf.py
        """
        pos, ref, alt = self.pos, self.hap1, self.hap2
        while len(ref) >= 1 and len(alt) >= 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]

        while len(ref) >= 1 and len(alt) >= 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1

        self.pos, self.hap1, self.hap2 = pos, ref, alt


def skip_variant(variant, min_switch_error, min_missmatch_error, include_pruned, summary):
    # Filter: Unphased variant
    if not variant.is_phased:
        summary["HapCUT2 file: Non-phased entries"] += 1
        return True

    # Filter: Only SNVs
    # TODO: Add support for phased SVs (1 bp deletion, insertion and bigger should currently not work)
    if not variant.is_snv():
        summary["HapCUT2 file: Phased SVs (not used)"] += 1
        return True

    # Filter: Phred score for switch error here
    if variant.switch_qual < min_switch_error:
        summary["HapCUT2 file: Phase entries under min quality"] += 1
        return True

    if variant.missmatch_qual < min_missmatch_error:
        summary["HapCUT2 file: Phase entries under min quality"] += 1
        return True

    # Filter: Phred score for mismatch at SNV
    if not include_pruned and variant.is_pruned:
        summary["HapCUT2 file: Pruned entries"] += 1
        return True

    return False


def phase_reads_and_molecules(bam_file, molecule_tag, phased_snv_dict, anomaly_file, min_mapq, summary):
    """
    Goes through a SAM file and matches reads at variants to one of the two haplotypes. Then sets the reads molecule
    phase accordingly.
    :param bam_file: file, SAM format
    :param molecule_tag: str, SAM tag
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info
    :anomaly_file: open tsv file for writing weird reads to
    :param min_mapq: int
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
            if skip_read(read, phased_snv_dict, min_mapq, summary):
                continue

            # Get variants at the current read and remove any variants whenever reads have passed their position
            if prev_read_pos != (read.reference_name, read.reference_start, read.reference_end):
                variants_at_read = get_phased_variants(read, phased_snv_dict)
                prev_read_pos = (read.reference_name, read.reference_start, read.reference_end)

            # Link molecule to phasing information
            if variants_at_read:
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
                    read_phase_dict[read.query_name] = (*read_haplotype, haplotype[read_haplotype])

                # Increment support for molecule belonging to haplotype (except whenever molecule is stripped, == -1)
                molecule = get_bamtag(read, molecule_tag)
                if molecule and molecule != -1:
                    summary["Phased reads with molecule info"] += 1
                    add_phase_info_to_molecules(molecule_phase_dict, molecule, haplotype)
                else:
                    summary["Phased reads without molecule info"] += 1

    # Phase molecules, decided by the most SNV concordance to a given haplotype
    molecule_phase_dict = decide_molecule_phase(molecule_phase_dict, anomaly_file, summary)

    return read_phase_dict, molecule_phase_dict


def skip_read(read, phased_snv_dict, min_mapq, summary):
    if read.is_unmapped or read.mapping_quality < min_mapq:
        return True

    if read.reference_name not in phased_snv_dict:
        summary["Reads without phased SNV in chr"] += 1
        return True
    # if first variant is after read, continue to next read
    if read.reference_end < next(iter(phased_snv_dict[read.reference_name]), read.reference_end + 1):
        return True


def get_phased_variants(read, phased_snv_dict):
    """
    Finds all called SNVs withing read.reference_start and read.reference_end. Removes variants upstream of
    read.reference_start from phased_snv_dict.

    :param read: pysam read alignment
    :param phased_snv_dict: dict, dict[chrom][pos][nucleotide] = phase_info, must be sorted within chromosomes
    :return: dict[pos][nucleotide] = phase_info
    """

    variants_at_read = dict()
    removal_list = list()

    for var_pos, haplotype_info in phased_snv_dict[read.reference_name].items():

        # Variant upstream of read => remove var from dict
        if var_pos < read.reference_start:
            removal_list.append(var_pos)
            continue
        # Variant at read => save pos, keep in dict
        # Note that reference_end points to one past the last aligned residue
        elif read.reference_start <= var_pos < read.reference_end:
            variants_at_read[var_pos] = haplotype_info
        # Variant downstream of read => keep in dict and return SNV positions found
        else:
            break

    # Remove variants upstream of reads
    for pos_to_remove in removal_list:
        # Remove using .pop() since this returns None on KeyError.
        phased_snv_dict[read.reference_name].pop(pos_to_remove, None)

    return variants_at_read


def phase_read(read, variants_at_read, summary):
    """
    Assigns a read a haplotype from the variants at that read
    :param read: pysam read alignment
    :param variants_at_read: dict, dict[ref_pos][phase_info] = allele_seq
    :param summary: Counter instance, from collections
    :return: str, phase_info for the variant with most support.
    """

    # Translate to read.seq positions to reference positions, where read.seq[seq_pos] is equal to <ref_seq>[ref_pos]
    pos_translations = translate_positions(read=read, ref_pos_iterator=iter(variants_at_read))
    read_haplotype = Counter()

    # Unable to translate positions
    if not pos_translations:
        summary["Reads not fitting called alleles"] += 1
        return

    for ref_pos, seq_pos in pos_translations.items():

        # TODO Handle multi-nucleotide variants.
        # Extract single nucleotide at position.
        nt = read.query_sequence[seq_pos]
        qual = read.query_qualities[seq_pos]

        if nt not in variants_at_read[ref_pos]:
            continue
        elif nt == variants_at_read[ref_pos][1]:
            haplotype = 1
        elif nt == variants_at_read[ref_pos][2]:
            haplotype = 2
        phase_set = variants_at_read[ref_pos][0]

        # TODO how to propertly summarize qualities for haplotype? WhatsHap haplotag does + v.s - for hap 1 vs 2.
        read_haplotype[(phase_set, haplotype)] += qual

    # Read did not fit any of the called haplotypes
    if not read_haplotype:
        summary["Reads not fitting called alleles"] += 1
        return

    return read_haplotype


def translate_positions(read, ref_pos_iterator):
    """
    Translates a reference position to the nucleotide position in an read.seq string, where insertion and deletions
    are accounted for.

    :param read: pysam read alignment
    :param ref_pos_iterator: iterator with positions in reference to be translated
    :return: dict, dict[ref_pos] = seq_pos
    """
    variant_ref_pos = next(ref_pos_iterator, None)
    pos_translations = dict()

    # Loop over pairs containing matched query and reference position.
    # NB! Does not work with padded reads.
    for read_seq_pos, read_ref_pos in read.get_aligned_pairs():
        # TODO Currently skips insertions (read_ref_pos is None) and deletions (read_seq_pos is None)
        if read_ref_pos is None or read_seq_pos is None:
            continue

        while variant_ref_pos < read_ref_pos:
            try:
                variant_ref_pos = next(ref_pos_iterator)
            except StopIteration:
                return pos_translations

        if read_ref_pos == variant_ref_pos:
            pos_translations[variant_ref_pos] = read_seq_pos

    return pos_translations


def add_phase_info_to_molecules(molecule_phase_dict, molecule, haplotype):
    """
    Add phasing information from read to molecule phase dict
    :param molecule_phase_dict:
    :param molecule:
    :param haplotype:
    """
    for hap, sup in haplotype.items():
        if molecule not in molecule_phase_dict:
            molecule_phase_dict[molecule] = Counter()
        molecule_phase_dict[molecule][hap] += sup


def decide_molecule_phase(molecule_phase_dict, anomaly_file, summary):
    """
    Compare support for haplotype and deside on final haplotype. If equal support, skip molecule
    :param molecule_phase_dict:
    :param anomaly_file:
    :param summary
    :return:
    """

    for molecule, phase_info in molecule_phase_dict.copy().items():
        # Check if there are two max values (equal support for top contending haplotypes)
        if len(phase_info.values()) >= 2 and len(set(sorted(phase_info.values(), reverse=True)[0:2])) == 1:
            line = "\t".join([f"MI:{molecule}", "h1=h2"])
            print(line, file=anomaly_file)
            summary["Moleules with equal support for haplotypes"] += 1
            del molecule_phase_dict[molecule]
            continue

        molecule_haplotype = max(phase_info, key=phase_info.get)
        molecule_phase_dict[molecule] = molecule_haplotype

    return molecule_phase_dict


def decide_haplotype(read_phase, molecule_phase, read, tag_priority, anomaly_file, summary):
    """
    Compare haplotype information between read and molecule and deside on final haplotype to tag.
    :param read_phase:
    :param molecule_phase:
    :param read:
    :param tag_priority:
    :param anomaly_file:
    :param summary:
    :return:
    """

    phase = None
    if molecule_phase and read_phase:
        concordant = True if molecule_phase == read_phase[:2] else False

        if concordant:
            summary["Read with concordant read/mol phase"] += 1
        else:
            summary["Read with discordant read/mol phase"] += 1

        if tag_priority == "CONCORDANT":

            if concordant:
                phase = read_phase
            else:
                line = "\t".join([read.query_name, "mol_hap!=read_hap"])
                print(line, file=anomaly_file)

        if tag_priority == "READ":
            summary["Read phased from read info"] += 1
            phase = read_phase

        if tag_priority == "MOLECULE":
            summary["Read phased from mol info"] += 1
            phase = molecule_phase

    elif read_phase:
        summary["Read phased from read info"] += 1
        phase = read_phase
    elif molecule_phase:
        summary["Read phased from mol info"] += 1
        phase = molecule_phase
    else:
        summary["Reads with no phasing info"] += 1

    return phase


def add_arguments(parser):
    arg = parser.add_argument
    arg("input_bam",
        help="BAM file. To read from stdin use '-'. Must be sorted.")
    arg("hapcut2_phase_file",
        help="Phaseblock file from HapCUT2.")

    arg("-o", "--output", default="-",
        help="Write output phased BAM to file rather then stdout.")
    arg("--anomaly-file-name", default="haplotag_anomalies.tsv",
        help="File to output information with anomalous reads, such as those with conflicting read/"
             "molecule phasing information. These will also be written to output but will not have "
             "any phasing information added to them. Default: %(default)s")
    arg("--min-switch-error", default=30, type=float, choices=range(0, 101), metavar="[0-100]",
        help="Minimum phred score for switch error at any given variant. Require HapCUT2 to have been"
             " run using the '--error_analysis_mode 1' option. Default: %(default)s.")
    arg("--min-missmatch-error", default=30, type=float, choices=range(0, 101), metavar="[0-100]",
        help="Minimum phred score for missmatch error at any given variant. Default: %(default)s.")
    arg("--include-pruned", default=False, action="store_true",
        help="Include phasing events marked as pruned by HapCUT2. Default: %(default)s.")
    arg("--discard-untagged", default=False, action="store_true",
        help="Discard alignments not tagged with haplotype info from output. Default: %(default)s.")
    arg("--tag-priority", choices=["READ", "MOLECULE", "CONCORDANT"], default="CONCORDANT", type=str.upper,
        help="How to prioritise between phasing information. 'READ': Read phase >> Molecule phase, 'MOLECULE': "
             "Molecule phase >> Read phase, 'CONCORDANT': Read phase == Molecule phase (agreement between phases). "
             "Default: %(default)s.")
    arg("--min-mapq", type=int, default=0,
        help="Minimum mapping-quality to include reads in analysis Default: %(default)s")

    arg("--molecule-tag", default="MI", help="Molecule SAM tag. Default: %(default)s.")
    arg("--haplotype-tag", default="HP", help="Haplotype SAM tag. Default: %(default)s")
    arg("--phase-set-tag", default="PS", help="Phase set SAM tag. Default: %(default)s")
