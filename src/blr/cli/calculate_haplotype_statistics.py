"""
Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.

Based on script: calculate_haplotype_statistisc.py form HapCUT2
https://github.com/vibansal/HapCUT2/blob/master/utilities/calculate_haplotype_statistics.py
"""

# Copyright (c) 2018, Bansal Lab
# All rights reserved.
# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

from collections import defaultdict
from functools import partial
import logging
from multiprocessing import Pool
import os
import statistics
import sys

from pysam import VariantFile

from blr.utils import smart_open

logger = logging.getLogger(__name__)


def main(args):
    logger.info("Starting analysis")
    if not args.vcf2:
        logger.info("No reference vcf provided - error rates will not be computed.")

    chromosomes = args.chromosomes.split(",") if args.chromosomes else None

    stats, chromosomes = vcf_vcf_error_rate(args.vcf1, args.vcf2, args.indels, chromosomes, args.threads)

    with smart_open(args.output) as file:
        if args.per_chrom:
            for c in chromosomes:
                print(f"----------- {c} -----------", file=file)
                print(stats[c].to_txt(), file=file)
            print("----------- All -----------", file=file)

        print(stats["all"].to_txt(), file=file)

    logger.info("Finished")


def chromosome_rank(chromosome: str):
    chromosome = chromosome.replace("chr", "")
    if chromosome.isnumeric():
        return (0, int(chromosome))
    elif chromosome == "X":
        return (0, 23)
    elif chromosome == "Y":
        return (0, 24)
    elif chromosome == "M":
        return (0, 25)
    else:
        return (1, chromosome)


def parse_variants(vcf, sample_name, indels=False):
    """Generator for heterozygous varinats from VCF records"""
    prev_chrom = None
    snp_ix = 0
    for rec in vcf:
        snp_ix += 1
        sample = rec.samples[sample_name]
        a0 = rec.ref
        a1, *a2 = rec.alts

        if len(a2) == 0:
            a2 = None
        elif len(a2) == 1:
            a2 = a2[0]
        else:
            continue

        genotype = sample["GT"]

        if len(genotype) != 2 or len(set(genotype) - {0, 1, 2}) or genotype[0] == genotype[1]:
            continue

        if not indels and any(len([a0, a1, a2][g]) != 1 for g in genotype):
            continue

        if rec.chrom != prev_chrom and prev_chrom is not None:
            snp_ix = 0

        prev_chrom = rec.chrom

        yield snp_ix, sample, genotype, rec, a0, a1, a2


def get_phaseblocks_chrom(chromosome, vcf_file, sample_name, indels=False):
    """Get chromsome phaseblocks from indexed VCF"""
    blocks = defaultdict(list)
    nr_het_var = 0
    with VariantFile(vcf_file) as vcf:
        for snp_ix, sample, genotype, rec, a0, a1, a2 in parse_variants(vcf.fetch(chromosome), sample_name, indels):
            nr_het_var += 1

            if not sample.phased:
                continue

            ps = sample.get("PS")
            if ps is not None:
                blocks[ps].append((snp_ix, rec.start, genotype[0], genotype[1], a0, a1, a2))

    logger.debug(f"Chromsome {chromosome} has {nr_het_var:,} heterozygous variants")
    blocklist = [v for k, v in sorted(list(blocks.items())) if len(v) > 1]
    return chromosome, blocklist, nr_het_var


def get_phaseblocks(vcf_file, sample_name, indels=False):
    """Get phaseblocks for all chromosomes from non-indexed VCF"""
    blocks = defaultdict(list)
    chrom_blocks = defaultdict(list)
    nr_het_var_per_chrom = defaultdict(int)
    prev_chrom = None
    with VariantFile(vcf_file) as vcf:
        for snp_ix, sample, genotype, rec, a0, a1, a2 in parse_variants(vcf, sample_name, indels):
            nr_het_var_per_chrom[rec.chrom] += 1

            if rec.chrom != prev_chrom and prev_chrom is not None:
                chrom_blocks[prev_chrom] = [v for k, v in sorted(list(blocks.items())) if len(v) > 1]
                blocks.clear()

            prev_chrom = rec.chrom

            if not sample.phased:
                continue

            ps = sample.get("PS")

            if ps is not None:
                blocks[ps].append((snp_ix, rec.start, genotype[0], genotype[1], a0, a1, a2))

    if blocks:
        chrom_blocks[prev_chrom] = [v for k, v in sorted(list(blocks.items())) if len(v) > 1]

    return chrom_blocks, nr_het_var_per_chrom


def parse_vcf_phase(vcf_file, indels=False, chromosomes=None, threads=1):
    with VariantFile(vcf_file) as open_vcf:
        if "PS" not in open_vcf.header.formats:
            logger.warning("PS flag is missing from VCF. Assuming that all phased variants are in the same phase "
                           "block.")

        if len(list(open_vcf.header.samples)) > 1:
            sys.exit("VCF file must be single-sample.")

        sample_name = open_vcf.header.samples[0]

        # Get blocks for all chromosomes if not specified
        chromosomes = chromosomes if chromosomes else open_vcf.header.contigs

    # Don't run in parallel if VCF not indexed
    is_indexed = os.path.isfile(vcf_file + ".tbi") or os.path.isfile(vcf_file + ".cbi")
    if not is_indexed and threads > 1:
        logger.warning(f"Cannot run multiple threads on non-indexed VCF '{vcf_file}'.")
        threads = 1

    if threads > 1:
        chrom_blocks = defaultdict(list)
        nr_het_var_per_chrom = defaultdict(int)
        func = partial(get_phaseblocks_chrom, vcf_file=vcf_file, sample_name=sample_name, indels=indels)
        with Pool(threads) as workers:
            for chromosome, blocks, nr_het_var in workers.imap_unordered(func, chromosomes):
                chrom_blocks[chromosome] = blocks
                nr_het_var_per_chrom[chromosome] = nr_het_var
    else:
        # If chromosomes are specified and the file is indexed we can fetch the blocks
        # directly for each chromosome for a significant speedup.
        if chromosomes and is_indexed:
            chrom_blocks = defaultdict(list)
            nr_het_var_per_chrom = defaultdict(int)
            for chromosome in chromosomes:
                _,  blocks, nr_het_var = get_phaseblocks_chrom(chromosome, vcf_file, sample_name, indels=indels)
                chrom_blocks[chromosome] = blocks
                nr_het_var_per_chrom[chromosome] = nr_het_var
        else:
            chrom_blocks, nr_het_var_per_chrom = get_phaseblocks(vcf_file, sample_name=sample_name, indels=indels)

    return chrom_blocks, nr_het_var_per_chrom


# this function is needed for "counting ahead" at beginning of blocks.
# error_rate() needs to properly combine switch errors into mismatches
# such that switch errors are minimized (basically, if a block begins
# with an odd number of consecutive switch errors, it should assume
# that this is a sequence of all mismatches and not a "1-less" sequence of mismatches
# with a switch error at the end of it.
def count_consecutive_switches(t1_dict, hap, allele):
    count = 0
    first_SNP = True
    switched = False

    for _, pos, a1, a2, _, _, _ in hap:
        x = t1_dict[pos]  # base in true haplotype
        y = a1 if allele == 0 else a2  # base in assembled haplotype
        if x == '-' or y == '-':
            if first_SNP:
                continue
            break
        elif first_SNP:
            switched = (t1_dict[pos] != y)
            first_SNP = False
        elif (x != y and not switched) or (x == y and switched):
            count += 1
            switched = not switched
        else:
            break
    return count


# combine two dicts
def merge_dicts(d1, d2):
    d3 = d2.copy()
    for k, v in d1.items():
        assert k not in d3
        d3[k] = v
    return d3


# the "ErrorResult" abstraction and its overloaded addition operator are handy
# for combining results for the same chromosome across blocks (when the "ground truth"
# is a set of blocks rather than trio), and combining results across different chromosomes
# into genome wide stats
class ErrorResult:
    def __init__(self, ref=None, switch_count=None, poss_sw=None, mismatch_count=None, poss_mm=None, flat_count=None,
                 poss_flat=None, phased_count=None, num_snps=None, maxblk_snps=None, AN50_spanlst=None,
                 N50_spanlst=None, switch_loc=None, mismatch_loc=None):

        def create_dict(val, d_type, ref):
            new_dict = defaultdict(d_type)
            if ref is not None and val is not None:
                new_dict[ref] = val
            return new_dict

        # set of references in this result (e.g. all chromosomes)
        self.ref = {ref} if ref is not None else set()

        # these are things that can be summed for the same reference,
        # e.g. switch counts for separate blocks are additive
        self.switch_count = create_dict(switch_count, int, ref)
        self.poss_sw = create_dict(poss_sw, int, ref)
        self.mismatch_count = create_dict(mismatch_count, int, ref)
        self.poss_mm = create_dict(poss_mm, int, ref)
        self.flat_count = create_dict(flat_count, int, ref)
        self.poss_flat = create_dict(poss_flat, int, ref)
        self.phased_count = create_dict(phased_count, int, ref)
        self.AN50_spanlst = create_dict(AN50_spanlst, list, ref)
        self.N50_spanlst = create_dict(N50_spanlst, list, ref)

        # these are things that are non-additive properties, because they
        # refer to the whole reference and would be double-counted
        # e.g. if we combine errors for two blocks, on same chromosome, we add their errors
        # but we can't just add "num_snps", their chromosomes' total snp counts
        # so we use dictionaries to make sure these properties aren't duplicated

        self.num_snps = create_dict(num_snps, int, ref)
        self.maxblk_snps = create_dict(maxblk_snps, int, ref)

        self.switch_loc = create_dict(switch_loc, list, ref)
        self.mismatch_loc = create_dict(mismatch_loc, list, ref)

    # combine two error rate results
    def __add__(self, other):
        new_err = ErrorResult()

        new_err.ref = self.ref.union(other.ref)
        new_err.switch_count = merge_dicts(self.switch_count, other.switch_count)
        new_err.poss_sw = merge_dicts(self.poss_sw, other.poss_sw)
        new_err.mismatch_count = merge_dicts(self.mismatch_count, other.mismatch_count)
        new_err.poss_mm = merge_dicts(self.poss_mm, other.poss_mm)
        new_err.flat_count = merge_dicts(self.flat_count, other.flat_count)
        new_err.poss_flat = merge_dicts(self.poss_flat, other.poss_flat)
        new_err.phased_count = merge_dicts(self.phased_count, other.phased_count)
        new_err.AN50_spanlst = merge_dicts(self.AN50_spanlst, other.AN50_spanlst)
        new_err.N50_spanlst = merge_dicts(self.N50_spanlst, other.N50_spanlst)
        new_err.num_snps = merge_dicts(self.num_snps, other.num_snps)
        new_err.maxblk_snps = merge_dicts(self.maxblk_snps, other.maxblk_snps)
        new_err.switch_loc = merge_dicts(self.switch_loc, other.switch_loc)
        new_err.mismatch_loc = merge_dicts(self.mismatch_loc, other.mismatch_loc)

        return new_err

    def get_switch_count(self):
        return sum(self.switch_count.values())

    def get_mismatch_count(self):
        return sum(self.mismatch_count.values())

    def get_flat_count(self):
        return sum(self.flat_count.values())

    def get_poss_sw(self):
        return sum(self.poss_sw.values())

    def get_poss_mm(self):
        return sum(self.poss_mm.values())

    def get_poss_flat(self):
        return sum(self.poss_flat.values())

    def get_num_snps(self):
        return sum(self.num_snps.values())

    def get_phased_count(self):
        return sum(self.phased_count.values())

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        poss_sw = self.get_poss_sw()
        if poss_sw > 0:
            return float(switch_count) / poss_sw
        return "n/a"

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(mismatch_count) / poss_mm
        return "n/a"

    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(self.get_switch_count() + self.get_mismatch_count()) / poss_mm
        return "n/a"

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return float(flat_count) / poss_flat
        return "n/a"

    def get_AN50(self):
        AN50_spanlst = [value for spanlst in self.AN50_spanlst.values() for value in spanlst]
        AN50_spanlst.sort(reverse=True)
        half_num_snps = self.get_num_snps() / 2.0
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > half_num_snps:
                return span
        return "n/a"

    def get_N50_phased_portion(self):
        N50_spanlst = [value for spanlst in self.N50_spanlst.values() for value in spanlst]
        N50_spanlst.sort(reverse=True)

        half_L = sum(N50_spanlst) / 2.0

        total = 0
        for span in N50_spanlst:
            total += span
            if total > half_L:
                return span
        return "n/a"

    def get_median_block_length(self):
        spanlst = [value for spanlst in self.N50_spanlst.values() for value in spanlst]
        return statistics.median(spanlst)

    def to_txt(self):

        s = f"switch rate:        {self.get_switch_rate()}\n" \
            f"mismatch rate:      {self.get_mismatch_rate()}\n" \
            f"flat rate:          {self.get_flat_error_rate()}\n" \
            f"phased count:       {self.get_phased_count()}\n" \
            f"AN50:               {self.get_AN50()}\n" \
            f"N50:                {self.get_N50_phased_portion()}\n" \
            f"num snps max blk:   {sum(self.maxblk_snps.values())}"

        return s

    def to_tsv(self):

        # Header. Sample name is included for MultiQC
        s = "Sample Name\tswitch rate\tmismatch rate\tflat rate\tphased count\tAN50 (Mbp)\tN50 (Mbp)\t" \
            "num snps max blk\n"

        values = [
            "",                                           # Need string for MultiQC
            self.get_switch_rate(),
            self.get_mismatch_rate(),
            self.get_flat_error_rate(),
            self.get_phased_count(),
            self.get_AN50()/1_000_000,                    # Show as Mbp
            self.get_N50_phased_portion()/1_000_000,      # Show as Mbp
            sum(self.maxblk_snps.values())
        ]

        s += "\t".join(list(map(str, values))) + "\n"

        return s


# compute haplotype error rates between 2 VCF files
def vcf_vcf_error_rate(assembled_vcf_file, reference_vcf_file, indels, input_chromosomes, threads):
    # parse and get stuff to compute error rates
    logger.info(f"Parsing {assembled_vcf_file}")
    chrom_a_blocklist, nr_het_var = parse_vcf_phase(assembled_vcf_file, indels, input_chromosomes, threads)
    chroms_in_a = [chrom for chrom, blocks in chrom_a_blocklist.items() if len(blocks) > 0]
    chroms_in_a.sort(key=chromosome_rank)
    logger.debug(f"Chromsomes in 'vcf1': {','.join(chroms_in_a)}")

    chromosomes = input_chromosomes if input_chromosomes else chroms_in_a

    if reference_vcf_file:
        logger.info(f"Parsing {reference_vcf_file}")
        chrom_t_blocklist, _ = parse_vcf_phase(reference_vcf_file, indels, chromosomes, threads)
        chroms_in_t = [chrom for chrom, blocks in chrom_t_blocklist.items() if len(blocks) > 0]
        chroms_in_t.sort(key=chromosome_rank)
        logger.debug(f"Chromsomes in 'vcf2': {','.join(chroms_in_t)}")
    else:
        chrom_t_blocklist = defaultdict(list)

    logger.info("Computing statistics")
    err = defaultdict(ErrorResult)
    for c in chromosomes:
        logger.debug(f"Current chromsome = {c}")
        err[c] = error_rate_calc(chrom_t_blocklist[c], chrom_a_blocklist[c], c, indels, num_snps=nr_het_var[c])
        err["all"] += err[c]
    return err, chromosomes


def error_rate_calc(t_blocklist, a_blocklist, ref_name, indels=False, num_snps=None):
    switch_count = 0
    mismatch_count = 0
    poss_sw = 0  # count of possible positions for switch errors
    poss_mm = 0  # count of possible positions for mismatches
    flat_count = 0
    different_alleles = 0
    switch_loc = []
    mismatch_loc = []

    AN50_spanlst, N50_spanlst, maxblk_snps, phased_count = parse_assembled_blocks(a_blocklist)

    for t_block in t_blocklist:
        # convert t_block to a dict for convenience
        a_dict, t1_dict, t2_dict = mapp_positions_to_block(t_block)
        t_block_start = t_block[0][1]
        t_block_end = t_block[-1][1]

        # Iterate over SNPs in the true and assembled haplotypes in parallel. i is the index of the current base.
        # x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:
            a_block_start = a_block[0][1]
            a_block_end = a_block[-1][1]

            # Skip comparison if t_block and a_block are not overlapping
            if a_block_start > t_block_end or t_block_start > a_block_end:
                continue

            blk_switches = [0, 0]
            blk_mismatches = [0, 0]
            blk_switchlist = [[], []]
            blk_mmlist = [[], []]

            # choose which allele to score. this only makes a difference for minimizing switch errors vs mismatches
            # in corner cases.
            for a in [0, 1]:
                switched = False
                last_base_was_switch = False
                first_SNP = True
                for blk_ix, (snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str) in enumerate(a_block):
                    y = a1 if a == 0 else a2
                    x = t1_dict[pos]

                    if x == '-' or y == '-':
                        continue

                    if {t1_dict[pos], t2_dict[pos]} != {a1, a2} or (ref_str, alt1_str, alt2_str) != a_dict[pos]:
                        if a == 0:
                            different_alleles += 1
                        continue

                    if first_SNP:
                        switched = (x != y)
                        if count_consecutive_switches(t1_dict, a_block[blk_ix:], a) % 2 == 1:
                            last_base_was_switch = True
                        else:
                            last_base_was_switch = False
                        first_SNP = False
                        continue

                    # if there is a mismatch against the true haplotype and we are in a normal state,
                    # or if there is a "match" that isn't a match because we are in a switched state,
                    # then we need to flip the state again and iterate the count
                    if (x != y and not switched) or (
                            x == y and switched):  # current base is mismatched, implying a switch
                        switched = not switched  # flip the "switched" status

                        if last_base_was_switch:  # then this is actually a single-base mismatch
                            # count the 2 switches as a single-base mismatch instead
                            blk_mismatches[a] += 1
                            blk_mmlist[a].append(pos)
                            blk_switches[a] -= 1  # undo count from last base switch
                            if len(blk_switchlist[a]) > 0:
                                blk_switchlist[a].pop()
                            if blk_switches[a] < 0:
                                blk_switches[a] = 0
                            last_base_was_switch = False

                        else:

                            blk_switches[a] += 1
                            blk_switchlist[a].append(pos)
                            last_base_was_switch = True

                    else:  # current base is not mismatched
                        last_base_was_switch = False

                # special case for switch on last base of previous a_block; should count as a mismatch
                if last_base_was_switch:
                    # count the switch as a single-base mismatch instead
                    blk_mismatches[a] += 1
                    blk_mmlist[a].append(pos)
                    blk_switches[a] -= 1
                    if len(blk_switchlist[a]) > 0:
                        blk_switchlist[a].pop()

                    if blk_switches[a] < 0:
                        blk_switches[a] = 0

            i = 0 if blk_switches[0] < blk_switches[1] else 1
            switch_count += blk_switches[i]
            mismatch_count += blk_mismatches[i]
            switch_loc += blk_switchlist[i]
            mismatch_loc += blk_mmlist[i]

            # tally up how many possible positions there are for switch errors and mismatches
            # count how many phased SNPs there are so we can calculate a rate of pruned SNPs
            flat_count1, flat_count2, phased_known = get_phased_pos_and_flat_count(a_block, a_dict, t1_dict, t2_dict)

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between SNPs counts, and -2 for the two ends.
            if phased_known >= 4:
                poss_sw += (phased_known - 3)
            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                poss_mm += phased_known

            flat_count += flat_count1 if flat_count1 < flat_count2 else flat_count2

        assert len(switch_loc) == switch_count
        assert len(mismatch_loc) == mismatch_count

    if different_alleles > 0:
        logger.warning(f"{different_alleles} positions had different ref,alt pairs and were skipped.")

    poss_flat = poss_mm

    if t_blocklist and poss_sw == 0 and poss_mm == 0:
        logger.warning('Possible switch positions and possible mismatch positions are both 0, it is likely that '
                       'something is very wrong.')

    total_error = ErrorResult(
        ref=ref_name,
        switch_count=switch_count,
        poss_sw=poss_sw,
        mismatch_count=mismatch_count,
        poss_mm=poss_mm,
        flat_count=flat_count,
        poss_flat=poss_flat,
        phased_count=phased_count,
        num_snps=num_snps,
        maxblk_snps=maxblk_snps,
        AN50_spanlst=AN50_spanlst,
        N50_spanlst=N50_spanlst,
        switch_loc=switch_loc,
        mismatch_loc=mismatch_loc
    )

    return total_error


def get_phased_pos_and_flat_count(a_block, a_dict, t1_dict, t2_dict):
    phased_known = 0
    flat_count1 = 0
    flat_count2 = 0
    for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in a_block:

        if {t1_dict[pos], t2_dict[pos]} != {a1, a2} or (ref_str, alt1_str, alt2_str) != a_dict[pos]:
            continue

        if t1_dict[pos] != '-' and a1 != '-':
            phased_known += 1

        if a1 == '-' or a2 == '-' or t1_dict[pos] == '-':
            continue

        if a1 != t1_dict[pos]:
            flat_count1 += 1
        if a2 != t1_dict[pos]:
            flat_count2 += 1
    return flat_count1, flat_count2, phased_known


def mapp_positions_to_block(t_block):
    t1_dict = defaultdict(lambda: '-')
    t2_dict = defaultdict(lambda: '-')
    a_dict = defaultdict(lambda: ('-', '-', '-'))
    for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in t_block:
        t1_dict[pos] = a1
        t2_dict[pos] = a2
        a_dict[pos] = (ref_str, alt1_str, alt2_str)
    return a_dict, t1_dict, t2_dict


def parse_assembled_blocks(a_blocklist):
    AN50_spanlst = []
    N50_spanlst = []
    phased_count = 0
    maxblk_snps = 0
    for blk in a_blocklist:

        first_pos = -1
        last_pos = -1
        first_SNP = -1
        last_SNP = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

            if a1 != '-':

                phased_count += 1

                blk_phased += 1
                if first_pos == -1:
                    first_pos = pos
                    first_SNP = snp_ix
                last_pos = pos
                last_SNP = snp_ix

        blk_total = last_SNP - first_SNP + 1

        AN50_spanlst.append(((last_pos - first_pos) * (float(blk_phased) / blk_total), blk_phased))
        N50_spanlst.append((last_pos - first_pos))

        maxblk_snps = max(blk_phased, maxblk_snps)

    return AN50_spanlst, N50_spanlst, maxblk_snps, phased_count


def add_arguments(parser):
    parser.add_argument('-v1', '--vcf1',
                        help="A phased, single sample VCF (uncompressed or bgzip) file to compute haplotype "
                             "statistics on.")
    parser.add_argument('-v2', '--vcf2',
                        help="A phased, single sample  VCF (uncompressed or bgzip) file to use as the 'ground truth' "
                             "haplotype.")
    parser.add_argument('-i', '--indels', action="store_true",
                        help='Use this flag to consider indel variants. Default: %(default)s', default=False)
    parser.add_argument('--per-chrom', action="store_true", default=False,
                        help="Include separate stats for each chromosome. Default: %(default)s")
    parser.add_argument('-c', '--chromosomes',
                        help="Name(s) of chromsome(s) to calculate stats for. Multiple chromsomes are joined through"
                             " commas. Default: use all chromosomes")
    parser.add_argument("-o", "--output", help="Output file name. Default: Print to stdout.")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for reading VCFs. Multithread parsing requires indexed VCFs (.cbi or "
                             ".tbi). Default: %(default)s.")
