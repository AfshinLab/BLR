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
import statistics
import logging

logger = logging.getLogger(__name__)


def main(args):
    logging.info("Starting analysis")
    print(vcf_vcf_error_rate(args.vcf1, args.vcf2, args.indels))
    logging.info("Finished")


def parse_vcf(vcf_file):
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue
            if len(el) != 10:
                logger.error("VCF file must be single-sample.")
                exit(1)
            yield el


def parse_vcf_phase(vcf_file, indels=False):
    PS_index = None
    blocks = defaultdict(list)
    chrom_blocks = dict()

    for el in parse_vcf(vcf_file):
        # get the index where the PS information is
        for i, f in enumerate(el[8].split(':')):
            if i == 0:
                assert (f == 'GT')
            if f == 'PS':
                if not PS_index:
                    PS_index = i
                else:
                    assert (PS_index == i)
                break

    if not PS_index:
        logger.warning("PS flag is missing from VCF. Assuming that all phased variants are in the same phase block.")

    prev_chrom = None
    snp_ix = 0
    for el in parse_vcf(vcf_file):
        consider = True

        phase_data = el[9]
        chrom = el[0]
        a0 = el[3]
        a1 = el[4]
        a2 = None

        if ',' in a1:
            alt_lst = a1.split(',')
            if len(alt_lst) == 2:
                a1, a2 = alt_lst
            else:
                consider = False

        dat = el[9].split(':')
        genotype = dat[0]

        if not (len(genotype) == 3 and genotype[0] in ['0', '1', '2'] and
                genotype[1] in ['|'] and genotype[2] in ['0', '1', '2']):
            consider = False

        if genotype[0] == genotype[2]:
            consider = False

        if consider and not indels and (('0' in genotype and len(a0) != 1) or
                                        ('1' in genotype and len(a1) != 1) or
                                        ('2' in genotype and len(a2) != 1)):
            consider = False

        ps = None
        if not PS_index:
            ps = 1  # put everything in one block
        elif consider and len(dat) > PS_index:
            ps = dat[PS_index]
            if ps == '.':
                consider = False

        if not prev_chrom:
            prev_chrom = chrom

        # If new chromosome, add blocks to chrom_blocks and reset
        if chrom != prev_chrom:
            chrom_blocks[prev_chrom] = [v for k, v in sorted(list(blocks.items())) if len(v) > 1]
            blocks = defaultdict(list)
            snp_ix = 0
            prev_chrom = chrom

        pos = int(el[1]) - 1
        if ps and consider and phase_data[1] == '|':
            blocks[ps].append((snp_ix, pos, phase_data[0:1], phase_data[2:3], a0, a1, a2))

        snp_ix += 1

    # Final
    chrom_blocks[prev_chrom] = [v for k, v in sorted(list(blocks.items())) if len(v) > 1]

    return chrom_blocks


# given a VCF file, simply count the number of heterozygous SNPs present.
def count_snps(vcf_file, ref_name, indels=False):
    count = 0
    for el in parse_vcf(vcf_file):
        if len(el) < 5:
            continue

        if el[0] != ref_name:
            if count == 0:
                continue
            else:
                break

        a0 = el[3]
        a1 = el[4]
        a2 = None

        if ',' in a1:
            alt_lst = a1.split(',')
            if len(alt_lst) == 2:
                a1, a2 = alt_lst
            else:
                continue

        genotype = el[9][:3]

        if not (len(genotype) == 3 and genotype[0] in ['0', '1', '2'] and
                genotype[1] in ['/', '|'] and genotype[2] in ['0', '1', '2']):
            continue

        if genotype[0] == genotype[2]:
            continue

        if (not indels) and (('0' in genotype and len(a0) != 1) or
                             ('1' in genotype and len(a1) != 1) or
                             ('2' in genotype and len(a2) != 1)):
            continue

        count += 1

    return count


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

    for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in hap:
        x = t1_dict[pos]  # base in true haplotype
        y = a1 if allele == 0 else a2  # base in assembled haplotype
        if x == '-' or y == '-':
            if first_SNP:
                continue
            else:
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
        assert (k not in d3)
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

        def create_dict(val, d_type):
            new_dict = defaultdict(d_type)
            if ref and val:
                new_dict[ref] = val
            return new_dict

        self.ref = set()  # set of references in this result (e.g. all chromosomes)
        if ref:
            self.ref.add(ref)

        # these are things that can be summed for the same reference,
        # e.g. switch counts for separate blocks are additive
        self.switch_count = create_dict(switch_count, int)
        self.poss_sw = create_dict(poss_sw, int)
        self.mismatch_count = create_dict(mismatch_count, int)
        self.poss_mm = create_dict(poss_mm, int)
        self.flat_count = create_dict(flat_count, int)
        self.poss_flat = create_dict(poss_flat, int)
        self.phased_count = create_dict(phased_count, int)
        self.AN50_spanlst = create_dict(AN50_spanlst, list)
        self.N50_spanlst = create_dict(N50_spanlst, list)

        # these are things that are non-additive properties, because they
        # refer to the whole reference and would be double-counted
        # e.g. if we combine errors for two blocks, on same chromosome, we add their errors
        # but we can't just add "num_snps", their chromosomes' total snp counts
        # so we use dictionaries to make sure these properties aren't duplicated

        self.num_snps = create_dict(num_snps, int)
        self.maxblk_snps = create_dict(maxblk_snps, int)

        self.switch_loc = create_dict(switch_loc, list)
        self.mismatch_loc = create_dict(mismatch_loc, list)

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
        else:
            return 0

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(mismatch_count) / poss_mm
        else:
            return 0

    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(self.get_switch_count() + self.get_mismatch_count()) / poss_mm
        else:
            return 0

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return float(flat_count) / poss_flat
        else:
            return 0

    def get_AN50(self):
        AN50 = 0
        AN50_spanlst = sum(self.AN50_spanlst.values(), [])
        AN50_spanlst.sort(reverse=True)
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > self.get_num_snps() / 2.0:
                AN50 = span
                break
        return AN50

    def get_N50_phased_portion(self):
        N50 = 0
        N50_spanlst = sum(self.N50_spanlst.values(), [])
        N50_spanlst.sort(reverse=True)

        L = sum(N50_spanlst)

        total = 0
        for span in N50_spanlst:
            total += span
            if total > L / 2.0:
                N50 = span
                break
        return N50

    def get_median_block_length(self):
        spanlst = sum(self.N50_spanlst.values(), [])
        return statistics.median(spanlst)

    def __str__(self):

        s = f"switch rate:        {self.get_switch_rate()}\n" \
            f"mismatch rate:      {self.get_mismatch_rate()}\n" \
            f"flat rate:          {self.get_flat_error_rate()}\n" \
            f"phased count:       {self.get_phased_count()}\n" \
            f"AN50:               {self.get_AN50()}\n" \
            f"N50:                {self.get_N50_phased_portion()}\n" \
            f"num snps max blk:   {sum(self.maxblk_snps.values())}"

        return s


# compute haplotype error rates between 2 VCF files
def vcf_vcf_error_rate(assembled_vcf_file, reference_vcf_file, indels):
    # parse and get stuff to compute error rates
    chrom_a_blocklist = parse_vcf_phase(assembled_vcf_file, indels)
    chrom_t_blocklist = parse_vcf_phase(reference_vcf_file, indels)

    chromosomes = sorted(chrom_a_blocklist)
    err = ErrorResult()
    for c in chromosomes:
        err += error_rate_calc(chrom_t_blocklist[c], chrom_a_blocklist[c], assembled_vcf_file, c, indels)
    return err


def error_rate_calc(t_blocklist, a_blocklist, vcf_file, ref_name, indels=False, phase_set=None):
    num_snps = count_snps(vcf_file, ref_name, indels)

    switch_count = 0
    mismatch_count = 0
    poss_sw = 0  # count of possible positions for switch errors
    poss_mm = 0  # count of possible positions for mismatches
    flat_count = 0
    phased_count = 0
    maxblk_snps = 0
    different_alleles = 0
    switch_loc = []
    mismatch_loc = []
    AN50_spanlst = []
    N50_spanlst = []

    for blk in a_blocklist:

        first_pos = -1
        last_pos = -1
        first_SNP = -1
        last_SNP = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

            if not (a1 == '-' or (phase_set and snp_ix not in phase_set)):

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

        if blk_phased > maxblk_snps:
            maxblk_snps = blk_phased

    for t_block in t_blocklist:
        # convert t_block to a dict for convenience
        t1_dict = defaultdict(lambda: '-')
        t2_dict = defaultdict(lambda: '-')
        a_dict = defaultdict(lambda: ('-', '-', '-'))
        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in t_block:
            t1_dict[pos] = a1
            t2_dict[pos] = a2
            a_dict[pos] = (ref_str, alt1_str, alt2_str)

        # Iterate over SNPs in the true and assembled haplotypes in parallel. i is the index of the current base.
        # x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

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

                    if x == '-' or y == '-' or (phase_set and pos not in phase_set):
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

            if blk_switches[0] < blk_switches[1]:
                switch_count += blk_switches[0]
                mismatch_count += blk_mismatches[0]
                switch_loc += blk_switchlist[0]
                mismatch_loc += blk_mmlist[0]

            else:
                switch_count += blk_switches[1]
                mismatch_count += blk_mismatches[1]
                switch_loc += blk_switchlist[1]
                mismatch_loc += blk_mmlist[1]

        assert len(switch_loc) == switch_count
        assert len(mismatch_loc) == mismatch_count

        # tally up how many possible positions there are for switch errors and mismatches
        # count how many phased SNPs there are so we can calculate a rate of pruned SNPs

        for blk in a_blocklist:
            phased_known = 0
            for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

                if {t1_dict[pos], t2_dict[pos]} != {a1, a2} or (ref_str, alt1_str, alt2_str) != a_dict[pos]:
                    continue

                if t1_dict[pos] != '-' and a1 != '-' and (not phase_set or pos in phase_set):
                    phased_known += 1

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between SNPs counts, and -2 for the two ends.
            if phased_known >= 4:
                poss_sw += (phased_known - 3)
            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                poss_mm += phased_known

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in
        # the assembled haplotype.

        for a_block in a_blocklist:

            flat_count1 = 0
            flat_count2 = 0
            for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in a_block:

                if {t1_dict[pos], t2_dict[pos]} != {a1, a2} or (ref_str, alt1_str, alt2_str) != a_dict[pos]:
                    continue

                if a1 == '-' or a2 == '-' or t1_dict[pos] == '-' or (phase_set and pos not in phase_set):
                    continue

                if a1 != t1_dict[pos]:
                    flat_count1 += 1
                if a2 != t1_dict[pos]:
                    flat_count2 += 1

            if flat_count1 < flat_count2:
                flat_count += flat_count1
            else:
                flat_count += flat_count2

    if different_alleles > 0:
        logger.warning(f"{different_alleles} positions had different ref,alt pairs and were skipped.")

    poss_flat = poss_mm

    if poss_sw == 0 and poss_mm == 0:
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


def add_arguments(parser):
    parser.add_argument('-v1', '--vcf1',
                        help="A phased, single sample VCF file to compute haplotype statistics on.")
    parser.add_argument('-v2', '--vcf2',
                        help="A phased, single sample  VCF file to use as the 'ground truth' haplotype.")
    parser.add_argument('-i', '--indels', action="store_true",
                        help='Use this flag to consider indel variants. Default: %(default)s', default=False)
