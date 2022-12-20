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

from blr.utils import smart_open, tqdm

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
        return 0, int(chromosome)
    elif chromosome == "X":
        return 0, 23
    elif chromosome == "Y":
        return 0, 24
    elif chromosome == "M":
        return 0, 25
    else:
        return 1, chromosome


def parse_variants(records, sample_name, indels=False):
    """Generator for heterozygous variants from VCF records"""
    prev_chrom = None
    variant_index = 0
    for record in records:
        variant_index += 1
        sample = record.samples[sample_name]
        alleles = (record.ref, *record.alts)
        if len(alleles) > 3:
            continue

        genotype = sample["GT"]

        # Skip non-diploid, homozygous genotypes
        if len(genotype) != 2 or len(set(genotype) - {0, 1, 2}) or genotype[0] == genotype[1]:
            continue

        # Skip INDELs if requested
        if not indels and any(len(alleles[g]) != 1 for g in genotype):
            continue

        # Reset index at start of new chromosome
        if record.chrom != prev_chrom and prev_chrom is not None:
            variant_index = 0

        prev_chrom = record.chrom

        yield variant_index, sample, genotype, record, alleles


def get_phaseblocks_chrom(chromosome, vcf_file, sample_name, indels=False):
    """Get chromsome phaseblocks from indexed VCF"""
    phaseset_to_block = defaultdict(list)
    variants_heterozygous = 0
    with VariantFile(vcf_file) as vcf:
        records = vcf.fetch(chromosome)
        for variant_index, sample, genotype, record, alleles in parse_variants(records, sample_name, indels):
            variants_heterozygous += 1

            if not sample.phased:
                continue

            ps = sample.get("PS")
            phaseset_to_block[ps].append((variant_index, record.start, genotype, alleles))

    logger.debug(f"Chromsome {chromosome} has {variants_heterozygous:,} heterozygous variants")

    # Remove block without phaseset if there are multiple blocks
    if None in phaseset_to_block and len(phaseset_to_block) > 1:
        del phaseset_to_block[None]

    blocks = [block for _, block in sorted(list(phaseset_to_block.items())) if len(block) > 1]
    return chromosome, blocks, variants_heterozygous


def get_phaseblocks(vcf_file, sample_name, indels=False):
    """Get phaseblocks for all chromosomes from non-indexed VCF"""
    phaseset_to_block = defaultdict(list)
    chromo_to_blocks = defaultdict(list)
    chrom_to_variants = defaultdict(int)
    prev_chrom = None
    with VariantFile(vcf_file) as vcf:
        for variant_index, sample, genotype, record, alleles in parse_variants(vcf, sample_name, indels):
            chrom_to_variants[record.chrom] += 1

            if record.chrom != prev_chrom and prev_chrom is not None:
                chromo_to_blocks[prev_chrom] = [v for k, v in sorted(list(phaseset_to_block.items())) if len(v) > 1]
                phaseset_to_block.clear()

            prev_chrom = record.chrom

            if not sample.phased:
                continue

            ps = sample.get("PS")
            phaseset_to_block[ps].append((variant_index, record.start, genotype, alleles))

    if phaseset_to_block:
        # Remove block without phaseset if there are multiple blocks
        if None in phaseset_to_block and len(phaseset_to_block) > 1:
            del phaseset_to_block[None]

        chromo_to_blocks[prev_chrom] = [block for _, block in sorted(list(phaseset_to_block.items())) if len(block) > 1]

    return chromo_to_blocks, chrom_to_variants


def parse_vcf_phase(vcf_file, indels=False, chromosomes=None, threads=1):
    with VariantFile(vcf_file) as open_vcf:
        if "PS" not in open_vcf.header.formats:
            logger.warning(f"PS flag is missing from {vcf_file}. Assuming that all phased variants are in the same"
                           " phase block.")

        if len(list(open_vcf.header.samples)) > 1:
            sys.exit("VCF file must be single-sample.")

        sample_name = open_vcf.header.samples[0]

        # Get blocks for all chromosomes if not specified
        chromosomes = chromosomes if chromosomes else open_vcf.header.contigs

    # Don't run in parallel if VCF not indexed
    is_indexed = os.path.isfile(vcf_file + ".tbi") or os.path.isfile(vcf_file + ".cbi")
    if not is_indexed and threads > 1:
        logger.warning(f"Cannot run multiple threads on non-indexed VCF '{vcf_file}'.")

    tqdm_kwargs = dict(total=len(chromosomes), desc="Chromosomes")
    if threads > 1 and is_indexed:
        chrom_blocks = defaultdict(list)
        chrom_to_variants = defaultdict(int)
        func = partial(get_phaseblocks_chrom, vcf_file=vcf_file, sample_name=sample_name, indels=indels)
        with Pool(threads) as workers:
            for chromosome, blocks, nr_het_var in tqdm(workers.imap_unordered(func, chromosomes), **tqdm_kwargs):
                chrom_blocks[chromosome] = blocks
                chrom_to_variants[chromosome] = nr_het_var
    else:
        # If chromosomes are specified and the file is indexed we can fetch the blocks
        # directly for each chromosome for a significant speedup.
        if chromosomes and is_indexed:
            chrom_blocks = defaultdict(list)
            chrom_to_variants = defaultdict(int)
            for chromosome in tqdm(chromosomes, **tqdm_kwargs):
                _,  blocks, nr_het_var = get_phaseblocks_chrom(chromosome, vcf_file, sample_name, indels=indels)
                chrom_blocks[chromosome] = blocks
                chrom_to_variants[chromosome] = nr_het_var
        else:
            chrom_blocks, chrom_to_variants = get_phaseblocks(vcf_file, sample_name=sample_name, indels=indels)

    return chrom_blocks, chrom_to_variants


# this function is needed for "counting ahead" at beginning of blocks.
# error_rate() needs to properly combine switch errors into mismatches
# such that switch errors are minimized (basically, if a block begins
# with an odd number of consecutive switch errors, it should assume
# that this is a sequence of all mismatches and not a "1-less" sequence of mismatches
# with a switch error at the end of it.
def count_consecutive_switches(position_to_genotype_ref, block_asm, allele):
    count = 0
    is_first = True
    switched = False

    for _, pos, genotype, *_ in block_asm:
        x = position_to_genotype_ref[pos][0]  # base in true haplotype
        y = genotype[allele]  # base in assembled haplotype
        if x == '-' or y == '-':
            if is_first:
                continue
            break
        elif is_first:
            switched = (x != y)
            is_first = False
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
    def __init__(
        self,
        ref=None,
        switch_count=None,
        switch_positions=None,
        mismatch_count=None,
        mismatch_positions=None,
        flat_count=None,
        flat_positions=None,
        phased_count=None,
        num_snps=None,
        maxblk_snps=None,
        AN50_spanlst=None,
        N50_spanlst=None,
        switch_loc=None,
        mismatch_loc=None
    ):

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
        self.switch_positions = create_dict(switch_positions, int, ref)
        self.mismatch_count = create_dict(mismatch_count, int, ref)
        self.mismatch_positions = create_dict(mismatch_positions, int, ref)
        self.flat_count = create_dict(flat_count, int, ref)
        self.flat_positions = create_dict(flat_positions, int, ref)
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
        new_err.switch_positions = merge_dicts(self.switch_positions, other.switch_positions)
        new_err.mismatch_count = merge_dicts(self.mismatch_count, other.mismatch_count)
        new_err.mismatch_positions = merge_dicts(self.mismatch_positions, other.mismatch_positions)
        new_err.flat_count = merge_dicts(self.flat_count, other.flat_count)
        new_err.flat_positions = merge_dicts(self.flat_positions, other.flat_positions)
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

    def get_switch_positions(self):
        return sum(self.switch_positions.values())

    def get_mismatch_positions(self):
        return sum(self.mismatch_positions.values())

    def get_flat_positions(self):
        return sum(self.flat_positions.values())

    def get_num_snps(self):
        return sum(self.num_snps.values())

    def get_phased_count(self):
        return sum(self.phased_count.values()) if self.phased_count.values() else "n/a"

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        switch_positions = self.get_switch_positions()
        if switch_positions > 0:
            return float(switch_count) / switch_positions
        return "n/a"

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        mismatch_positions = self.get_mismatch_positions()

        if mismatch_positions > 0:
            return float(mismatch_count) / mismatch_positions
        return "n/a"

    def get_switch_mismatch_rate(self):
        mismatch_positions = self.get_mismatch_positions()

        if mismatch_positions > 0:
            return float(self.get_switch_count() + self.get_mismatch_count()) / mismatch_positions
        return "n/a"

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        flat_positions = self.get_flat_positions()
        if flat_positions > 0:
            return float(flat_count) / flat_positions
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

    def get_auN(self):
        # Calculate auN = 'Area under the Nx curve'
        # see https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
        spans = [value for spanlst in self.N50_spanlst.values() for value in spanlst]
        span_sum = sum(spans)
        span_sq_sum = sum(s**2 for s in spans)
        try:
            return span_sq_sum / span_sum
        except ZeroDivisionError:
            return "n/a"

    def get_median_block_length(self):
        spanlst = [value for spanlst in self.N50_spanlst.values() for value in spanlst]
        return statistics.median(spanlst)

    def get_num_snps_max_blk(self):
        return sum(self.maxblk_snps.values()) if self.maxblk_snps.values() else "n/a"

    def to_txt(self):

        s = f"switch rate:        {self.get_switch_rate()}\n" \
            f"mismatch rate:      {self.get_mismatch_rate()}\n" \
            f"flat rate:          {self.get_flat_error_rate()}\n" \
            f"phased count:       {self.get_phased_count()}\n" \
            f"AN50:               {self.get_AN50()}\n" \
            f"N50:                {self.get_N50_phased_portion()}\n" \
            f"auN:                {self.get_auN()}\n" \
            f"num snps max blk:   {self.get_num_snps_max_blk()}"

        return s

    def to_tsv(self):

        # Header. Sample name is included for MultiQC
        s = "Sample Name\tswitch rate\tmismatch rate\tflat rate\tphased count\tAN50 (Mbp)\tN50 (Mbp)\t" \
            "auN\tnum snps max blk\n"

        values = [
            "",                                           # Need string for MultiQC
            self.get_switch_rate(),
            self.get_mismatch_rate(),
            self.get_flat_error_rate(),
            self.get_phased_count(),
            self.get_AN50()/1_000_000,                    # Show as Mbp
            self.get_N50_phased_portion()/1_000_000,      # Show as Mbp
            self.get_auN() / 1_000_000,                   # Show as Mbp
            sum(self.maxblk_snps.values())
        ]

        s += "\t".join(list(map(str, values))) + "\n"

        return s


# compute haplotype error rates between 2 VCF files
def vcf_vcf_error_rate(assembled_vcf, reference_vcf, keep_indels, input_chromosomes, threads):
    # parse and get stuff to compute error rates
    logger.info(f"Parsing {assembled_vcf}")
    chrom_to_block_asm, chrom_to_variants = parse_vcf_phase(assembled_vcf, keep_indels, input_chromosomes, threads)

    chroms_in_asm = [chrom for chrom, blocks in chrom_to_block_asm.items() if len(blocks) > 0]
    chroms_in_asm.sort(key=chromosome_rank)
    logger.debug(f"Chromsomes in 'vcf1': {','.join(chroms_in_asm)}")

    chromosomes = input_chromosomes if input_chromosomes else chroms_in_asm

    chrom_to_block_ref = defaultdict(list)
    if reference_vcf:
        logger.info(f"Parsing {reference_vcf}")
        chrom_to_block_ref, _ = parse_vcf_phase(reference_vcf, keep_indels, chromosomes, threads)

        chroms_in_ref = [chrom for chrom, blocks in chrom_to_block_ref.items() if len(blocks) > 0]
        chroms_in_ref.sort(key=chromosome_rank)
        logger.debug(f"Chromsomes in 'vcf2': {','.join(chroms_in_ref)}")

    logger.info("Computing statistics")
    chrom_to_err_result = defaultdict(ErrorResult)
    if threads > 1:
        iter_blocks = (
            (chrom_to_block_ref[chrom], chrom_to_block_asm[chrom], chrom, keep_indels, chrom_to_variants[chrom])
            for chrom in chromosomes
        )
        with Pool(threads) as workers:
            for err_result, chromosome in workers.imap_unordered(error_rate_calc_parallel, iter_blocks):
                chrom_to_err_result[chromosome] = err_result
                chrom_to_err_result["all"] += err_result
    else:
        for chrom in tqdm(chromosomes):
            logger.debug(f"Current chromsome = {chrom}")
            chrom_to_err_result[chrom] = error_rate_calc(
                chrom_to_block_ref[chrom], chrom_to_block_asm[chrom], chrom, keep_indels, chrom_to_variants[chrom]
            )
            chrom_to_err_result["all"] += chrom_to_err_result[chrom]
    return chrom_to_err_result, chromosomes


def error_rate_calc_parallel(args):
    return error_rate_calc(*args), args[2]


def error_rate_calc(blocks_ref, blocks_asm, ref_name, indels=False, num_snps=None):
    switch_count = 0
    mismatch_count = 0
    switch_positions = 0  # count of possible positions for switch errors
    mismatch_positions = 0  # count of possible positions for mismatches
    flat_count = 0
    different_alleles = 0
    switch_loc = []
    mismatch_loc = []

    AN50_spanlst, N50_spanlst, maxblk_snps, phased_count = parse_assembled_blocks(blocks_asm)

    for block_ref in blocks_ref:
        # convert block_ref to a dict for convenience
        position_to_alleles_ref, position_to_genotype_ref = mapp_positions_to_block(block_ref)
        block_ref_start = block_ref[0][1]
        block_ref_end = block_ref[-1][1]

        # Iterate over SNPs in the true and assembled haplotypes in parallel. i is the index of the current base.
        # x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for block_asm in blocks_asm:
            block_asm_start = block_asm[0][1]
            block_asm_end = block_asm[-1][1]

            # Skip comparison if block_ref and block_asm are not overlapping
            if block_asm_start > block_ref_end or block_ref_start > block_asm_end:
                continue

            allele_switches = [0, 0]
            allele_mismatches = [0, 0]
            allele_switch_loc = [[], []]
            allele_mismatch_loc = [[], []]

            # choose which allele to score. this only makes a difference for minimizing switch errors vs mismatches
            # in corner cases.
            for a in [0, 1]:
                switched = False
                last_base_was_switch = False
                is_first = True
                for block_asm_index, (variant_index, position, genotype, alleles) in enumerate(block_asm):
                    genotype_ref = position_to_genotype_ref[position]

                    if genotype_ref[0] == '-':
                        continue

                    alleles_ref = position_to_alleles_ref[position]
                    if set(genotype_ref) != set(genotype) or alleles_ref != alleles:
                        if a == 0:
                            different_alleles += 1
                            logger.debug(f"Different alleles for position {ref_name}:{position}")
                            logger.debug(f"Ref {set(genotype_ref)} != {set(genotype)} Query")
                            logger.debug(f"Ref {alleles_ref} != {alleles} Query")
                        continue

                    if is_first:
                        switched = (genotype_ref[0] != genotype[a])
                        consecutive_switches = count_consecutive_switches(
                            position_to_genotype_ref, block_asm[block_asm_index:], a
                        )
                        last_base_was_switch = consecutive_switches % 2 == 1
                        is_first = False
                        continue

                    # if there is a mismatch against the true haplotype and we are in a normal state,
                    # or if there is a "match" that isn't a match because we are in a switched state,
                    # then we need to flip the state again and iterate the count
                    if (genotype_ref[0] != genotype[a] and not switched) or (
                            genotype_ref[0] == genotype[a] and switched):  # current base is mismatched, implying a switch
                        switched = not switched  # flip the "switched" status

                        if last_base_was_switch:  # then this is actually a single-base mismatch
                            # count the 2 switches as a single-base mismatch instead
                            allele_mismatches[a] += 1
                            allele_mismatch_loc[a].append(position)
                            allele_switches[a] -= 1  # undo count from last base switch
                            if len(allele_switch_loc[a]) > 0:
                                allele_switch_loc[a].pop()
                            if allele_switches[a] < 0:
                                allele_switches[a] = 0
                            last_base_was_switch = False

                        else:

                            allele_switches[a] += 1
                            allele_switch_loc[a].append(position)
                            last_base_was_switch = True

                    else:  # current base is not mismatched
                        last_base_was_switch = False

                # special case for switch on last base of previous block_asm; should count as a mismatch
                if last_base_was_switch:
                    # count the switch as a single-base mismatch instead
                    allele_mismatches[a] += 1
                    allele_mismatch_loc[a].append(position)
                    allele_switches[a] -= 1
                    if len(allele_switch_loc[a]) > 0:
                        allele_switch_loc[a].pop()

                    if allele_switches[a] < 0:
                        allele_switches[a] = 0

            # select allele with fewer switches
            i = 0 if allele_switches[0] < allele_switches[1] else 1
            switch_count += allele_switches[i]
            mismatch_count += allele_mismatches[i]
            switch_loc += allele_switch_loc[i]
            mismatch_loc += allele_mismatch_loc[i]

            # tally up how many possible positions there are for switch errors and mismatches
            # count how many phased variants there are so we can calculate a rate of pruned variants
            flat_count_block,  phased_known = get_phased_pos_and_flat_count(
                block_asm, position_to_alleles_ref, position_to_genotype_ref
            )

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between variants counts, and -2 for the two ends.
            if phased_known >= 4:
                switch_positions += (phased_known - 3)

            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                mismatch_positions += phased_known

            flat_count += flat_count_block

        assert len(switch_loc) == switch_count
        assert len(mismatch_loc) == mismatch_count

    if different_alleles > 0:
        logger.warning(f"On {ref_name}: {different_alleles} positions had different ref,alt pairs and were skipped.")

    flat_positions = mismatch_positions

    if blocks_ref and switch_positions == 0 and mismatch_positions == 0:
        logger.warning('Possible switch positions and possible mismatch positions are both 0, it is likely that '
                       'something is very wrong.')

    total_error = ErrorResult(
        ref=ref_name,
        switch_count=switch_count,
        switch_positions=switch_positions,
        mismatch_count=mismatch_count,
        mismatch_positions=mismatch_positions,
        flat_count=flat_count,
        flat_positions=flat_positions,
        phased_count=phased_count,
        num_snps=num_snps,
        maxblk_snps=maxblk_snps,
        AN50_spanlst=AN50_spanlst,
        N50_spanlst=N50_spanlst,
        switch_loc=switch_loc,
        mismatch_loc=mismatch_loc
    )

    return total_error


def get_phased_pos_and_flat_count(block_asm, position_to_alleles_ref, position_to_genotype_ref):
    phased_known = 0
    flat_count1 = 0
    flat_count2 = 0
    for variant_index, pos, genotype, alleles in block_asm:
        alleles_ref = position_to_alleles_ref[pos]
        if set(position_to_genotype_ref[pos]) != set(genotype) or alleles != alleles_ref:
            continue

        if position_to_genotype_ref[pos][0] != '-':
            phased_known += 1

        if position_to_genotype_ref[pos][0] == '-':
            continue

        if genotype[0] != position_to_genotype_ref[pos][0]:
            flat_count1 += 1
        if genotype[1] != position_to_genotype_ref[pos][0]:
            flat_count2 += 1

    # Pick the lowest value
    flat_count = flat_count1 if flat_count1 < flat_count2 else flat_count2
    return flat_count, phased_known


def mapp_positions_to_block(block_ref):
    position_to_genotype = defaultdict(lambda: ('-', '-'))
    position_to_alleles = defaultdict(lambda: ('-', '-', '-'))
    for variant_index, position, genotype, alleles in block_ref:
        position_to_genotype[position] = genotype
        position_to_alleles[position] = alleles
    return position_to_alleles, position_to_genotype


def parse_assembled_blocks(blocks_asm):
    adjusted_block_lengths = []
    block_lengths = []
    variants_phased = 0
    variants_phased_block_max = 0
    for block in blocks_asm:

        position_first = -1
        position_last = -1
        variant_index_first = -1
        variant_index_last = -1
        variants_phased_block = 0

        for variant_index, position, genotype, *_ in block:
            if '-' in genotype:
                continue

            variants_phased += 1

            variants_phased_block += 1
            if position_first == -1:
                position_first = position
                variant_index_first = variant_index
            position_last = position
            variant_index_last = variant_index

        block_index_span = variant_index_last - variant_index_first + 1
        block_length = position_last - position_first

        adjusted_block_lengths.append(
            (block_length * (variants_phased_block / block_index_span), variants_phased_block)
        )
        block_lengths.append(block_length)

        variants_phased_block_max = max(variants_phased_block, variants_phased_block_max)

    return adjusted_block_lengths, block_lengths, variants_phased_block_max, variants_phased


def add_arguments(parser):
    parser.add_argument(
        '-v1', '--vcf1', required=True,
        help="A phased, single sample VCF (uncompressed or bgzip) file to compute haplotype statistics on."
    )
    parser.add_argument(
        '-v2', '--vcf2',
        help="A phased, single sample  VCF (uncompressed or bgzip) file to use as the 'ground truth' haplotype."
    )
    parser.add_argument(
        '-i', '--indels', action="store_true",
        help='Use this flag to consider indel variants. Default: %(default)s.', default=False
    )
    parser.add_argument(
        '--per-chrom', action="store_true", default=False,
        help="Include separate stats for each chromosome. Default: %(default)s."
    )
    parser.add_argument(
        '-c', '--chromosomes',
        help="Name(s) of chromsome(s) to calculate stats for. Multiple chromsomes are joined through commas. "
             "Default: use all chromosomes."
    )
    parser.add_argument(
        "-o", "--output", help="Output file name. Default: Print to stdout."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of threads for reading VCFs. Multithread parsing requires indexed VCFs (.cbi or .tbi). "
             "Default: %(default)s."
    )
