"""
Rules related to phasing of variants (called or reference set)
"""
import os


rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = temp("{base}.calling.unlinked.txt")
    input:
        bam = "{base}.calling.bam",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.calling.unlinked.txt.log"
    params:
        indels = "1" if config["phase_indels"] else "0",
        reference = config["genome_reference"],
    shell:
        "extractHAIRS"
        " --10X 1"
        " --indels {params.indels}"
        # TODO - test to include option `--triallelic` to include GT 1/2 
        " --realign_variants 1"  # Improves overall error-rate
        " --ref {params.reference}"
        " --bam {input.bam}"
        " --VCF {input.vcf}"
        " --out {output.unlinked} 2> {log}"


rule extract_long_read_chunk:
    output:
        bam = temp("chunks/{chunk}.long_read.bam")
    input:
        bed = "chunks/{chunk}.bed"
    params:
        bam = config["long_read_bam"],
    shell:
        "samtools view"
        " -bh"
        " -q 20" # Only MAPQ >= 20 reads are used for phasing
        " -o {output.bam}"
        " -M -L {input.bed}"
        " {params.bam}"


# Map long_read_type to required extractHAIRS flag
long_read_flags = {"pacbio": "--pacbio 1", "ont": "--ont 1"}


rule hapcut2_extracthairs_long:
    """Extract heterozygous variants covered by long read alignments in BAM"""
    output:
        txt = temp("{base}.calling.linked_long.txt")
    input:
        bam = "{base}.long_read.bam",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.calling.linked_long.txt.log"
    params:
        long_read_flag = long_read_flags.get(config["long_read_type"], ""),
        reference = config["genome_reference"],
    shell:
        "extractHAIRS"
        " --indels 0"
        " {params.long_read_flag}"
        " --ep 1"  # Estimate parameters 
        " --new_format 1"
        " --ref {params.reference}"
        " --bam {input.bam}"
        " --VCF {input.vcf}"
        " --out {output.txt}"
        " 2> {log}"


rule hapcut2_linkfragments:
    """Link heterozygous variants together using barcode information"""
    output:
        linked = temp("{base}.calling.linked.txt")
    input:
        bam = "{base}.calling.bam",
        bai = "{base}.calling.bam.bai",
        vcf = "{base}.phaseinput.vcf",
        unlinked = "{base}.calling.unlinked.txt"
    log: "{base}.calling.linked.txt.log"
    params:
        window = config["window_size"],
    shell:
        "LinkFragments.py"
        " --bam {input.bam}"
        " -v {input.vcf}"
        " --fragments {input.unlinked}"
        " --out {output.linked}"
        " --distance {params.window} &> {log}"


rule merge_fragments:
    """Merge fragments from linked and long reads."""
    output: 
        txt = "{base}.calling.linked_merged.txt"
    input:
        txt1 = "{base}.calling.linked.txt",
        txt2 = "{base}.calling.linked_long.txt"
    shell:
        "cat {input.txt1} {input.txt2} > {output.txt}"


def get_fragments(wildcards):
    if config["long_read_bam"] is not None:
        return "{base}.calling.linked_merged.txt"
    return "{base}.calling.linked.txt"


rule hapcut2_phasing:
    """Phase heterozygous variants using HapCUT2. Output phased VCF"""
    output:
        phase = "{base}.calling.phase",
        phased_vcf = "{base}.calling.phased.vcf.gz"
    input:
        txt = get_fragments,
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.calling.phase.log"
    shell:
        "hapcut2"
        " --nf 1"
        " --fragments {input.txt}"
        " --vcf {input.vcf}"
        " --out {output.phase}"
        " --error_analysis_mode 1"
        " --outvcf 1 &> {log}"
        " && "
        # Compress phased VCF
        "bgzip -c {output.phase}.phased.VCF > {output.phased_vcf}"
        " && "
        # Remove uncompressed VCF
        " rm -f {output.phase}.phased.VCF"


rule hapcut2_stats:
    """Collect phasing statistics."""
    output:
        stats = "final.phasing_stats.txt",
        plot_stats = "final.phasing_stats.plots.txt"
    input:
        vcf1 = "final.phased.vcf.gz",
        vcf1_index = "final.phased.vcf.gz.tbi",
    params:
        vcf2 = f" -v2 {config['phasing_ground_truth']}" if config['phasing_ground_truth'] else "",
        indels = " --indels" if config["phase_indels"] else "",
        reference_lengths = config["genome_reference"] + ".fai",
    log: "final.phasing_stats.txt.log"
    threads: 20
    shell:
        "blr calculate_haplotype_statistics"
        " -v1 {input.vcf1}"
        " {params.vcf2}"
        " {params.indels}"
        " --per-chrom"
        " --reference-lengths {params.reference_lengths}"
        " --stats {output.plot_stats}"
        " --threads {threads}"
        " -o {output.stats} 2> {log}"


rule haplotag:
    """
    Transfer haplotype information from the phased VCF file to the bam file.
    Adds HP tag with haplotype (1 or 2) and PS tag with phase set information.
    """
    output:
        bam = "{base}.calling.phased.bam"
    input:
        bam = "{base}.calling.bam",
        vcf = "{base}.calling.phased.vcf.gz",
        vcf_index = "{base}.calling.phased.vcf.gz.tbi"
    log: "{base}.calling.phased.bam.log"
    params:
        ignore_readgroups = "--ignore-read-groups" if config["reference_variants"] else "",
        window = config["window_size"],
        reference = config["genome_reference"],
    shell:
        "whatshap haplotag"
        " {input.vcf}"
        " {input.bam}"
        " --linked-read-distance-cutoff {params.window}"
        " --reference {params.reference}"
        " {params.ignore_readgroups}"
        " 2> {log}"
        " |"
        " samtools view -bh -o {output.bam} -"  # Faster than using whatshap for compression



