"""
Rules related to phasing of variants (called or reference set)
"""
import os


def get_linked_vcf(wildcards):
    """Include linked file in input until symlink is no longer used"""
    if config["reference_variants"]:
        return []
    else:
        if config["filter_variants"]:
            return f"{wildcards.base}.variants.called.filtered.vcf"
        else:
            return f"{wildcards.base}.variants.called.vcf"


rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = temp("{base}.calling.unlinked.txt")
    input:
        bam = "{base}.calling.bam",
        vcf = "{base}.phaseinput.vcf",
        vcf_link = get_linked_vcf
    log: "{base}.calling.unlinked.txt.log"
    params:
        indels = "1" if config["phase_indels"] else "0"
    shell:
        "extractHAIRS"
        " --10X 1"
        " --indels {params.indels}"
        " --realign_variants 1"  # Improves overall error-rate
        " --ref {config[genome_reference]}"
        " --bam {input.bam}"
        " --VCF {input.vcf}"
        " --out {output.unlinked} 2> {log}"


rule hapcut2_linkfragments:
    """Link heterozygous variants together using barcode information"""
    output:
        linked = temp("{base}.calling.linked.txt")
    input:
        bam = "{base}.calling.bam",
        bai = "{base}.calling.bam.bai",
        vcf = "{base}.phaseinput.vcf",
        vcf_link = get_linked_vcf,
        unlinked = "{base}.calling.unlinked.txt"
    log: "{base}.calling.linked.txt.log"
    shell:
        "LinkFragments.py"
        " --bam {input.bam}"
        " -v {input.vcf}"
        " --fragments {input.unlinked}"
        " --out {output.linked}"
        " --distance {config[window_size]} &> {log}"


rule hapcut2_phasing:
    """Phase heterozygous variants using HapCUT2. Output phased VCF"""
    output:
        phase = "{base}.calling.phase",
        phased_vcf = temporary("{base}.calling.phased.vcf")
    input:
        linked = "{base}.calling.linked.txt",
        vcf = "{base}.phaseinput.vcf",
        vcf_link = get_linked_vcf
    log: "{base}.calling.phase.log"
    shell:
        "hapcut2"
        " --nf 1"
        " --fragments {input.linked}"
        " --vcf {input.vcf}"
        " --out {output.phase}"
        " --error_analysis_mode 1"
        " --outvcf 1 &> {log}"
        " && "
        "mv {output.phase}.phased.VCF {output.phased_vcf}"


rule hapcut2_stats:
    """Get phasing statistics relative the ground truth. See https://github.com/vibansal/HapCUT2/tree/master/utilities
    for details. """
    output:
        stats = "final.phasing_stats.txt"
    input:
        vcf1 = "final.phased.vcf.gz",
        vcf1_index = "final.phased.vcf.gz.tbi",
    params:
        vcf2 = f" -v2 {config['phasing_ground_truth']}" if config['phasing_ground_truth'] else "",
        indels = " --indels" if config["phase_indels"] else ""
    log: "final.phasing_stats.txt.log"
    threads: 20
    shell:
        "blr calculate_haplotype_statistics"
        " -v1 {input.vcf1}"
        " {params.vcf2}"
        " {params.indels}"
        " --per-chrom"
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
        ignore_readgroups = "--ignore-read-groups" if config["reference_variants"] else ""
    shell:
        "whatshap haplotag"
        " {input.vcf}"
        " {input.bam}"
        " --linked-read-distance-cutoff {config[window_size]}"
        " --reference {config[genome_reference]}"
        " -o {output.bam}"
        " {params.ignore_readgroups}"
        " 2> {log}"



def bams_for_lsv_calling(wildcards):
    phased = ".phased" if wildcards.chunk in {c[0].name for c in chunks["phased"]} else ""
    return {
        "bam": f"chunks/{{chunk}}.calling{phased}.bam",
        "bai": f"chunks/{{chunk}}.calling{phased}.bam.bai"
    }


rule build_config:
    """
    Builds a config file required for running NAIBR.
    """
    output:
        config = "chunks/{chunk}.naibr.config"
    input: unpack(bams_for_lsv_calling)
    log: "chunks/{chunk}.naibr.config.log"
    params:
        cwd = os.getcwd(),
        blacklist = f"--blacklist {config['naibr_blacklist']}" if config['naibr_blacklist'] else ""
    shell:
        "blr naibrconfig"
        " --bam-file {input.bam}"
        " --outdir {params.cwd}/chunks/{wildcards.chunk}_naibr"
        " --distance 10000"
        " --min-mapq {config[naibr_min_mapq]}"
        " --min-sv 1000"
        " --threads 1"
        " --min-overlaps 3"
        " {params.blacklist}"
        " --output {output.config}"
        " 2> {log}"


rule get_naibr_path:
    """Symlink existing NAIBR repo if specified or clone the latest version"""
    output:
        naibr_path = temp(directory("NAIBR"))
    run:
        if config["naibr_path"]:
            shell("ln -s {config[naibr_path]} {output.naibr_path}")
        else:
            shell("git clone https://github.com/raphael-group/NAIBR.git {output.naibr_path}")


rule lsv_calling:
    """
    Runs NAIBR for LSV calling. This involves activating a python2 env, changing wd, running and changing back wd and
    env.
    """
    output:
        results = "{base}.naibr_sv_calls.tsv"
    input:
        config = "{base}.naibr.config",
        naibr_path = "NAIBR"
    log: "{base}.naibr_sv_calls.tsv.log"
    threads: 2
    conda: "../naibr-environment.yml"
    params:
        cwd = os.getcwd(),
        outdir = lambda wildcards: f"{wildcards.base}_naibr"
    shell:
        "mkdir -p {params.outdir}"
        " &&"
        " cd {input.naibr_path}"
        " &&"
        " python"
        " NAIBR.py"
        " {params.cwd}/{input.config}"
        " > {params.cwd}/{log}"
        " &&"
        " cd -"
        " &&"
        " mv {params.outdir}/NAIBR_SVs.bedpe {output.results}"
        " &&"
        " rm -rf {params.outdir}"


rule format_naibr_bedpe:
    output:
        bedpe = "final.naibr_sv_calls.bedpe"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    script:
        "../scripts/format_naibr_bedpe.py"


rule aggregate_sv_sizes:
    output:
        tsv = "final.sv_sizes.tsv"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    script:
        "../scripts/aggregate_sv_sizes.py"
