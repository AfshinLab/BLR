import os
import pandas as pd


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
        blacklist = f"--blacklist {config['naibr_blacklist']}" if config['naibr_blacklist'] else "",
        min_mapq = config["naibr_min_mapq"],
    shell:
        "blr naibrconfig"
        " --bam-file {input.bam}"
        " --outdir {params.cwd}/chunks/{wildcards.chunk}_naibr"
        " --distance 10000"
        " --min-mapq {params.min_mapq}"
        " --min-sv 1000"
        " --threads 1"
        " --min-overlaps 3"
        " {params.blacklist}"
        " --output {output.config}"
        " 2> {log}"


rule lsv_calling:
    """
    Runs NAIBR for LSV calling. This involves activating a python2 env, changing wd, running and changing back wd and
    env.
    """
    output:
        tsv = "chunks/{chunk}.naibr_sv_calls.tsv",
        bedpe = "chunks/{chunk}.naibr_sv_calls.bedpe",
        vcf = "chunks/{chunk}.naibr_sv_calls.vcf"
    input:
        config = "chunks/{chunk}.naibr.config",
    log: "chunks/{chunk}.naibr_sv_calls.tsv.log"
    threads: 2
    conda: "../envs/naibr.yml"
    params:
        cwd = os.getcwd(),
        outdir = lambda wildcards: f"chunks/{wildcards.chunk}_naibr"
    shell:
        "naibr {input.config}"
        " &> {log}"
        " &&"
        " mv {params.outdir}/NAIBR_SVs.bedpe {output.tsv}"
        " &&"
        " mv {params.outdir}/NAIBR_SVs.reformat.bedpe {output.bedpe}"
        " &&"
        " mv {params.outdir}/NAIBR_SVs.vcf {output.vcf}"
        " &&"
        " rm -rf {params.outdir}"


def concat_lsv_calls_bedpe_input(wildcards):
    return expand(
        "chunks/{chunk[0].name}.naibr_sv_calls.{filetype}", 
        chunk=chunks["phased"], 
        filetype=wildcards.filetype
    )

rule concat_lsv_calls:
    output:
        file = "final.naibr_sv_calls.{filetype,(bedpe|tsv)}"
    input:
        concat_lsv_calls_bedpe_input
    run:
        dfs = list()
        for file in input:
            try:
                dfs.append(pd.read_csv(file, sep="\t"))
            except pd.errors.EmptyDataError:
                continue

        concat = pd.concat(dfs, ignore_index=True)
        concat.to_csv(output.file, sep="\t", index=False)


def concat_lsv_calls_vcf_input(wildcards):
    return expand(
        "chunks/{chunk[0].name}.naibr_sv_calls.vcf", 
        chunk=chunks["phased"], 
    )


rule concat_lsv_calls_vcf:
    output:
        vcf = "final.naibr_sv_calls.vcf"
    input:
        concat_lsv_calls_vcf_input
    shell:
        "bcftools concat"
        " {input}"
        " | "
        "bcftools sort"
        " -o {output.vcf}"
        " -"


rule aggregate_sv_sizes:
    output:
        tsv = "final.sv_sizes.tsv"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    script:
        "../scripts/aggregate_sv_sizes.py"
