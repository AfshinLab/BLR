import os


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


rule get_naibr_path:
    """Symlink existing NAIBR repo if specified or clone the latest version"""
    output:
        naibr_path = temp(directory("NAIBR"))
    log: "NAIBR.log"
    run:
        if config["naibr_path"]:
            shell("ln -s {config[naibr_path]} {output.naibr_path}")
        else:
            # Using forked branch with fix for NAIBR repo.
            # TODO Use original repo as below when https://github.com/raphael-group/NAIBR/pull/20 and 
            # https://github.com/raphael-group/NAIBR/pull/21 is merged.
            # shell("git clone https://github.com/raphael-group/NAIBR.git {output.naibr_path}")
            shell("git clone --branch fix-coverage https://github.com/pontushojer/NAIBR.git {output.naibr_path} &> {log}")


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
    conda: "../envs/naibr.yml"
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
        " &> {params.cwd}/{log}"
        " &&"
        " cd - > /dev/null"
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


rule format_naibr_vcf:
    output:
        vcf = temp("final.naibr_sv_calls.vcf")
    input:
        tsv = "final.naibr_sv_calls.tsv", 
        fai = config["genome_reference"] + ".fai"
    script:
        "../scripts/format_naibr_vcf.py"


rule aggregate_sv_sizes:
    output:
        tsv = "final.sv_sizes.tsv"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    script:
        "../scripts/aggregate_sv_sizes.py"
