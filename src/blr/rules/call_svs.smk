import os
import pandas as pd
from snakemake.script import get_source


config["mtglink"] = {}
config["mtglink"]["min_contigsize"] = 10_000
config["mtglink"]["ins_extension"] = 50
config["mtglink"]["flanksize"] = 10_000
config["mtglink"]["minbarcocc"] = 2
config["mtglink"]["extsize"] = 500
config["mtglink"]["maxlength"] = 50_000
config["mtglink"]["kmers"] = "61 51 41 31 21"


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


rule manta:
    "Run manta for short SV calling"
    output:
        vcf = "final.manta_sv_calls.vcf.gz",
        cand_small = "final.manta_cand_small.vcf.gz",
        cand_small_idx = "final.manta_cand_small.vcf.gz.tbi",
        cand_sv = "final.manta_cand_sv.vcf.gz",
        cand_sv_idx = "final.manta_cand_sv.vcf.gz.tbi",
        rundir = directory("manta_tmp"),
        runscript = temp("manta_tmp/runWorkflow.py"),
    input:
        phased_cram = "final.phased.cram",
        phased_crai = "final.phased.cram.crai",
        bed = "primary.bed.gz",
        bed_tbi = "primary.bed.gz.tbi",
    conda: "../envs/manta.yml"
    log: "manta.log"
    threads: 20
    params:
        reference = config["genome_reference"],
        resultsdir = lambda wc, output: f"{output.rundir}/results/variants",
    shell:
        # 1. Setup manta
        "configManta.py"
        " --bam {input.phased_cram}"
        " --referenceFasta {params.reference}"
        " --callRegions {input.bed}"
        " --runDir {output.rundir} &> {log}"
        " && "
        # 2. Run manta
        "python2 {output.runscript} --jobs {threads} &>> {log}"
        " && "
        # 3. Run convertInversion.py and compress output
        "convertInversion.py"
        " $(which samtools)"
        " {params.reference}"
        " {params.resultsdir}/diploidSV.vcf.gz"
        " 2> {log}"
        " | "
        " bgzip -c"
        " > {params.resultsdir}/diploidSV.convertInversion.vcf.gz"
        " && "
        # 4. Copy outputs to workdir
        "cp {params.resultsdir}/diploidSV.convertInversion.vcf.gz {output.vcf}"
        " && "
        "cp {params.resultsdir}/candidateSmallIndels.vcf.gz {output.cand_small}"
        " && "
        "cp {params.resultsdir}/candidateSmallIndels.vcf.gz.tbi {output.cand_small_idx}"
        " && "
        "cp {params.resultsdir}/candidateSV.vcf.gz {output.cand_sv}"
        " && "
        "cp {params.resultsdir}/candidateSV.vcf.gz.tbi {output.cand_sv_idx}"


rule get_large_insertions:
    # https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#insertions-with-incomplete-insert-sequence-assembly
    input:
        vcf = "final.manta_sv_calls.vcf.gz",
    output:
        vcf = "final.manta_large_insertions.vcf"
    shell:
        "bcftools view"
        " -f 'PASS'"
        " -i 'SVTYPE=\"INS\" && (LEFT_SVINSSEQ!=\".\" || RIGHT_SVINSSEQ!=\".\")'"
        " -o {output.vcf}"
        " {input.vcf}"


rule create_insertion_gfa:
    # https://github.com/anne-gcd/MTG-Link/blob/851e4c14034a8a11be371c8e87561c222b4702df/utils/README.md#vcf2gfapy
    input:
        vcf = "final.manta_large_insertions.vcf",
    output:
        gfa = "mtglink.gfa",
        fastas = directory("mtglink_tmp/insertion_fastas/")
    log: "mtglink.gfa.log"
    conda: "../envs/mtglink.yml"
    params:
        vcf_abs = lambda wc, input: os.path.join(os.getcwd(), input.vcf),
        gfa_abs = lambda wc, output: os.path.join(os.getcwd(), output.gfa),
        reference = os.path.abspath(config["genome_reference"]),
        min_contigsize = config["mtglink"]["min_contigsize"], # Minimum size of the flanking contigs of the insertion to treat as a target
        extension = config["mtglink"]["ins_extension"], # Number of bases to extend each insertion
    shell: 
        # vcf2gfa.py outputs FASTAs with flankin contigs of the insertions in the
        # CWD so we need to change to the output directory
        "mkdir -p {output.fastas}"
        " && "
        "cd {output.fastas}"
        " && "
        "vcf2gfa.py"
        " -vcf {params.vcf_abs}"
        " -fa {params.reference}"
        " -out {output.gfa}"
        " -contigs {params.min_contigsize}"
        " -extension {params.extension}"
        " &> ../../{log}"
        " && "
        "mv {output.gfa} {params.gfa_abs}"
        " && "
        "cd -"


rule get_barcode_reads_bam:
    # Extract reads with the barcode tag from the phased CRAM file
    # TODO - this could maybe be skipped since we use the -bxuDir option in MTG-Link
    #        however extract_barcode_fastq may be faster with this prefiltering
    input:
        cram = "final.phased.cram",
        crai = "final.phased.cram.crai",
    output:
        bam = "final.barcodes.bam",
    params:
        barcode_tag = config["cluster_tag"],
    threads: 20
    shell:
        "samtools view"
        " -@ {threads}"
        " -bh"
        " -F 0x400" # No duplicates
        " -d {params.barcode_tag}"
        " -o {output.bam}"
        " {input.cram}"


checkpoint extract_barcodes:
    # Extract barcodes from the flanking regions of the insertions
    # TODO - this could maybe be done on the final.phased.cram directly
    input:
        gfa = ancient("mtglink.gfa"),
        bam = "final.barcodes.bam",
        bai = "final.barcodes.bam.bai",
    output:
        dir = directory("mtglink_tmp/read_subsampling_pre"),
    log: "mtglink_tmp/read_subsampling_pre.log",
    conda: "../envs/mtglink.yml"
    params:
        flanksize = config["mtglink"]["flanksize"],
        minbarcocc = config["mtglink"]["minbarcocc"],
    script:
        "../scripts/extract_barcodes.py"


rule extract_barcode_fastq:
    # Extract FASTQ with barcodes found in the list from the BAM file
    # TODO - this could maybe be done on the final.phased.cram directly
    input:
        lst = "mtglink_tmp/read_subsampling_pre/{file}.bxu",
        bam = "final.barcodes.bam",
        bai = "final.barcodes.bam.bai",
    output:
        bam = "mtglink_tmp/read_subsampling_pre/{file}.bxu.fastq",
    log: "mtglink_tmp/read_subsampling_pre/{file}.bxu.bam.log",
    threads: 2
    params:
        barcode_tag = config["cluster_tag"],
    shell:
        "samtools view"
        " {input.bam}"
        " -D {params.barcode_tag}:{input.lst}"
        " -@ {threads}"
        " -uh"
        " | "
        "samtools fastq"
        " -@ {threads}"
        " -T {params.barcode_tag}"
        " -"
        " > {output.bam}"
        " 2> {log}"


def aggregate_extracts_input(wildcards):
    checkpoint_out = checkpoints.extract_barcodes.get(**wildcards).output[0]
    lst = expand(os.path.join(checkpoint_out, "{file}.bxu.fastq"),
                  file=glob_wildcards(os.path.join(checkpoint_out, "{file}.bxu")).file)
    return lst


rule aggregate_extracts:
    # Required to make sure all the extract_barcode_fastq jobs are done
    input:
        aggregate_extracts_input,
    output:
        flag = "aggregate_extracts.done"
    run:
        for file in input:
            exists = os.path.exists(file)
            if not exists:
                raise Exception(f"File {file} does not exist")
        
        open(output.flag, "w").close()


rule touch_mtglink_input:
    # This FASTQ and index is not needed as we extract the reads from the BAM file
    # directly and used the -bxuDir option in MTG-Link.
    output:    
        fastq = touch("final.barcodes.fastq.gz"),
        bci = touch("final.barcodes.fastq.gz.bci"),


rule mtglink:
    # Run MTG-Link to assemble the insertions
    # https://github.com/anne-gcd/MTG-Link#options
    input:
        gfa = "mtglink.gfa", # Name must not contain any additional '.' except for the one before the extension
        bam = "final.barcodes.bam", # Does not work with CRAM, but if we use the -bxuDir option it should not be needed
        bai = "final.barcodes.bam.bai",
        fastq = "final.barcodes.fastq.gz",
        bci = "final.barcodes.fastq.gz.bci",
        bxudir = "mtglink_tmp/read_subsampling_pre",
        flag = "aggregate_extracts.done",
    output:
        fasta = "final.mtglink_insertions.fasta",
        bad_fasta = "final.mtglink_bad_insertions.fasta",
        read_subsampling = directory("mtglink_tmp/read_subsampling"),
    log: "mtglink.log",
    threads: 20
    conda: "../envs/mtglink.yml"
    params:
        fasta = "mtglink_tmp/mtglink.assembled_sequences.fasta",
        bad_fasta = "mtglink_tmp/bad_solutions.fasta",
        outdir = "mtglink_tmp",
        contigs = "mtglink_tmp/contigs",
        local_assembly = "mtglink_tmp/local_assembly",
        qual_evaluation = "mtglink_tmp/qual_evaluation",
        read_subsampling = "mtglink_tmp/read_subsampling",
        flanksize = config["mtglink"]["flanksize"],
        minbarcocc = config["mtglink"]["minbarcocc"],
        extsize = config["mtglink"]["extsize"],
        maxlength = config["mtglink"]["maxlength"],
        # minlength should be the length of the flanking regions + the length of 
        # the extended insert and have a minimum of insert size of 50
        minlength = 2*(config["mtglink"]["extsize"] + config["mtglink"]["ins_extension"]) + 50,
        k = config["mtglink"]["kmers"],
        max_mem = lambda wc, threads: threads * config["heap_space"]
    shell:
        # Remove previous runs
        "rm -rf"
        " {params.fasta}"
        " {params.bad_fasta}"
        " {params.outdir}/final.manta_large_insertions_mtglink.gfa"
        " {params.contigs}"
        " {params.local_assembly}"
        " {params.qual_evaluation}"
        " {params.read_subsampling}"
        " && " 
        "mtglink.py"
        " DBG"
        " -gfa {input.gfa}"
        " -bam {input.bam}"
        " -fastq {input.fastq}"
        " -index {input.bci}"
        " -out {params.outdir}"
        " -bxuDir {input.bxudir}"
        " -t 1"                         # Number of threads to use for the Read Subsampling step [default: 1]
        " -flank {params.flanksize}"    # Flanking sequences' size (bp) [default: 10000]
        " -occ {params.minbarcocc}"     # Minimum number of occurrences in target flanking regions for a barcode to be retained in the union set [default: 2]
        " -ext {params.extsize}"        # Size of the extension of the target on both sides (bp); determine start/end of local assembly [default: 500]
        " -l {params.maxlength}"        # Maximum assembly length (bp) [default: 10000]
        " -m {params.minlength}"        # Minimum assembly length (bp), by default 2*(-ext) bp [default: 1000]
        " -nb-cores {threads}"          # Number of cores to use for the Local Assembly step (DBG assembly) [default: 1]
        " -max-memory {params.max_mem}" # Maximum memory for graph building (in MBytes) [default: 0]
        " -verbose 0"                   # Verbosity level for DBG assembly [default: 0]
        " -k {params.k}"                # k-mer sizes for DBG assembly
        " &> {log}"
        " && "
        "cp {params.fasta} {output.fasta}"
        " && "
        "cp {params.bad_fasta} {output.bad_fasta}"


rule create_mtglink_vcf:
    input:
        fasta = "final.mtglink_insertions.fasta",
        vcf = "final.manta_large_insertions.vcf.gz",
        tbi = "final.manta_large_insertions.vcf.gz.tbi",
    output:
        vcf = "final.mtglink_insertions.vcf",
    log: "final.mtglink_insertions.vcf.log"
    params:
        flank = config["mtglink"]["ins_extension"] + config["mtglink"]["extsize"],
    script: 
        "../scripts/create_mtglink_vcf.py"
