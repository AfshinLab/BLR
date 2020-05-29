"""
Rules for trimming stLFR type reads


READ1 LAYOUT

5'-NNNN...NNNN (CTGTCTCTTATACACATCT)-3'
   <---gDNA--> <---------h3-------->
               Presence depends on
               insert length

READ 2 LAYOUT

5'-NNNN...NNNN (CTGTCTCTTATACACATCT)-3'
   <---gDNA--> <---------h3-------->
               Presence depends on
               insert length


"""

rule link_to_barcodes:
    output: "barcodes.txt"
    shell: "ln -s {config[stlfr_barcodes]} {output}"


rule trim_stlfr:
    """Trim away possible 3' adapter."""
    output:
        interleaved_fastq = pipe("trimmed.fastq")
    input:
        r1_fastq = "reads.1.fastq.gz",
        r2_fastq = "reads.2.fastq.gz",
        barcodes = "barcodes.txt"
    log: "cutadapt_trim.log",
    threads: workflow.cores - 1
    shell:
        "cutadapt "
        " -a {config[stlfr_adapter]}"
        " -A {config[stlfr_adapter]}"
        " -e 0.2"
        " -j {threads}"
        " --pair-filter 'both'"
        " -m 25"
        " --interleaved"
        " -o {output.interleaved_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"


rule tag_stlfr:
    """Modify header for downstream analysis."""
    output:
        r1_fastq = "trimmed.barcoded.1.fastq.gz",
        r2_fastq = "trimmed.barcoded.2.fastq.gz"
    input:
        interleaved_fastq = "trimmed.fastq",
        barcodes = "barcodes.txt"
    log: "process_stlfr.log"
    shell:
        "blr process_stlfr"
        " --o1 {output.r1_fastq}"
        " --o2 {output.r2_fastq}"
        " -b {config[cluster_tag]}"
        " --mapper {config[read_mapper]}"
        " {input.barcodes}"
        " {input.interleaved_fastq}"
        " 2> {log}"
