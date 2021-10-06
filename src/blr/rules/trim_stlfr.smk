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

localrules: trim_stlfr, tag_stlfr

rule trim_stlfr:
    """Trim away possible 3' adapter."""
    output:
        interleaved_fastq = pipe("trimmed.fastq")
    input:
        r1_fastq = "reads.1.fastq.gz",
        r2_fastq = "reads.2.fastq.gz",
    log: "trimmed.fastq.log",
    threads: workflow.cores - 1
    params:
        adapter = config["stlfr_adapter"],
    shell:
        "cutadapt "
        " -a {params.adapter}"
        " -A {params.adapter}"
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
        r2_fastq = "trimmed.barcoded.2.fastq.gz",
        tranlations = "process_stlfr.barcode_translations.csv",
    input:
        interleaved_fastq = "trimmed.fastq",
    params:
        barcodes = "" if config["stlfr_barcodes"] is None else f"--barcodes {config['stlfr_barcodes']}",
        barcode_tag = config["cluster_tag"],
        mapper = config["read_mapper"],
        sample_nr = config["sample_nr"],
    log: "process_stlfr.log"
    shell:
        "blr process_stlfr"
        " --o1 {output.r1_fastq}"
        " --o2 {output.r2_fastq}"
        " -b {params.barcode_tag}"
        " --mapper {params.mapper}"
        " --sample-nr {params.sample_nr}"
        " --output-translations {output.tranlations}"
        " {params.barcodes}"
        " {input.interleaved_fastq}"
        " 2> {log}"
