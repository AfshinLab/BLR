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


if config["read_mapper"] == "ema" and config["fastq_bins"] > 1:
    tag_output = temp(expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"]))
    output_cmd = f" --output-bins {config['_ema_bins_dir']} --nr-bins {config['fastq_bins']}"
    ruleorder: merge_bins > tag_stlfr
else:
    tag_output = expand("trimmed.barcoded.{nr}.fastq.gz", nr=["1", "2"])
    output_cmd = f" --output1 {tag_output[0]} --output2 {tag_output[1]}"
    ruleorder: tag_stlfr > merge_bins

if config["read_mapper"] == "ema":
    # Add non barcoded reads to output
    tag_output += expand("trimmed.non_barcoded.{nr}.fastq.gz", nr=["1", "2"])
    output_cmd += f" --output-nobc1 {tag_output[-2]} --output-nobc2 {tag_output[-1]}"


rule tag_stlfr:
    """Modify header for downstream analysis."""
    output:
        tag_output
    input:
        interleaved_fastq = "trimmed.fastq",
    params:
        output = output_cmd,
        barcode_tag = config["cluster_tag"],
        mapper = config["read_mapper"],
        sample_nr = config["sample_nr"],
    log:
        log = "process_stlfr.log",
        csv = "process_stlfr.barcode_translations.csv"
    shell:
        "blr process_stlfr"
        " {params.output}"
        " -b {params.barcode_tag}"
        " --mapper {params.mapper}"
        " --sample-nr {params.sample_nr}"
        " --output-translations {log.csv}"
        " {input.interleaved_fastq}"
        " 2> {log.log}"


rule merge_bins:
    """Merge bins of trimmed and barcoded reads together"""
    output:
        interleaved_fastq=temp("trimmed.barcoded.fastq"),
    input:
        bins = expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"]),
    shell:
        "cat {input.bins} |"
        """ awk -F " " 'BEGIN{{OFS="\\n"}} {{print $2":"$1" BX:Z:"$1"-1",$3,"+",$4,$2":"$1" BX:Z:"$1"-1",$5,"+",$6}}'"""
        " > {output.interleaved_fastq}"


rule split_pairs:
    """Split into read pairs."""
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz"
    input:
        interleaved_fastq="trimmed.barcoded.fastq",
    shell:
        "paste - - - - - - - - < {input.interleaved_fastq} |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq}) |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"
