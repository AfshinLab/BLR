"""
Rules for trimming and demultiplexing of raw TELL-seq fastq files.
"""

localrules: tellseq_link_barcodes, tellseq_barcodes_correction, tag_tellseq_reads, merge_bins


rule tellseq_link_barcodes:
    output:
        "barcodes.fastq.gz"
    params:
        index = config["tellseq_index"]
    shell:
        "ln -s {params.index} {output}"


rule tellseq_barcodes_correction:
    """Correct barcodes"""
    output:
        temp("barcodes.clstr")
    input:
        "barcodes.fastq.gz"
    threads: 20 if config["tellseq_correction"] == "cluster" else 1
    log: "barcodes.clstr.log"
    run:
        command = {
            "cluster":
                "pigz -cd {input} |"
                " starcode"
                " -o {output}"
                " -t {threads}"
                " -d 1"
                " -r 2"
                " --print-clusters",
            "correct_singles":
                "blr correctbc"
                " {input}"
                " -o {output}"
        }[config["tellseq_correction"]]
        shell(f"{command} 2> {log}")


if config["read_mapper"] == "ema" and config["fastq_bins"] > 1:
    tag_output = temp(expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"]))
    output_cmd = f" --output-bins {config['_ema_bins_dir']} --nr-bins {config['fastq_bins']}"
    ruleorder: merge_bins > tag_tellseq_reads
else:
    tag_output = expand("trimmed.barcoded.{nr}.fastq.gz", nr=["1", "2"])
    output_cmd = f" --output1 {tag_output[0]} --output2 {tag_output[1]}"
    ruleorder: tag_tellseq_reads > merge_bins

if config["read_mapper"] == "ema":
    # Add non barcoded reads to output
    tag_output += expand("trimmed.non_barcoded.{nr}.fastq.gz", nr=["1", "2"])
    output_cmd += f" --output-nobc1 {tag_output[-2]} --output-nobc2 {tag_output[-1]}"


rule tag_tellseq_reads:
    """Tag reads with uncorrected and corrected barcode."""
    output:
        tag_output
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        uncorrected_barcodes="barcodes.fastq.gz",
        corrected_barcodes="barcodes.clstr.gz"
    log: "tagfastq.log"
    threads: 1
    params:
        output = output_cmd,
        barcode_tag = config["cluster_tag"],
        sequence_tag = config["sequence_tag"],
        mapper = config["read_mapper"],
        pattern = config["tellseq_barcode"],
        sample_nr = config["sample_nr"],
    shell:
        "blr tagfastq"
        " {params.output}"
        " -b {params.barcode_tag}"
        " -s {params.sequence_tag}"
        " --mapper {params.mapper}"
        " --pattern-match {params.pattern}"
        " --min-count 2"
        " --sample-nr {params.sample_nr}"
        " {input.uncorrected_barcodes}"
        " {input.corrected_barcodes}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " 2> {log}"


rule merge_bins:
    """Merge bins of trimmed and barcoded reads together"""
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz"
    input:
        bins = expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"]),
    params:
        modify_header = "" if config["read_mapper"]  == "ema" else " | tr ' ' '_' "
    shell:
        "cat {input.bins}"
        "{params.modify_header}"
        " |"
        " paste - - - - - - - -"
        " |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq})"
        " |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"
