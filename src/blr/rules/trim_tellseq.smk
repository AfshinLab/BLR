"""
Rules for trimming and demultiplexing of raw TELL-seq fastq files.
"""

rule tellseq_link_barcodes:
    output: "barcodes.fastq.gz"
    shell: "ln -s {config[tellseq_index]} {output}"


rule tellseq_barcodes_correction:
    """Correct barcodes"""
    output:
        "barcodes.clstr"
    input:
        "barcodes.fastq.gz"
    threads: 20 if config["tellseq_correction"] == "cluster" else 1
    log: "barcodes.clstr.log"
    run:
        commands = {
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
        }
        shell(commands[config["tellseq_correction"]] + " 2> {log}")


rule tag_tellseq_reads:
    """Tag reads with uncorrected and corrected barcode."""
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz"
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        uncorrected_barcodes="barcodes.fastq.gz",
        corrected_barcodes="barcodes.clstr"
    log: "tag_tellseq_reads.log"
    threads: 1
    shell:
        "blr tagfastq"
        " --o1 {output.r1_fastq}"
        " --o2 {output.r2_fastq}"
        " -b {config[cluster_tag]}"
        " -s {config[sequence_tag]}"
        " --mapper {config[read_mapper]}"
        " --pattern-match {config[tellseq_barcode]}"
        " --skip-singles"
        " {input.uncorrected_barcodes}"
        " {input.corrected_barcodes}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " 2> {log}"
