"""
Rules for trimming and demultiplexing of raw BLR FASTQ files.

READ1 LAYOUT

5'-CAGTTGATCATCAGCAGGTAATCTGG BDVHBDVHBDVHBDVHBDVH CATGACCTCTTGGAACTGTCAGATGTGTATAAGAGACAG NNNN...NNNN (CTGTCTCTTATACACATCT)-3'
   <------------h1----------> <-------DBS--------> <-----------------h2------------------> <---gDNA--> <---------h3-------->
    h1 may inlude frameshift                                                                           Presence depends on
    oligos varying from 0-4                                                                            insert length
    extra oligos in 5' end

READ 2 LAYOUT

5'-NNNN...NNNN (CTGTCTCTTATACACATCT)-3'
   <---gDNA--> <---------h3-------->
               Presence depends on
               insert length
"""

barcode_placeholder = "N"*len(config["barcode"])
trim_len = sum(map(len, [config["h1"], barcode_placeholder, config["h2"]]))
extract_len = len(config["h1"])


rule trim:
    """Trim away 5' and possible 3' handles on read1 and trim possible 3' handles on read2."""
    output:
        interleaved_fastq=pipe("trimmed.fastq")
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
    log: "cutadapt_trim.log"
    threads: workflow.cores - 1  # rule tag needs one thread
    shell:
        "cutadapt"
        " -g 'XNNN{config[h1]}{barcode_placeholder}{config[h2]};min_overlap={trim_len}...{config[h3]};optional'"
        " -A {config[h3]}"
        " --pair-filter 'any'"
        " -e 0.2"
        " -j {threads}"
        " -m 25"
        " --interleaved"
        " -o {output.interleaved_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"


if config["read_mapper"] == "ema" and config["fastq_bins"] > 1:
    output_name = os.path.join(config['ema_bins_dir'], "ema-bin-{nr}")
    output_nrs = [str(i).zfill(3) for i in range(config['fastq_bins'])]
    output_cmd = f" --output-bins {config['ema_bins_dir']} --nr-bins {config['fastq_bins']}"
    ruleorder: merge_bins > tag
else:
    output_name = "trimmed.barcoded.{nr}.fastq.gz"
    output_nrs = ["1", "2"]
    output_cmd = f" --o1 {output_name.format(nr=output_nrs[0])} --o2 {output_name.format(nr=output_nrs[1])}"
    ruleorder: tag > merge_bins


rule tag:
    """Tag reads with uncorrected and corrected barcode."""
    output:
        expand(output_name, nr=output_nrs)
    input:
        interleaved_fastq="trimmed.fastq",
        uncorrected_barcodes="barcodes.fasta.gz",
        corrected_barcodes="barcodes.clstr"
    log: "tag_fastq.log"
    threads: 1
    shell:
        "blr tagfastq"
        " {output_cmd}"
        " -b {config[cluster_tag]}"
        " -s {config[sequence_tag]}"
        " --mapper {config[read_mapper]}"
        " --pattern-match {config[barcode]}"
        " --sample-nr {config[sample_nr]}"
        " {input.uncorrected_barcodes}"
        " {input.corrected_barcodes}"
        " {input.interleaved_fastq}"
        " 2> {log}"


rule extract_DBS:
    """Extract barcode sequence from read1 FASTQ"""
    output:
        fastq="barcodes.fasta.gz"
    input:
        fastq="reads.1.fastq.gz"
    log: "cutadapt_extract_DBS.log"
    threads: 20
    shell:
        "cutadapt"
        " -g 'XNNN{config[h1]};min_overlap={extract_len}...{config[h2]}'"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m 19"
        " -M 21"
	    " --max-n 0"
        " -o {output.fastq}"
        " {input.fastq}"
        " > {log}"


rule starcode_clustering:
    """Cluster DBS barcodes using starcode"""
    output:
        "barcodes.clstr"
    input:
        "barcodes.fasta.gz"
    threads: 20
    log: "starcode_clustering.log"
    shell:
        "pigz -cd {input} |"
        " starcode"
        " -o {output}"
        " -t {threads}"
        " -d {config[barcode_max_dist]}"
        " -r {config[barcode_ratio]}"
        " --print-clusters"
        " 2> {log}"


rule merge_bins:
    """Merge bins of trimmed and barcoded reads together"""
    output:
        r1_fastq="trimmed.barcoded.1.fastq.gz",
        r2_fastq="trimmed.barcoded.2.fastq.gz"
    input:
        bins = expand(os.path.join(config['ema_bins_dir'], "ema-bin-{nr}"),
                      nr=[str(i).zfill(3) for i in range(config['fastq_bins'])]),
    run:
        modify_header = "" if config["read_mapper"]  == "ema" else " | tr ' ' '_' "
        shell(
            "cat {input.bins}" +
            modify_header +
            " |"
            " paste - - - - - - - -"
            " |"
            " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq})"
            " |"
            " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"
            " &&"
            " rm {input.bins}"
        )
