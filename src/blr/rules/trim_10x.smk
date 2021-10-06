"""
Rules to process 10xGenomics Chromium linked-read raw FASTQs.

READ1 LAYOUT

5'-NNNNNNNNNNNNNNNN NNNNNNN NNN...NNN-3'
   <----barcode---> <-UMI-> <-gDNA-->
         16 nt        7 nt

READ2 LAYOUT

5'-NNNN...NNNN-3'
   <---gDNA-->

Processing is partly based on the 10x end-to-end workflow described for the EMA aligner. See the docs on their github
https://github.com/arshajii/ema#end-to-end-workflow-10x
"""

localrules: link_to_whitelist, count_10x, preproc_10x, merge_bins, split_pairs, split_nobc_reads

# If not mapping with ema it is unneccessary to generate a lot of bins.
if config["read_mapper"] != "ema":
    config["fastq_bins"] = min(config["fastq_bins"], workflow.cores)
    config["_fastq_bin_nrs"] = [str(i).zfill(3) for i in range(config['fastq_bins'])]


rule link_to_whitelist:
    output: "barcodes_whitelist.txt"
    shell: "ln -s {config[barcode_whitelist]} {output}"


rule count_10x:
    """Create list of per-barcode count"""
    output:
        counts_ncnt = temp("reads.ema-ncnt"),
        counts_fcnt = temp("reads.ema-fcnt")
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        whitelist="barcodes_whitelist.txt"
    log: "ema_count.log"
    shell:
        "paste <(pigz -c -d {input.r1_fastq} | paste - - - -) <(pigz -c -d {input.r2_fastq} | paste - - - -) |"
        " tr '\t' '\n' |"
        " ema count"
        " -w {input.whitelist}"
        " -o reads 2> {log}"


rule preproc_10x:
    """Trim reads and bin reads containing the same barcode together. Reads missing barcodes outputed to ema-nobc."""
    output:
        bins = temp(expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"])),
        nobc = temp(config['_ema_bins_dir'] / "ema-nobc")
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
        counts_ncnt = "reads.ema-ncnt",
        counts_fcnt = "reads.ema-fcnt",
        whitelist = "barcodes_whitelist.txt",
    log: "ema_preproc.log"
    threads: 20
    params:
        hamming_correction = "" if not config["apply_hamming_correction"] else " -h",
        bins = config["fastq_bins"],
        bins_dir = config["_ema_bins_dir"]
    shell:
        "paste <(pigz -c -d {input.r1_fastq} | paste - - - -) <(pigz -c -d {input.r2_fastq} | paste - - - -) |"
        " tr '\t' '\n' |"
        " ema preproc"
        " -w {input.whitelist}"
        " -n {params.bins}"
        "{params.hamming_correction}"
        " -t {threads}"
        " -o {params.bins_dir} {input.counts_ncnt} 2>&1 | tee {log}"


rule merge_bins:
    """Merge bins of trimmed and barcoded reads together"""
    output:
        interleaved_fastq=temp("trimmed.barcoded.fastq"),
    input:
        bins = expand(config['_ema_bins_dir'] / "ema-bin-{nr}", nr=config["_fastq_bin_nrs"]),
        nobc = config['_ema_bins_dir'] / "ema-nobc"
    params:
        mapper = config["read_mapper"],
    run:
        # Format of rows is <barcode> <r1_name> <r1_seq> <r1_qual> <r2_seq> <r2_qual>
        if params.mapper == "ema":
            # Header is <name>:<barcode> BX:Z:<barcode>-1
            shell(
                "cat {input.bins} |"
              """ awk -F " " 'BEGIN{{OFS="\\n"}} {{print $2":"$1" BX:Z:"$1"-1",$3,"+",$4,$2":"$1" BX:Z:"$1"-1",$5,"+",$6}}'"""
                " > {output.interleaved_fastq}"
            )
        else:
            # Header is <name>_BX:Z:<barcode>-1
            shell(
                "cat {input.bins} |"
              """ awk -F " " 'BEGIN{{OFS="\\n"}} {{print $2"_BX:Z:"$1"-1",$3,"+",$4,$2"_BX:Z:"$1"-1",$5,"+",$6}}'"""
                " > {output.interleaved_fastq}"
                " &&"
                " cat {input.nobc} >> {output.interleaved_fastq}" # Include non-barcoded reads.
            )


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


rule split_nobc_reads:
    """Split non-barcoded reads into read pairs."""
    output:
        r1_fastq="trimmed.non_barcoded.1.fastq.gz",
        r2_fastq="trimmed.non_barcoded.2.fastq.gz",
    input:
        fastq = config["_ema_bins_dir"] / "ema-nobc"
    shell:
        "paste - - - - - - - - < {input} |"
        " tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > {output.r1_fastq}) |"
        " cut -f 5-8 | tr '\t' '\n' | pigz -c > {output.r2_fastq}"
