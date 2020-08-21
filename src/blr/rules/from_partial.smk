"""
Prepare input for running pipeline starting previously generated run.
"""
import sys
import glob
import pandas as pd

from blr.utils import  chromosome_chunks, parse_fai, symlink_relpath

# Generate chunks
if config["genome_reference"] is not None:
    try:
        with open(config["genome_reference"] + ".fai") as f:
            chunks = list(chromosome_chunks(parse_fai(f), size=config["chunk_size"]))
    except FileNotFoundError as e:
        sys.exit(
            f"The genome index file {config['genome_reference']} is missing. "
            "Please create it with 'samtools faidx'"
        )
else:
    chunks = []


rule all:
    input:
        "setup_from_input.done"


rule final:
    output:
        touch("setup_from_input.done")
    input:
        "trimmed.barcoded.1_fastqc.html",
        "trimmed.barcoded.2_fastqc.html",
        "barcodes.clstr",
        "final.molecule_stats.filtered.tsv",
        "unmapped.bam",
        expand("chunks/{chunk[0].name}.sorted.tag.bcmerge.mkdup.mol.filt.bam", chunk=chunks),
        "final.bam"


rule index_bam:
    output:
        bai = "{base}.bam.bai"
    input:
        bam = "{base}.bam"
    shell:
        "samtools index"
        " {input.bam}"
        " {output.bai}"


rule make_chunk_beds:
    output:
        expand("chunks/{chunk[0].name}.bed", chunk=chunks)
    run:
        for chunk in chunks:
            with open(f"chunks/{chunk[0].name}.bed", "w") as f:
                for chromosome in chunk:
                    print(chromosome.name, 0, chromosome.length, sep="\t", file=f)


rule split_input_into_chunks:
    output:
        bam = "chunks/{chunk}.sorted.tag.bcmerge.mkdup.mol.filt.bam"
    input:
        bam = "final.bam",
        bai = "final.bam.bai",
        bed = "chunks/{chunk}.bed",
    shell:
        "samtools view -M -L {input.bed} -o {output.bam} {input.bam}"


rule get_unmapped_reads_from_input:
    output:
        bam = "unmapped.bam"
    input:
        bam = "final.bam"
    shell:
        "samtools view -bhf 13 {input.bam} > {output.bam}"


rule touch_files:
    output:
        touch("trimmed.barcoded.1_fastqc.html"),
        touch("trimmed.barcoded.2_fastqc.html")


rule merge_or_link_input_bams:
    output:
        bam = "final.bam"
    input:
        dir = "inputs"
    run:
        input_bams = glob.glob(input.dir + "/*." + output.bam)
        print(input_bams)
        if len(input_bams) == 1:
            symlink_relpath(input_bams[0], output.bam)
        else:
            # TODO Run multi-threaded when https://github.com/snakemake/snakemake/issues/208 fixed or switched from
            #  subworkflow to 'include' format.
            shell("samtools merge -p {output.bam} {input_bams}")


rule concat_or_link_input_molecule_stats:
    output:
        tsv = "final.molecule_stats.filtered.tsv"
    input:
        dir = "inputs"
    run:
        input_files = glob.glob(input.dir + "/*." + output.tsv)
        if len(input_files) == 1:
            symlink_relpath(input_files[0], output.tsv)
        else:
            dfs = list()
            for file in input_files:
                try:
                    df = pd.read_csv(file, sep="\t")
                except pd.errors.EmptyDataError:
                    continue

                dfs.append(df)

            concat = pd.concat(dfs, ignore_index=True)
            concat.to_csv(output.tsv, sep="\t", index=False)


rule concat_barcode_clstrs:
    output:
        clstr = "barcodes.clstr"
    input:
        dir = "inputs"
    run:
        input_clstrs = glob.glob(input.dir + "/*." + output.clstr)
        if len(input_clstrs) == 1:
            symlink_relpath(input_clstrs[0], output.clstr)
        else:
            # TODO Tag barcodes based on input before concat
            shell("cat {input_clstrs} > {output.clstr}")
