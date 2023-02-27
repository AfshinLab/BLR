"""
Prepare input for running pipeline starting from previously generated run(s).
"""
import pandas as pd
from snakemake.utils import validate

from blr.utils import  generate_chunks, symlink_relpath

configfile: "blr.yaml"
validate(config, "config.schema.yaml")

# Generate chunks
chunks = generate_chunks(reference=config["genome_reference"], 
                         size=config["chunk_size"], 
                         phasing_contigs_string=config["phasing_contigs"],
                         contigs_skipped=config["contigs_skipped"])

final_input = [
    "trimmed.barcoded.1_fastqc.html",
    "trimmed.barcoded.2_fastqc.html",
    "final.molecule_stats.filtered.tsv",
    "unmapped.bam",
    expand("chunks/{chunk[0].name}.calling.bam", chunk=chunks["all"]),
    "final.bam"
]
if config["library_type"] in {"dbs", "blr", "tellseq"}:
    final_input.append("barcodes.clstr.gz")


rule final:
    input:
        final_input


input_bams = glob_wildcards("inputs/{name}.bam").name
input_crams = glob_wildcards("inputs/{name}.cram").name
input_tsvs = glob_wildcards("inputs/{name}.tsv").name
input_clstrs = glob_wildcards("inputs/{name}.clstr.gz").name


rule index_bam:
    output:
        bai = "{base}.bam.bai"
    input:
        bam = "{base}.bam"
    threads: workflow.cores
    priority: 50
    shell:
        "samtools index"
        " -@ {threads}"
        " {input.bam}"
        " {output.bai}"


rule index_cram:
    output:
        crai = "{base}.cram.crai"
    input:
        cram = "{base}.cram"
    threads: workflow.cores
    shell:
        "samtools index -@ {threads} {input.cram} {output.crai}"


rule make_chunk_beds:
    output:
        expand("chunks/{chunk[0].name}.bed", chunk=chunks["all"])
    run:
        for chunk in chunks["all"]:
            with open(f"chunks/{chunk[0].name}.bed", "w") as f:
                for chromosome in chunk:
                    print(chromosome.name, 0, chromosome.length, sep="\t", file=f)


rule split_input_into_chunks:
    output:
        bam = "chunks/{chunk}.calling.bam"
    input:
        bam = "final.bam",
        bai = "final.bam.bai",
        bed = "chunks/{chunk}.bed",
    shell:
        "samtools view"
        # Strip phasing tags in case that input BAM was symlinked
        " -x PC -x HP -x PS"
        " -M -L {input.bed}"
        "-o {output.bam}"
        " {input.bam}"


rule get_unmapped_reads_from_input:
    output:
        bam = "unmapped.bam"
    input:
        bam = "final.bam",
        bai = "final.bam.bai"
    shell:
        "samtools view"
        # Strip phasing tags in case that input BAM was symlinked
        " -x PC -x HP -x PS"  
        " -bh"
        " -o {output.bam}"
        " {input.bam}"
        " '*'"


rule touch_files:
    output:
        touch("trimmed.barcoded.1_fastqc.html"),
        touch("trimmed.barcoded.2_fastqc.html")


rule merge_or_link_input_bams:
    output:
        bam = "final.bam"
    input:
        bams = expand("inputs/{name}.bam", name=input_bams),
        bais = expand("inputs/{name}.bam.bai", name=input_bams),
        crams = expand("inputs/{name}.cram", name=input_crams),
        crais = expand("inputs/{name}.cram.crai", name=input_crams),
    threads: 20
    priority: 50
    run:
        # If the input is a single BAM its ok to symlink here. CRAM or 
        # mulitple CRAM/BAM files need to be merged
        if len(input.bams) == 1:
            symlink_relpath(input.tsvs[0], output.tsv)
        else:
            shell(
                "samtools merge"
                " -p"
                " -@ {threads}"
                " - {input.bams} {input.crams}"
                " |"
                " samtools view"
                " -x PC -x HP -x PS"  # Strip phasing tags 
                " -o {output.bam}"
                " -"
            )


rule concat_or_link_input_molecule_stats:
    output:
        tsv = "final.molecule_stats.filtered.tsv"
    input:
        tsvs = expand("inputs/{name}.tsv", name=input_tsvs),
    run:
        if len(input.tsvs) == 1:
            symlink_relpath(input.tsvs[0], output.tsv)
        else:
            dfs = list()
            for file in input.tsvs:
                try:
                    df = pd.read_csv(file, sep="\t")
                except pd.errors.EmptyDataError:
                    continue

                dfs.append(df)

            concat = pd.concat(dfs, ignore_index=True)
            concat.to_csv(output.tsv, sep="\t", index=False)


rule concat_barcode_clstrs:
    output:
        clstr = "barcodes.clstr.gz"
    input:
        clstrs =  expand("inputs/{name}.clstr.gz", name=input_clstrs)
    run:
        if len(input.clstrs) == 1:
            symlink_relpath(input.clstrs[0], output.clstr)
        else:
            # TODO Tag barcodes based on input before concat
            shell("cat {input.clstrs} > {output.clstr}")
