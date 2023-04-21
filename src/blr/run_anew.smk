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
]
if config["library_type"] in {"dbs", "blr", "tellseq"}:
    final_input.append("barcodes.clstr.gz")


rule final:
    input:
        final_input


input_bams = glob_wildcards("inputs/{name}.final.bam").name
input_crams = glob_wildcards("inputs/{name}.final.cram").name
input_tsvs = glob_wildcards("inputs/{name}.final.molecule_stats.filtered.tsv").name
input_clstrs = glob_wildcards("inputs/{name}.barcodes.clstr.gz").name


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
        bams = expand("inputs/{name}.final.bam", name=input_bams),
        bais = expand("inputs/{name}.final.bam.bai", name=input_bams),
        crams = expand("inputs/{name}.final.cram", name=input_crams),
        crais = expand("inputs/{name}.final.cram.crai", name=input_crams),
        bed = "chunks/{chunk}.bed",
    run:
        inputs = [*input.bams, *input.crams]
        if len(inputs) == 1:
            shell(
                "samtools view"
                " -x PC -x HP -x PS"  # Strip phasing tags 
                f" -o {output.bam}"
                f" -M -L {input.bed}"
                f" {inputs[0]}"
            )
        else:
            shell(
                "samtools merge"
                " -p"
                " -u"  # Output uncompressed BAM
                f" -L {input.bed}"
                f" - {inputs}"
                " |"
                " samtools view"
                " -x PC -x HP -x PS"  # Strip phasing tags 
                f" -o {output.bam}"
                " -"
            )

rule get_unmapped_reads_from_input:
    output:
        bam = "unmapped.bam",
        tmp_bams = temp(expand("inputs/{name}.unmapped.bam", name=input_bams+input_crams)),
    input:
        bams = expand("inputs/{name}.final.bam", name=input_bams),
        bais = expand("inputs/{name}.final.bam.bai", name=input_bams),
        crams = expand("inputs/{name}.final.cram", name=input_crams),
        crais = expand("inputs/{name}.final.cram.crai", name=input_crams),
    run:
        inputs = [*input.bams, *input.crams]
        if len(inputs) == 1:
            shell(f"samtools view -x PC -x HP -x PS {inputs[0]} '*' -o {output.bam}")
            shell("touch {output.tmp_bams}")
        else:
            for in_bam, tmp_bam in zip(inputs, output.tmp_bams):
                # Write temporary uncompressed BAM
                shell(f"samtools view -x PC -x HP -x PS -u {in_bam} '*' -o {tmp_bam}")

            # Concatenate BAMs
            shell(f"samtools cat -o {output.bam} {output.tmp_bams}")


rule touch_files:
    output:
        touch("trimmed.barcoded.1_fastqc.html"),
        touch("trimmed.barcoded.2_fastqc.html")


rule concat_or_link_input_molecule_stats:
    output:
        tsv = "final.molecule_stats.filtered.tsv"
    input:
        tsvs = expand("inputs/{name}.final.molecule_stats.filtered.tsv", name=input_tsvs),
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
        clstrs =  expand("inputs/{name}.barcodes.clstr.gz", name=input_clstrs)
    run:
        if len(input.clstrs) == 1:
            symlink_relpath(input.clstrs[0], output.clstr)
        else:
            # TODO Tag barcodes based on input before concat
            shell("cat {input.clstrs} > {output.clstr}")
