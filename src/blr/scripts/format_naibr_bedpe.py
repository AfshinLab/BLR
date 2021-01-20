"""
Fromat BEDPE from NAIBR output TSV
"""
import pandas as pd

# Translation based on this: https://github.com/raphael-group/NAIBR/issues/11
# Note that this is not entirely accurate as the more complex variants are possible.
sv_types = {"+-": "DEL", "++": "INV", "--": "INV", "-+": "DUP"}
zygosity = {"1,1": "HOM", "2,2": "HOM", "1,2": "HET",
            "2,1": "HET"}  # From https://github.com/raphael-group/NAIBR/issues/10
names = ["Chr1", "Break1", "Chr2", "Break2", "SplitMolecules", "DiscordantReads", "Orientation",
         "Haplotype", "Score", "PassFilter"]
data = pd.read_csv(snakemake.input.tsv, sep="\t", header=0, names=names)
with open(snakemake.output.bedpe, "w") as file:
    for nr, row in enumerate(data.itertuples(index=False)):
        # Use BEDPE format according to 10x: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bedpe
        # Supported by IGV (see https://github.com/igvteam/igv/wiki/BedPE-Support)
        print(
            f"chr{row.Chr1}",
            row.Break1,
            row.Break1,
            f"chr{row.Chr2}",
            row.Break2,
            row.Break2,
            f"call_{nr}",
            row.Score,
            "+",
            "+",
            "." if row.PassFilter == "PASS" else "Filtered",
            "Type={};Zygosity={};Split_molecules={};Discordant_reads={};Orientation={};Haplotype={}".format(
                sv_types[row.Orientation],
                zygosity.get(row.Haplotype, "UNK"),
                row.SplitMolecules,
                row.DiscordantReads,
                row.Orientation,
                row.Haplotype,
            ),
            sep="\t",
            file=file
        )
