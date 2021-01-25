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
data = pd.read_csv(snakemake.input.tsv, sep="\t", header=0, names=names)  # noqa: F821
with open(snakemake.output.bedpe, "w") as file:  # noqa: F821
    # Header
    print("#chrom1", "start1", "stop1", "chrom2", "start2", "stop2", "sv_type", "sv_id", "sv_length", "qual_score",
          "filter", "info", sep="\t", file=file)

    for nr, row in enumerate(data.itertuples(index=False)):
        # Use BEDPE format according to LinkedSV: https://github.com/WGLab/LinkedSV#-sv-call-file
        print(
            f"chr{row.Chr1}",
            row.Break1,
            row.Break1+1,
            f"chr{row.Chr2}",
            row.Break2,
            row.Break2+1,
            sv_types[row.Orientation],
            f"ID{nr:04}",
            abs(row.Break1 - row.Break2),
            row.Score,
            row.PassFilter,
            "ZYGOSITY={};NUM_SPLIT_MOLECULES={};NUM_DISCORDANT_READS={};ORIENTATION={};HAPLOTYPE={}".format(
                zygosity.get(row.Haplotype, "UNK"),
                row.SplitMolecules,
                row.DiscordantReads,
                row.Orientation,
                row.Haplotype,
            ),
            sep="\t",
            file=file
        )
