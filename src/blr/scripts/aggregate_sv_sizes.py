"""
Break down NAIBR SVs into type and count their lengths in select bins.

Use indepentant from snakemake with command:

    python aggregate_sv_sizes.py input.tsv output.tsv

"""
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import sys


def main(input_tsv, output_tsv):
    sv_trans = {
        "+-": "DEL",
        "++": "INV",
        "--": "INV",
        "-+": "DUP"
    }

    sv_types = ["DEL", "INV", "DUP"]
    bins = [0, 1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000]
    bin_labels = ["0-1k", "1k-10k", "10k-100k", "100k-1M", "1M-10M", "10M-100M", "100M-1G"]

    names = [
        "Chr1",
        "Break1",
        "Chr2",
        "Break2",
        "SplitMolecules",
        "DiscordantReads",
        "Orientation",
        "Haplotype",
        "Score",
        "PassFilter"
    ]

    data = pd.read_csv(input_tsv, sep="\t", header=0, names=names)  # noqa: F821
    counts = Counter()
    lengths = defaultdict(list)
    with open(output_tsv, "w") as file:  # noqa: F821
        # Remove filtered SVs
        data = data[data["PassFilter"] == "PASS"]

        for nr, row in enumerate(data.itertuples(index=False)):
            sv_type = sv_trans[row.Orientation]
            if row.Chr1 != row.Chr2:
                counts[("Interchrom", sv_type)] += 1
            else:
                lengths[sv_type].append(abs(row.Break1 - row.Break2))

        for sv_type in sv_types:
            bin_counts, _ = np.histogram(lengths[sv_type], bins=bins)

            for count, label in zip(bin_counts, bin_labels):
                counts[(label, sv_type)] += count

        print("Size", *sv_types, sep="\t", file=file)
        for label in ["Interchrom"] + bin_labels:
            print(label, *[counts[(label, sv_type)] for sv_type in sv_types], sep="\t", file=file)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if len(sys.argv) == 3:
            input_tsv = sys.argv[1]
            output_tsv = sys.argv[2]
        else:
            print(__doc__)
            sys.exit(1)
    else:
        input_tsv = snakemake.input.tsv
        output_tsv = snakemake.output.tsv

    main(input_tsv, output_tsv)
