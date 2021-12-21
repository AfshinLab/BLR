"""
Break down NAIBR SVs into type and count their lengths in select bins.

Use indepentant from snakemake with command:

    python aggregate_sv_sizes.py input.tsv output.tsv

"""
from collections import Counter, defaultdict
import numpy as np
import sys

from blr.utils import parse_naibr_tsv


def main(input_tsv, output_tsv):
    counts = Counter()
    lengths = defaultdict(list)
    with open(input_tsv) as infile:
        for nr, sv in enumerate(parse_naibr_tsv(infile)):
            if sv.pass_filter != "PASS":
                continue

            if sv.chr1 != sv.chr2:
                counts[("Interchrom", sv.svtype())] += 1
            else:
                lengths[sv.svtype()].append(len(sv))

    # Aggregate counts over bins
    sv_types = ["DEL", "INV", "DUP"]
    bins = [0, 1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000]
    bin_labels = ["0-1k", "1k-10k", "10k-100k", "100k-1M", "1M-10M", "10M-100M", "100M-1G"]
    for sv_type in sv_types:
        bin_counts, _ = np.histogram(lengths[sv_type], bins=bins)

        for count, label in zip(bin_counts, bin_labels):
            counts[(label, sv_type)] += count

    with open(output_tsv, "w") as outfile:
        print("Size", *sv_types, sep="\t", file=outfile)
        for label in ["Interchrom"] + bin_labels:
            print(label, *[counts[(label, sv_type)] for sv_type in sv_types], sep="\t", file=outfile)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if len(sys.argv) == 3:
            input_tsv = sys.argv[1]
            output_tsv = sys.argv[2]
        else:
            print(__doc__)
            sys.exit(1)
    else:
        input_tsv = snakemake.input.tsv  # noqa: F821
        output_tsv = snakemake.output.tsv  # noqa: F821

    main(input_tsv, output_tsv)
