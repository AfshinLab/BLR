"""
Compile stats from molecule TSV file(s)

Files should be TSV with at least columns: MoleculeID, Barcode, Reads, Length, BpCovered, ChunkID.

Use indepentant from snakemake with command:

    python molecules_stats.py output.txt [input1.tsv [input2.tsv [...]]]

"""
import pandas as pd
import sys
import numpy as np
from collections import OrderedDict, Counter

from blr import __version__
from blr.utils import smart_open, calculate_N50


def main(tsvs, output):
    # Process data
    all_data = []
    for nr, tsv in enumerate(tsvs):
        data = pd.read_csv(tsv, sep="\t")
        if not set(data.columns).issuperset({"MoleculeID", "Barcode", "Reads", "Length", "BpCovered", "ChunkID"}):
            print(f"# WARNING: File {tsv} does not have the expected columns.")
            continue

        data["FileNr"] = nr
        all_data.append(data)

    if len(all_data) == 0:
        print("ERROR: No data for stats compliation.")
        sys.exit(1)

    all_data = pd.concat(all_data)
    all_data["Coverage"] = all_data["BpCovered"] / all_data["Length"]

    # Collect stats
    stats = OrderedDict()
    stats["Barcodes final"] = len(all_data["Barcode"].unique())
    stats["N50 reads per molecule"] = calculate_N50(all_data["Reads"])
    stats["Mean reads per molecule"] = all_data["Reads"].mean()
    stats["Median reads per molecule"] = all_data["Reads"].median()
    stats["Mean molecule length"] = float(all_data["Length"].mean())
    stats["Median molecule length"] = float(all_data["Length"].median())

    total_dna = sum(all_data["Length"])
    stats["DNA in molecules >20 kbp (%)"] = 100 * sum(all_data[all_data["Length"] > 20_000]["Length"]) / total_dna
    stats["DNA in molecules >100 kbp (%)"] = 100 * sum(all_data[all_data["Length"] > 100_000]["Length"]) / total_dna
    stats["Weighted mean length"] = float(np.average(all_data["Length"].values, weights=all_data["Length"].values))

    molecule_count = all_data.groupby(["Barcode", "FileNr"]).count()["MoleculeID"]
    stats["Mean molecule count"] = float(molecule_count.mean())
    stats["Median molecule count"] = float(molecule_count.median())
    stats["Single molecule droplets (%)"] = float(100 * sum(molecule_count.values == 1) / len(molecule_count))

    dna_per_barcode = all_data.groupby(["Barcode", "FileNr"]).sum()["Length"]
    stats["Mean DNA per barcode"] = float(dna_per_barcode.mean())
    stats["Median DNA per barcode"] = float(dna_per_barcode.median())

    # Molecule read coverage
    coverage_bins = np.array(range(0, 101)) / 100
    all_data["Bin"] = pd.cut(all_data["Coverage"], bins=coverage_bins, labels=coverage_bins[:-1])
    binned_coverage = all_data.groupby("Bin", as_index=False)["Coverage"].count()
    binned_coverage = binned_coverage[binned_coverage["Coverage"] > 0]  # Remove zero entries

    # Molecules per barcode
    mols_per_bc = list(Counter(all_data.groupby(["Barcode", "FileNr"])["Barcode"].count().values).items())

    # Reads per barcode
    readcounts = all_data.groupby(["Barcode", "FileNr"], as_index=False)["Reads"].sum()
    read_bins = list(range(0, max(readcounts["Reads"])+2, 2))
    readcounts["Bin"] = pd.cut(readcounts["Reads"], bins=read_bins, labels=read_bins[:-1], right=False)
    binned_counts = readcounts.groupby("Bin", as_index=False)["Reads"].count()
    binned_counts = binned_counts[binned_counts["Reads"] > 0]  # Remove zero entries

    # Write output. Format is based on `samtools stats`
    with smart_open(output) as f:
        print(f"# Stats compiled from molecule_stats.py ({__version__})", file=f)
        print(f"# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.", file=f)
        for stat, value in stats.items():
            print("SN", stat, value, sep="\t", file=f)

        print("# Molecule coverage. Use `grep ^MC | cut -f 2-` to extract this part.", file=f)
        print("# Columns are: Molecule coverage bin (numbers refer to lower threshold), Count.", file=f)
        for row in binned_coverage.itertuples():
            print("MC", row.Bin, row.Coverage, sep="\t", file=f)

        print("# Molecules per barcode. Use `grep ^MB | cut -f 2-` to extract this part.", file=f)
        print("# Columns are: Molecules per barcode, Count.", file=f)
        for count, freq in sorted(mols_per_bc):
            print("MB", count, freq, sep="\t", file=f)

        print("# Reads per barcode. Use `grep ^RB | cut -f 2-` to extract this part", file=f)
        print("# Columns are: Reads bin (numbers relate to the lower threshold for the bin), Nr of barcodes", file=f)
        for row in binned_counts.itertuples():
            print("RB", row.Bin, row.Reads, sep="\t", file=f)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        tsvs = [snakemake.input.tsv]  # noqa: F821
        output = snakemake.output.txt  # noqa: F821
    elif len(sys.argv) >= 3:
        tsvs = sys.argv[2:]
        output = sys.argv[1]
    else:
        print(__doc__)
        sys.exit(1)

    main(tsvs, output)
