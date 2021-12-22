"""
Compile stats from barcode CLSTR file(s)

Files should be TSV with three columns: Canonical sequence, Read count, Comma-separated
Component sequences. E.g.

    TAAACATGCTCAGGAGCTAA	18	TAAACATGCTCAGGAACTAA,TAAACATGCTCAGGAGCTAA

Use indepentant from snakemake with command:

    python barcode_stats.py output.txt [input1.clstr [input2.clstr [...]]]
"""
import pandas as pd
import sys
from collections import OrderedDict

from blr import __version__
from blr.utils import smart_open


def main(clstrs, output):
    # Process data
    all_data = []
    for nr, clstr in enumerate(clstrs):
        data = pd.read_csv(clstr, sep="\t", names=["Canonical", "Reads", "Components"])
        data["FileNr"] = nr
        all_data.append(data)

    all_data = pd.concat(all_data)
    all_data["Size"] = all_data["Components"].apply(lambda x: len(x.split(',')))
    all_data["SeqLen"] = all_data["Canonical"].apply(len)

    # Collect stats
    stats = OrderedDict()
    stats["Barcodes raw"] = sum(all_data["Size"])
    stats["Barcodes corrected"] = len(all_data)
    stats["Barcodes corrected with > 3 read-pairs"] = len(all_data[all_data["Reads"] > 3])
    stats["Maximum reads per barcode"] = all_data["Reads"].max()
    stats["Mean reads per barcode"] = all_data["Reads"].mean()
    stats["Median reads per barcode"] = all_data["Reads"].median()

    reads_per_barcode = []
    total_reads = all_data["Reads"].sum()
    for reads, group in all_data.groupby("Reads"):
        nr_barcodes = len(group)
        total = reads * nr_barcodes
        density = total / total_reads
        reads_per_barcode.append((reads, nr_barcodes, total, density))

    components_per_barcode = OrderedDict(all_data.groupby("Size").apply(len))

    # Write output. Format is based on `samtools stats`
    with smart_open(output) as f:
        print(f"# Stats compiled from barcode_stats.py ({__version__})", file=f)
        print("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.", file=f)
        for stat, value in stats.items():
            print("SN", stat, value, sep="\t", file=f)

        print("# Reads per Barcode. Use `grep ^RB | cut -f 2-` to extract this part.", file=f)
        print("# Columns are: Reads, Nr of Barcodes, Total reads, Density.", file=f)
        for reads, count, total_reads, density in reads_per_barcode:
            print("RB", reads, count, total_reads, density, sep="\t", file=f)

        print("# Components per barcode. Use `grep ^CB | cut -f 2-` to extract this part.", file=f)
        print("# Columns are: Nr of components, Count.", file=f)
        for size, count in components_per_barcode.items():
            print("CB", size, count, sep="\t", file=f)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        clstrs = [snakemake.input.clstr]  # noqa: F821
        output = snakemake.output.txt  # noqa: F821
    elif len(sys.argv) >= 3:
        clstrs = sys.argv[2:]
        output = sys.argv[1]
    else:
        print(__doc__)
        sys.exit(1)

    main(clstrs, output)
