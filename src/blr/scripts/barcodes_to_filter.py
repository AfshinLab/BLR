"""
Generate a list of barcodes to filter out
"""
import pandas as pd
import sys
from scipy.special import gammaln
import math
from functools import lru_cache
import argparse

from blr.utils import parse_fai, smart_open


@lru_cache(maxsize=None)
def probability_collision(r: int, N: int):
    """Probablity for collision when choosing r samples from N possible choises."""
    # Based on: https://dzone.com/articles/calculating-birthday-paradox-scipy
    if N >= r:
        return 1 - math.exp(gammaln(N+1) - gammaln(N-r+1) - r*math.log(N))
    return 1


def main(
    input,
    output,
    threshold,
    reference,
    window
):
    # Process data
    barcodes_to_filter = []
    if threshold >= 1:
        # Count nr of molecules per barcode
        molecules = pd.read_csv(input, sep="\t")
        molecules_per_barcode = molecules.groupby("Barcode", sort=False)["Barcode"].count()
        barcodes_to_filter = molecules_per_barcode[molecules_per_barcode > threshold].index.to_list()

    elif threshold > 0:
        # Calculate the probability for molecule collisions for each barcode
        # Output barcodes above the definede probability threshold.

        with open(reference + ".fai") as f:
            genome_size = sum(chromosome.length for chromosome in parse_fai(f))

        # Divide genome into a number of bins based on the total size.
        total_bins = (genome_size // window) + 1

        molecules = pd.read_csv(input, sep="\t")

        # Caculate the number of windows required to cover each molecule
        # The +3 is so that each molecule is expanded both left and right by one window as
        # read within this span would be assigned together in `buildmolecules` and thus seem
        # to belong to the same molecule.
        molecules["NrBinsCovered"] = (molecules["Length"] // window) + 3

        # Sum the number of windows for each barcodes, this is number of samples we us for calculating the
        # probability for collisions.
        barcodes = molecules.groupby("Barcode", sort=False, as_index=False).agg({"NrBinsCovered": "sum"})

        # Calculate the probability for collision
        barcodes["p_value"] = barcodes["NrBinsCovered"].apply(lambda x: probability_collision(x, total_bins))
        barcodes_to_filter = barcodes[barcodes["p_value"] > threshold]["Barcode"].to_list()

    with smart_open(output) as file:
        print("\n".join(barcodes_to_filter), file=file)


if __name__ == "__main__":
    if len(sys.argv) == 1 and "snakemake" in globals():
        main(
            input=snakemake.input.tsv,
            output=snakemake.output.tsv,
            threshold=snakemake.params.threshold,
            reference=snakemake.params.reference,
            window=snakemake.params.window
        )
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("input")
        parser.add_argument("-o", "--output", default="-")
        parser.add_argument("-t", "--threshold", type=float)
        parser.add_argument("-r", "--reference")
        parser.add_argument("-w", "--window", type=int)
        args = parser.parse_args()
        main(
            input=args.input,
            output=args.output,
            threshold=args.threshold,
            reference=args.reference,
            window=args.window
        )
