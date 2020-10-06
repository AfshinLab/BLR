"""
Plot data from
"""
import logging
import pandas as pd
import numpy as np
from collections import OrderedDict, defaultdict, Counter
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from pathlib import Path
import os

from blr.utils import print_stats, calculate_N50
from blr import __version__

logger = logging.getLogger(__name__)

SIZE_WIDE = (10, 6)


def main(args):
    summary = Counter()

    # Fix for running matplotlib in background on SSH.
    # https://stackoverflow.com/questions/2443702/problem-running-python-matplotlib-in-background-after-ending-ssh-session
    matplotlib.use("Agg")

    name_to_function = {
        "molecule_stats": process_molecule_stats,   # Files containing "molecule_stats" and ending with "tsv"
        "barcode_clstrs": process_barcode_clstr     # Files named "barcodes.clstr"
    }

    # Make output directory if not present.
    args.output_dir.mkdir(exist_ok=True)

    matched = defaultdict(list)
    # Loop over file add assign to process.
    for filepath in args.input:
        filename = filepath.name
        if filename == "barcodes.clstr":
            matched["barcode_clstrs"].append(filepath)
            continue

        if "molecule_stats" in filename and filename.endswith(".tsv"):
            matched["molecule_stats"].append(filepath)
            continue

        logger.info(f"File '{filename}' does not match possible inputs. Skipping from analysis.")

    print(f"# Stats compiled from {__name__} ({__version__})")
    for func_name, files in matched.items():
        proc_func = name_to_function.get(func_name)
        proc_func(files, args.output_dir, summary)

    if summary:
        print_stats(summary, __name__)


def process_barcode_clstr(files, directory: Path, summary):
    if len(files) == 1:
        data = pd.read_csv(files[0], sep="\t", names=["Canonical", "Reads", "Components"])
        data["Size"] = data["Components"].apply(lambda x: len(x.split(',')))
        data["SeqLen"] = data["Canonical"].apply(len)
        plot_barcode_clstr(data, directory)
    else:
        logger.warning("Cannot handle multiple 'barcode.clstr' files.")


def plot_barcode_clstr(data: pd.DataFrame, directory: Path):
    # Histogram over reads per barcode cluster
    # - x = reads per cluster
    # - y = frequency
    with Plot("Reads per barcode cluster histogram", output_dir=directory) as (fig, ax):
        bins = np.logspace(0, np.log10(data['Reads'].max()), num=10)
        data["Reads"].plot(ax=ax, bins=bins, logy=True, logx=True, kind="hist")
        ax.set_xlabel("Reads per cluster")

    # Histogram over components per barcode cluster
    # - x = number of components per cluster
    # - y = frequency
    with Plot("Nr of components per barcode cluster histogram", output_dir=directory, figsize=SIZE_WIDE) as (fig, ax):
        data["Size"].plot(ax=ax, logy=True, bins=data["Size"].max(), kind="hist")
        ax.set_xlabel("Nr components per cluster")


def process_multiple(files, func):
    d = list()
    for nr, file in enumerate(files):
        d.append(func(file, nr=nr))
    return pd.concat(d)


def process_molecule_stats(files, directory: Path, summary):
    if len(files) == 1:
        data = process_molecule_stats_file(files[0])
    else:
        data = process_multiple(files, process_molecule_stats_file)

    summary["Barcodes"] = len(data["Barcode"].unique())
    summary["N50 reads per molecule"] = calculate_N50(data["Reads"])
    summary["Mean molecule length"] = float(data["Length"].mean())
    summary["Median molecule length"] = float(data["Length"].median())
    summary["DNA in molecules >20 kbp (%)"] = 100 * len(data[data["Length"] > 20_000]) / len(data)
    summary["DNA in molecules >100 kbp (%)"] = 100 * len(data[data["Length"] > 100_000]) / len(data)

    molecule_count = data.groupby("Barcode").count()["MoleculeID"]
    summary["Mean molecule count"] = float(molecule_count.mean())
    summary["Median molecule count"] = float(molecule_count.median())

    plot_molecule_stats(data, directory)


def process_molecule_stats_file(file, nr=None):
    data = pd.read_csv(file, sep="\t")
    if nr is not None and "ChunkID" not in data.columns.values:
        data["ChunkID"] = nr
    return data


def plot_molecule_stats(data: pd.DataFrame, directory: Path):
    # Histogram of molecule length distribution
    # - x = molecule length in kbp
    # - y = sum of lengths in bin
    with Plot("Molecule length histogram", output_dir=directory, figsize=SIZE_WIDE) as (fig, ax):
        bins, weights = bin_sum(data["Length"], binsize=2000)

        plt.hist(bins[:-1], bins, weights=list(weights))
        ax.set_ylabel("Total DNA mass")
        ax.set_yticklabels([])
        ax.set_xlabel("Molecule length (kbp)")
        ax.set_xticklabels(map(int, plt.xticks()[0] / 1000))

    # Read count per molecule vs molecule length
    # - x = molecule length in kbp
    # - y = read count for molecule
    with Plot("Molecule read count vs length", output_dir=directory) as (fig, ax):
        data.plot.hexbin(ax=ax, x="Length", y="Reads",  gridsize=20,
                         norm=LogNorm())
        ax.set_xlabel("Molecule length (kbp)")
        ax.set_xticklabels(map(int, plt.xticks()[0] / 1000))
        ax.set_ylabel("Reads per molecule")

    # Histogram of molecule read coverage
    # - x = molecule read coverage rate
    # - y = frequency
    coverage = data["BpCovered"] / data["Length"]
    with Plot("Molecule read coverage histogram", output_dir=directory) as (fig, ax):
        coverage.plot(ax=ax, kind="hist", xlim=(0, 1))
        ax.set_xlabel("Molecule read coverage")

    # Molecule coverage vs length
    # - x = molecule length in kbp
    # - y = molecule coverage in kbp
    with Plot("Molecule coverage vs length", output_dir=directory) as (fig, ax):
        data.plot.hexbin(ax=ax, x="Length", y="BpCovered", gridsize=20,
                         norm=LogNorm())
        ax.set_xlabel("Molecule length (kbp)")
        ax.set_xticklabels(map(int, plt.xticks()[0] / 1000))
        ax.set_ylabel("Molecule coverage (kbp)")
        ax.set_yticklabels(map(int, plt.yticks()[0] / 1000))

    # Molecules per barcode
    mols_per_bc = list(Counter(data.groupby("Barcode")["Barcode"].count().values).items())
    print("# Molecules per barcode. Use `grep ^MB | cut -f 2-` to extrat this part. Columns are: Molecules per "
          "barcode, Count.")
    for count, freq in sorted(mols_per_bc):
        print("MB", count, freq, sep="\t")


def bin_sum(data, binsize=2000):
    bins = range(0, max(data) + binsize, binsize)
    weights = OrderedDict({b: 0 for b in bins[:-1]})
    for value in data:
        current_bin = int(value / binsize) * binsize
        weights[current_bin] += value
    return bins, weights.values()


class Plot:
    """
    Plotting class for automatic filename generation, logging and file output. Using the defaults a PNG file with the
    suffix '_mqc.png' is outputed where 'mqc' makes the file detectable as custom content by MultiQC.
    """
    def __init__(self, title: str, output_dir: Path, figsize=(6.4, 4.8)):
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.title = title
        self.filename = self._make_filename()
        if output_dir:
            self.filepath = output_dir / self.filename
        else:
            self.filepath = self.filename

        logger.info(f"Ploting figure: {self.filename}")

    def _make_filename(self, suffix="_mqc.png"):
        """
        Create similar plot filenames with MultiQC suffix e.g "My title" -> 'My_Title_mqc.png'
        """
        capitalized = map(str.capitalize, self.title.lower().split(" "))
        return Path("_".join(capitalized) + suffix)

    def __enter__(self):
        return self.fig, self.ax

    def __exit__(self, exc_type, exc_val, exc_tb):
        plt.title(self.title)
        plt.tight_layout()
        plt.savefig(self.filepath)


def add_arguments(parser):
    parser.add_argument(
        "input", nargs="+", type=Path,
        help="Path to data files from pipeline run accepted file are. Currently accepted files are: barcodes.clstr, "
             "*molecule_stats.tsv."
    )
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=Path(os.getcwd()),
        help="Output directory for plot images. Default: CWD e.g %(default)s"
    )
