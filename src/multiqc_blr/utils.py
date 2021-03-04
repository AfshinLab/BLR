"""
Utility functions
"""

from collections import OrderedDict
import pandas as pd


def bin_sum(data, binsize=2000, normalize=False):
    """
    Sum values over bins and return a list of bin edges and bin weights.
    """
    bins = range(0, max(data) + 2*binsize, binsize)
    weights = OrderedDict({b: 0 for b in bins[:-1]})
    for value in data:
        b = int(value / binsize) * binsize
        weights[b] += value

    weights = list(weights.values())
    if normalize:
        weights = [w / sum(weights) for w in weights]
    return bins, weights


def get_tail_x(data, threshold=0.999):
    """
    Get x value from dict of x: y values for lineplot for which to set the xmax to avoid long tails
    """
    d = pd.DataFrame(data.items(), index=range(len(data)))
    d.columns = ["x", "y"]
    d.sort_values(by="x", ascending=True, inplace=True)
    d["y_cumsum"] = d["y"].cumsum()
    d["p_of_tot"] = d["y_cumsum"] / d["y"].sum()
    return int(d.loc[d["p_of_tot"].gt(threshold).idxmax(), "x"])
