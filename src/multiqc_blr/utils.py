"""
Utility functions
"""

from collections import OrderedDict


def bin_sum(data, binsize=2000, normalize=False):
    """
    Sum values over bins and return a list of bin edges and bin weights.
    """
    bins = range(0, max(data) + binsize, binsize)
    weights = OrderedDict({b: 0 for b in bins[:-1]})
    for value in data:
        b = int(value / binsize) * binsize
        weights[b] += value

    weights = list(weights.values())
    if normalize:
        weights = [w / sum(weights) for w in weights]
    return bins, weights
