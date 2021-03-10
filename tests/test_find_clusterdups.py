import pytest

from blr.utils import ACCEPTED_LIBRARY_TYPES
from blr.cli.find_clusterdups import get_non_acceptable_overlap_func


@pytest.mark.parametrize("library_type", ACCEPTED_LIBRARY_TYPES)
def test_get_non_acceptable_overlap_func(library_type):
    test_values = [100, 0, -1, -5, -9, -100]
    acceptable = {
        "blr":     [0, 0, 1, 1, 0, 1],
        "stlfr":   [0, 0, 1, 1, 0, 1],
        "10x":     [0, 0, 0, 0, 0, 0],
        "tellseq": [0, 0, 1, 0, 1, 1],
    }
    non_acceptable_overlap = get_non_acceptable_overlap_func(library_type)
    for x, ref in zip(test_values, acceptable[library_type]):
        assert non_acceptable_overlap(x) == bool(ref)
