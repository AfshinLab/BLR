import pytest

from blr.utils import ACCEPTED_LIBRARY_TYPES
from blr.cli.find_clusterdups import get_non_acceptable_overlap_func, UnionFind


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


def test_union_find_generation():
    uf = UnionFind()
    uf.union("A", "B")
    uf.union("B", "C")

    assert frozenset(list(uf.connected_components())[0]) == frozenset({"A", "B", "C"})
    assert uf.same_component("A", "C")


def test_union_find_update():
    uf1 = UnionFind()
    uf1.union("A", "B")

    uf2 = UnionFind()
    uf2.union("B", "C")

    uf1.update(uf2)
    assert uf1.same_component("A", "C")


def test_union_find_from_dict():
    uf = UnionFind.from_dict({"A": "A", "B": "A", "C": "A"})
    assert uf.same_component("A", "C")
