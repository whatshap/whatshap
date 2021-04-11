import pytest

from whatshap.core import Read, Variant
from whatshap.variants import merge_two_reads, merge_reads, ReadSetReader


@pytest.mark.parametrize("merge", [merge_two_reads, merge_reads])
def test_merge_pair_without_shared_positions(merge):
    empty1 = Read("Name1")
    empty2 = Read("Name2")
    assert merge(empty1, empty2).name == "Name1"
    assert merge(empty2, empty1).name == "Name2"

    # add_variant parameters are: (position, allele, quality)
    left = Read("Name1")
    left.add_variant(100, 0, 31)
    left.add_variant(200, 0, 32)
    right = Read("Name2")
    right.add_variant(300, 1, 41)
    right.add_variant(400, 1, 42)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 0, 32),
        Variant(300, 1, 41),
        Variant(400, 1, 42),
    ]
    assert expected == list(merge(left, right))
    assert expected == list(merge(right, left))

    outer = Read("Name1")
    outer.add_variant(100, 0, 31)
    outer.add_variant(400, 1, 42)
    inner = Read("Name2")
    inner.add_variant(200, 0, 32)
    inner.add_variant(300, 1, 41)
    assert expected == list(merge(inner, outer))
    assert expected == list(merge(outer, inner))


@pytest.mark.parametrize("merge", [merge_two_reads, merge_reads])
def test_merge_pair_with_shared_positions(merge):
    left = Read("Name1")
    left.add_variant(100, 0, 31)
    left.add_variant(200, 0, 32)
    left.add_variant(300, 0, 33)
    right = Read("Name2")
    right.add_variant(200, 0, 41)  # alleles disagree
    right.add_variant(300, 1, 42)  # alleles agree
    right.add_variant(400, 1, 43)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 0, 32 + 41),
        Variant(300, 1, 42),
        Variant(400, 1, 43),
    ]
    assert expected == list(merge(left, right))
    assert expected == list(merge(right, left))


def test_merge_many_reads():
    reads = [
        Read("Name1"),
        Read("Name2"),
        Read("Name3"),
    ]
    reads[0].add_variant(100, 0, 31)
    reads[0].add_variant(200, 1, 32)
    reads[0].add_variant(300, 0, 33)

    reads[1].add_variant(200, 1, 41)
    reads[1].add_variant(400, 0, 42)
    reads[1].add_variant(500, 0, 43)

    reads[2].add_variant(200, 0, 51)
    reads[2].add_variant(500, 0, 52)
    reads[2].add_variant(600, 0, 53)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 1, 73),  # see note: this depends on order of reads
        Variant(300, 0, 33),
        Variant(400, 0, 42),
        Variant(500, 0, 43 + 52),
        Variant(600, 0, 53),
    ]
    assert expected == list(merge_reads(*reads))

    # TODO merging should not depend on the order of reads
    expected[1] = Variant(200, 0, 51)
    assert expected == list(merge_reads(*reads[::-1]))


def test_split_by_distance():
    reads = [
        Read("Name1"),
        Read("Name2"),
        Read("Name3"),
        Read("Name4"),
    ]
    reads[0].add_variant(100, 0, 31)
    reads[0].add_variant(200, 0, 32)

    reads[1].add_variant(1000, 0, 41)
    reads[1].add_variant(1100, 0, 42)
    reads[1].add_variant(2000, 0, 42)

    # No variants on third read

    reads[3].add_variant(3000, 0, 51)
    reads[3].add_variant(3100, 0, 52)

    for distance, expected in [
        (10, [["Name1"], ["Name2", "Name3"], ["Name4"]]),
        (100, [["Name1"], ["Name2", "Name3"], ["Name4"]]),
        (799, [["Name1"], ["Name2", "Name3"], ["Name4"]]),
        (800, [["Name1", "Name2", "Name3"], ["Name4"]]),
        (999, [["Name1", "Name2", "Name3"], ["Name4"]]),
        (1000, [["Name1", "Name2", "Name3", "Name4"]]),
        (2000, [["Name1", "Name2", "Name3", "Name4"]]),
    ]:
        result = [
            [r.name for r in reads]
            for reads in ReadSetReader._split_by_distance(reads, max_distance=distance)
        ]
        assert expected == result, distance
