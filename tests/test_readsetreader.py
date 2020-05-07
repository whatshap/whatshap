from whatshap.core import Read, Variant
from whatshap.variants import merge_two_reads


def test_merge_pair_without_shared_positions():
    empty1 = Read("Name1")
    empty2 = Read("Name2")
    assert merge_two_reads(empty1, empty2).name == "Name1"
    assert merge_two_reads(empty2, empty1).name == "Name2"

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
    assert expected == list(merge_two_reads(left, right))
    assert expected == list(merge_two_reads(right, left))

    outer = Read("Name1")
    outer.add_variant(100, 0, 31)
    outer.add_variant(400, 1, 42)
    inner = Read("Name2")
    inner.add_variant(200, 0, 32)
    inner.add_variant(300, 1, 41)
    assert expected == list(merge_two_reads(inner, outer))
    assert expected == list(merge_two_reads(outer, inner))


def test_merge_pair_with_shared_positions():
    left = Read("Name1")
    left.add_variant(100, 0, 31)
    left.add_variant(200, 0, 32)
    left.add_variant(300, 0, 33)
    right = Read("Name2")
    right.add_variant(200, 0, 41)  # alleles disagree
    right.add_variant(300, 1, 41)  # alleles agree
    right.add_variant(400, 1, 42)

    expected = [
        Variant(100, 0, 31),
        Variant(200, 0, 32 + 41),
        Variant(300, 1, 41),
        Variant(400, 1, 42),
    ]
    assert expected == list(merge_two_reads(left, right))
    assert expected == list(merge_two_reads(right, left))
