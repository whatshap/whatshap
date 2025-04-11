from whatshap.core import Read
from whatshap.variants import ReadSetReader as Reader, AlignedRead
import pytest


def test_supplementary_alignment():
    r = Read("S1", 60, 0, 0, 10, "")
    r.add_variant(10, 0, 60)
    ret = Reader.create_read_from_group([AlignedRead(r, True, False, 10, 20)], 10)
    assert ret is None


def test_primary_alignment():
    r = Read("P1", 60, 0, 0, 10)
    r.add_variant(10, 0, 60)
    ret = Reader.create_read_from_group([AlignedRead(r, False, False, 10, 20)], 10)
    assert len(ret) == 1


@pytest.mark.parametrize("rev1,rev2", [(False, False), (False, True), (True, False), (True, True)])
def test_two_primary_alignment(rev1, rev2):
    r1 = Read("P1", 60, 0, 0, 10)
    r1.add_variant(10, 0, 60)
    r2 = Read("P1", 60, 0, 0, 10)
    r2.add_variant(15, 1, 60)
    ret = Reader.create_read_from_group(
        [AlignedRead(r1, False, rev1, 10, 20), AlignedRead(r2, False, rev2, 10, 20)],
        distance_threshold=10,
    )
    assert len(ret) == 2


def test_three_primary_alignment():
    r1 = Read("P1", 60, 0, 0, 10)
    r1.add_variant(10, 0, 60)
    r2 = Read("P1", 60, 0, 0, 10)
    r2.add_variant(15, 1, 60)
    r3 = Read("P1", 60, 0, 0, 10)
    r3.add_variant(20, 1, 60)
    ret = Reader.create_read_from_group(
        [
            AlignedRead(r1, False, False, 10, 30),
            AlignedRead(r2, False, False, 10, 30),
            AlignedRead(r3, False, False, 10, 30),
        ],
        distance_threshold=10,
    )
    assert ret is None


def test_two_alignments_same_orientation():
    primary = Read("P1", 60, 0, 0, 10)
    primary.add_variant(10, 0, 60)
    supplementary = Read("S1", 60, 0, 0, 10)
    supplementary.add_variant(10, 0, 60)
    supplementary.add_variant(20, 0, 60)
    ret = Reader.create_read_from_group(
        [AlignedRead(primary, False, True, 10, 20), AlignedRead(supplementary, True, True, 10, 30)],
        100,
    )
    assert len(ret) == 2


def test_two_alignments_different_orientation():
    primary = Read("P1", 60, 0, 0, 10)
    primary.add_variant(10, 0, 60)
    supplementary = Read("S1", 60, 0, 0, 10)
    supplementary.add_variant(10, 0, 60)
    supplementary.add_variant(20, 0, 60)
    ret = Reader.create_read_from_group(
        [
            AlignedRead(primary, False, True, 10, 20),
            AlignedRead(supplementary, True, False, 10, 30),
        ],
        100,
    )
    assert len(ret) == 1


def test_distance():
    primary = Read("P1", 60, 0, 0, 10)
    primary.add_variant(10, 0, 60)
    supplementary = Read("S1", 60, 0, 0, 10)
    supplementary.add_variant(10, 0, 60)
    supplementary.add_variant(20, 0, 60)
    ret = Reader.create_read_from_group(
        [
            AlignedRead(primary, False, True, 10, 11),
            AlignedRead(supplementary, True, True, 20, 30),
        ],
        5,
    )
    assert len(ret) == 1
