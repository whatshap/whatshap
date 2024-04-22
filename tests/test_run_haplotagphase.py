"""
Tests for whatshap haplotagphase module
"""

from whatshap.cli.haplotagphase import compute_votes, run_haplotagphase
from whatshap.core import Read


def test_haplotagphase():
    run_haplotagphase(
        variant_file="tests/data/pacbio/variants.vcf",
        alignment_file="tests/data/pacbio/haplotagged.bam",
        reference="tests/data/pacbio/reference.fasta",
        output="/dev/null",
    )


def test_compute_votes():
    a = Read("a", 60, 0, 0, 0, "", 1, 1)
    a.add_variant(1, 0, 30)
    a.add_variant(2, 0, 10)
    a.add_variant(3, 0, 50)
    b = Read("b", 60, 0, 0, 0, "", 2, 1)
    b.add_variant(1, 1, 20)
    b.add_variant(2, 0, 30)
    b.add_variant(3, 0, 90)
    c = Read("c", 60, 0, 0, 0, "", 1, 2)
    c.add_variant(1, 1, 20)
    c.add_variant(3, 0, 10)
    d = Read("d", 60, 0, 0, 0, "", 0, 2)
    d.add_variant(1, 0, 30)
    d.add_variant(2, 0, 10)
    d.add_variant(3, 0, 50)
    e = Read("d", 60, 0, 0, 0, "", 1, 0)
    e.add_variant(1, 0, 30)
    e.add_variant(2, 0, 10)
    e.add_variant(3, 0, 50)
    expected_votes = {
        1: {(0, 0): 50, (0, 1): 0, (1, 1): 20, (1, 0): 0},
        2: {(0, 0): 10, (0, 1): 30},
    }
    votes = compute_votes({1: False, 2: False, 3: True}, [a, b, c])
    assert votes == expected_votes
