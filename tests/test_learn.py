"""
Tests for whatshap learn module

"""

from whatshap.cli.learn import learn
from whatshap.core import Caller
import filecmp


def test_learn():
    expected = "tests/data/short-genome/learn-data/expected.txt"
    observed = "tests/data/short-genome/learn-data/observed.txt"
    learn(
        reference="tests/data/short-genome/learn-data/short_ref.fasta",
        bam="tests/data/short-genome/learn-data/short-reads.bam",
        vcf="tests/data/short-genome/learn-data/variant.vcf",
        kmer=7,
        window=25,
        output=observed,
    )
    assert filecmp.cmp(expected, observed, shallow=True)
