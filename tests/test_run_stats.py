"""
Tests for 'whatshap stats'
"""
from collections import namedtuple
from whatshap.cli.stats import run_stats


def test_stats1(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_stats(
        vcf="tests/data/phased1.vcf",
        tsv=outtsv,
        sample="sample1",
        chr_lengths="tests/data/chr-lengths.txt",
    )
    lines = [l.split("\t") for l in open(outtsv)]
    assert len(lines) == 4
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chrA, entry_chrB, entry_all = [Fields(*l) for l in lines[1:]]

    assert entry_chrA.chromosome == "chrA"
    assert entry_chrA.variants == "8"
    assert entry_chrA.phased == "7"
    assert entry_chrA.unphased == "1"
    assert entry_chrA.blocks == "2"
    assert entry_chrA.variant_per_block_sum == "7"
    assert entry_chrA.bp_per_block_sum == "551"
    assert entry_chrA.block_n50[:-1] == "101"

    assert entry_chrB.chromosome == "chrB"
    assert entry_chrB.variants == "2"
    assert entry_chrB.phased == "2"
    assert entry_chrB.unphased == "0"
    assert entry_chrB.blocks == "1"
    assert entry_chrB.bp_per_block_sum == "50"
    assert entry_chrB.variant_per_block_sum == "2"
    assert entry_chrB.block_n50[:-1] == "0"

    assert entry_all.chromosome == "ALL"
    assert entry_all.variants == "10"
    assert entry_all.phased == "9"
    assert entry_all.unphased == "1"
    assert entry_all.blocks == "3"
    assert entry_all.bp_per_block_sum == "601"
    assert entry_all.variant_per_block_sum == "9"
    assert entry_all.block_n50[:-1] == "0"


def test_stats2(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_stats(
        vcf="tests/data/phased3.vcf",
        tsv=outtsv,
        sample="sample1",
        chr_lengths="tests/data/chr-lengths.txt",
    )
    lines = [l.split("\t") for l in open(outtsv)]
    assert len(lines) == 4
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chrA, entry_chrB, entry_all = [Fields(*l) for l in lines[1:]]

    assert entry_chrA.chromosome == "chrA"
    assert entry_chrA.variants == "9"
    assert entry_chrA.phased == "4"
    assert entry_chrA.unphased == "5"
    assert entry_chrA.blocks == "1"
    assert entry_chrA.variant_per_block_sum == "4"
    assert entry_chrA.bp_per_block_sum == "350"
    assert entry_chrA.block_n50[:-1] == "0"

    assert entry_chrB.chromosome == "chrB"
    assert entry_chrB.variants == "4"
    assert entry_chrB.phased == "4"
    assert entry_chrB.unphased == "0"
    assert entry_chrB.blocks == "1"
    assert entry_chrB.variant_per_block_sum == "4"
    assert entry_chrB.bp_per_block_sum == "400"
    assert entry_chrB.block_n50[:-1] == "400"

    assert entry_all.chromosome == "ALL"
    assert entry_all.variants == "13"
    assert entry_all.phased == "8"
    assert entry_all.unphased == "5"
    assert entry_all.blocks == "2"
    assert entry_all.variant_per_block_sum == "8"
    assert entry_all.bp_per_block_sum == "750"
    assert entry_all.block_n50[:-1] == "350"
