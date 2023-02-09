"""
Tests for 'whatshap stats'
"""
from collections import namedtuple
from whatshap.cli.stats import run_stats, unpack_chromosomes


def test_stats1(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_stats(
        vcf="tests/data/phased1.vcf",
        tsv=outtsv,
        sample="sample1",
        chr_lengths="tests/data/chr-lengths.txt",
    )
    with open(outtsv) as f:
        lines = [l.split("\t") for l in f]
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
    with open(outtsv) as f:
        lines = [l.split("\t") for l in f]
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


def test_overlapping_phaseblocks(tmp_path):
    """
    We have three phaseblocks on chrA which is 1000 bp

        chrA:100-700 --> 600 bp
        chrA:410-470 --> 60 bp
        chrA:800-950 --> 150 bp

    Total block sum should be 600 + 60 + 150 = 810 bp

    For NG50 the first block is split since the second block overlaps, now we have four blocks

        chrA:100-350 --> 250 bp
        chrA:410-470 --> 60 bp
        chrA:500-700 --> 200 bp
        chrA:800-950 --> 150 bp

    Total block sum is now 250 + 60 + 200 + 150 = 660 bp
    Half of the total length is 1000 * 0.5 = 500 bp.
    Let's calculate NG50 by adding block lengths in descending order until we exceed 500 bp

        block   length  total   >500
        1       250     250     no
        2       200     450     no
        3       150     600     yes ->  NG50 = 150 bp
    """

    outtsv = tmp_path / "output.tsv"
    run_stats(
        vcf="tests/data/phased_overlapping.vcf",
        tsv=outtsv,
        sample="sample1",
    )
    with open(outtsv) as f:
        lines = [l.split("\t") for l in f]
    assert len(lines) == 2
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry = Fields(*lines[1])

    assert entry.chromosome == "chrA"
    assert entry.blocks == "3"
    assert entry.bp_per_block_sum == "660"
    assert entry.block_n50[:-1] == "150"


def test_unphased_vcf(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_stats(
        vcf="tests/data/unphased.vcf",
        tsv=outtsv,
        sample="sample1",
        chromosomes=["chrA"],
    )
    with open(outtsv) as f:
        lines = [l.strip().split("\t") for l in f]
    assert len(lines) == 2
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry = Fields(*lines[1])

    assert entry.variants == "4"
    assert entry.heterozygous_variants == "3"
    assert entry.phased == "0"
    assert entry.unphased == "3"
    assert entry.bp_per_block_avg == "nan"
    assert entry.block_n50 == "nan"


def test_unpack_chromosomes_comma():
    unpack_comma = unpack_chromosomes(["chrA,chrB,chrC"])
    assert unpack_comma == ["chrA", "chrB", "chrC"]


def test_unpack_chromosomes_mixed():
    unpack_comma = unpack_chromosomes(["chrA,chrB,", "chrC"])
    assert unpack_comma == ["chrA", "chrB", "chrC"]
