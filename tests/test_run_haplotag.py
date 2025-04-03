import os
from collections import defaultdict
import shutil
import pysam
import pytest

from whatshap.cli.haplotag import run_haplotag, SupplementaryHaplotaggingStrategy
from whatshap.cli import CommandLineError
import argparse


def test_haplotag(tmp_path):
    outbam1 = tmp_path / "output1.bam"
    outbam2 = tmp_path / "output2.bam"
    outlist1 = tmp_path / "list1.tsv"
    outlist2 = tmp_path / "list2.tsv"

    # run haplotag with two vcfs containing opposite phasings (i.e. 1|0 - 0|1 ..)
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist1,
        output=outbam1,
    )
    run_haplotag(
        variant_file="tests/data/haplotag_2.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist2,
        output=outbam2,
    )
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP"):
            assert a2.has_tag("HP")
            assert a1.get_tag("HP") != a2.get_tag("HP")
    for n, (line1, line2) in enumerate(zip(open(outlist1), open(outlist2))):
        fields1 = line1.split(sep="\t")
        fields2 = line2.split(sep="\t")
        assert len(fields1) == len(fields2) == 4
        if n == 0:
            continue
        queryname1, haplotype1, phaseset1, chromosome1 = fields1
        queryname2, haplotype2, phaseset2, chromosome2 = fields2
        assert queryname1 == queryname2
        assert (haplotype1 == haplotype2 == "none") or (haplotype1 != haplotype2)
        assert chromosome1 == chromosome2
    assert n == 20


@pytest.mark.parametrize(
    "vcf_path",
    [
        "tests/data/haplotag_2.vcf.gz",
        "tests/data/haplotag_with_csi_index.vcf.gz",
        "tests/data/haplotag_2.bcf",
    ],
)
def test_haplotag2(tmp_path, vcf_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(variant_file=vcf_path, alignment_file="tests/data/haplotag.bam", output=outbam)
    ps_count = 0
    for alignment in pysam.AlignmentFile(outbam):
        if alignment.has_tag("PS"):
            ps_count += 1
        if alignment.has_tag("HP"):
            # simulated bam, we know from which haplotype each read originated (given in read name)
            true_ht = int(alignment.query_name[-1])
            assert true_ht == alignment.get_tag("HP")
    assert ps_count > 0


def test_haplotag_fails_if_index_missing(tmp_path):
    outbam = tmp_path / "output.bam"
    vcf_path = tmp_path / "vcf_without_index.vcf.gz"
    shutil.copy("tests/data/haplotag_1.vcf.gz", vcf_path)
    with pytest.raises(CommandLineError):
        run_haplotag(variant_file=vcf_path, alignment_file="tests/data/haplotag.bam", output=outbam)


def test_haplotag_cli_parser(tmp_path):
    """
    This test captures an error in the parser of the cli haplotag
    module - a wrong default value of "[]" instead of "None" for the
    "--regions" option leads to an empty output
    :return:
    """
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    outbam = tmp_path / "output.bam"
    parser = argparse.ArgumentParser(description="haplotag_test_parser", prog="whatshap_pytest")
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [
            "--no-reference",
            "--output",
            str(outbam),
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
        ]
    )
    haplotag_args.reference = False
    del haplotag_args.no_reference
    run_haplotag(**vars(haplotag_args))
    ps_count = 0
    for alignment in pysam.AlignmentFile(outbam):
        if alignment.has_tag("PS"):
            ps_count += 1
        if alignment.has_tag("HP"):
            # simulated BAM, we know from which haplotype each read originated (given in read name)
            true_ht = int(alignment.query_name[-1])
            assert true_ht == alignment.get_tag("HP")
    assert ps_count > 0


@pytest.mark.parametrize(
    "supplementary_strategy_cli_flag",
    [
        ("", SupplementaryHaplotaggingStrategy.SKIP),
        ("--tag-supplementary", SupplementaryHaplotaggingStrategy.COPY_PRIMARY),
        ("--tag-supplementary=skip", SupplementaryHaplotaggingStrategy.SKIP),
        ("--tag-supplementary=copy-primary", SupplementaryHaplotaggingStrategy.COPY_PRIMARY),
        (
            "--tag-supplementary=independent-or-skip",
            SupplementaryHaplotaggingStrategy.INDEPENDENT_OR_SKIP,
        ),
        (
            "--tag-supplementary=independent-or-copy-primary",
            SupplementaryHaplotaggingStrategy.INDEPENDENT_OR_COPY_PRIMARY,
        ),
    ],
)
def test_haplotag_cli_parser_supplementary_strategy_strategy_cli_flag(
    tmp_path,
    supplementary_strategy_cli_flag,
):
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    parser = argparse.ArgumentParser()
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [x for x in supplementary_strategy_cli_flag[0].split("=") if len(x) > 0]
        + [
            "--no-reference",
            "--output",
            os.devnull,
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
        ]
    )
    assert haplotag_args.supplementary_strategy == supplementary_strategy_cli_flag[1]


@pytest.mark.parametrize(
    "supplementary_distance_cli_flag",
    [("", 100_000), ("--supplementary-distance=100", 100)],
)
def test_haplotag_cli_parser_supplementary_strategy_distance_cli_flag(
    tmp_path,
    supplementary_distance_cli_flag,
):
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    parser = argparse.ArgumentParser()
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [x for x in supplementary_distance_cli_flag[0].split("=") if len(x) > 0]
        + [
            "--no-reference",
            "--output",
            os.devnull,
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
        ]
    )
    assert haplotag_args.supplementary_distance_threshold == supplementary_distance_cli_flag[1]


@pytest.mark.parametrize(
    "supplementary_strands_cli_flag",
    [
        ("", True),
        ("--no-supplementary-strand-match", False),
    ],
)
def test_haplotag_cli_parser_supplementary_strategy_strands_cli_flag(
    tmp_path,
    supplementary_strands_cli_flag,
):
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    parser = argparse.ArgumentParser()
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [x for x in supplementary_strands_cli_flag[0].split("=") if len(x) > 0]
        + [
            "--no-reference",
            "--output",
            os.devnull,
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
        ]
    )
    assert haplotag_args.supplementary_strand_match == supplementary_strands_cli_flag[1]


def test_haplotag_cli_parser_supplementary_distance_threshold(tmp_path):
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    outbam = tmp_path / "output.bam"
    parser = argparse.ArgumentParser()
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [
            "--no-reference",
            "--output",
            str(outbam),
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
            "--tag-supplementary=skip",
        ]
    )
    assert haplotag_args.supplementary_distance_threshold == 100_000
    haplotag_args = parser.parse_args(
        [
            "--no-reference",
            "--output",
            str(outbam),
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
            "--tag-supplementary=skip",
            "--supplementary-distance=100",
        ]
    )
    assert haplotag_args.supplementary_distance_threshold == 100


def test_haplotag_cli_parser_supplementary_strand_match_requirement(tmp_path):
    from whatshap.cli.haplotag import add_arguments as haplotag_add_arguments

    outbam = tmp_path / "output.bam"
    parser = argparse.ArgumentParser()
    haplotag_add_arguments(parser)
    haplotag_args = parser.parse_args(
        [
            "--no-reference",
            "--output",
            str(outbam),
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
            "--tag-supplementary=skip",
        ]
    )
    assert haplotag_args.supplementary_strand_match
    haplotag_args = parser.parse_args(
        [
            "--no-reference",
            "--output",
            str(outbam),
            "tests/data/haplotag_2.vcf.gz",
            "tests/data/haplotag.bam",
            "--tag-supplementary=skip",
            "--supplementary-distance=100",
            "--no-supplementary-strand-match",
        ]
    )
    assert not haplotag_args.supplementary_strand_match


"""
For the following tests we cover cases of supplementary alignments haplotagging strategy
The idea is cover the use case of having a vcf produced and/or phased with long reads in a matching normal sample
and have derived tumor sample reads haplotagged with respective phased vcf.
Due to potential rearrangements in tumor, respective reads may have multiple supplementary alignments, that
fall into various "germline" phase blocks.

supplementary_strategy_test.grch38.bam -- alignment of 2 fake reads, imitating multy supplementary alignments of long reads
supplementary_strategy_test.grch38.vcf.gz -- phased snps that span the alignment regions of 2 reads in question

chr1_PS1
    region: chr1:17,985,758-17,997,194  (~11Kbp)
    PS_id: 16849384
chr1_PS1_sub
    region: chr1:17,986,527-17,989,576 (~3Kbp)
    PS_id:  16849384
chr1_NPS
    region: chr1:18,063,682-18,082,783 (~19Kbp)
    PS_id: NA
chr1_NP_sub
    region: chr1:18,071,841-18,074,275 (~2.5Kbp)
    PS_id: NA
chr1 PS2:
    region:  chr1:18,130,745-18,132,827 (~2Kbp)
    PS_id: 18103117

chr2_PS1:
    region: chr2:27,018,266-27,023,671 (~5Kbp)
    PS_id: 26802880
chr2_NPS:
    region: chr2:28,310,196-28,312,671 (~2Kbp)
    PS_id: NA
chr2_PS2:
    region: chr2:28,458,793-28,462,863 (~4Kbp)
    PS_id: 28342675

read R1 is represented by the following string of reference segment (rc == reverse complement; sub == subregion)
chr1_PS1_H1 -> chr1_NPS_H1 -> chr1_PS2_H1 -> chr1_rcPS1_H2 -> chr2_PS1_H1 -> chr1_rcPS1_H1

chr1_PS1_H1         -- supplem, cigar: 3442M4D7991M49443S, flag: 2048
chr1_NPS_H1         -- primary, cigar: 11433S6757M1D8551M6D3789M30346S, flag: 0.
chr1_PS2_H1         -- supplem, cigar: 30528S1424M4D655M28269S, flag: 2048
chr1_rcPS1_H2       -- supplem, cigar: 16839S3442M2D6222M5D1768M32605S, flag: 2064
chr2_PS1_H1         -- supplem, cigar: 44037S5406M11433S, flag: 2048
chr1_rcPS1_H1       -- supplem, cigar: 3442M4D7991M49443S, flag: 2064

for read R1 the primary alignment is the chr1_NPS_H1 segment that fall into a non-vcf-phased region and can't be assigned to any haplotype

read R2 is represented by the following string of reference segment (rc == reverse complement; sub == subregion)
chr1_PS1_H1 -> chr1_NPS_sub_H1 -> chr1_PS2_H2 -> chr1_rcPS1_sub_H2 -> chr2_PS1_H1 -> chr1_rPS1_sub_H1

chr1_PS1_H1         -- primary, cigar: 3442M4D7991M16014S, flag: 0
chr1_NPS_sub_H1     -- supplem, cigar: 11432S2445M13570S, flag: 2048
chr1_PS2_H2         -- supplem, cigar: 13868S1424M4D655M11500S, flag = 2048
chr1_rcPS1_sub_H2   -- supplem, cigar: 8452S2673M2D375M15947S, flag = 2064
chr2_PS1_H1         -- supplem, cigar: 18995S5407M3045S, flag = 2048
chr1_rPS1_sub_H1    -- supplem, cigar: 2673M4D373M24401S, flag = 2064
"""


def test_run_haplotag_supplementary_skip(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_default_strategy = tmp_path / "output.default_haplotag_strategy.bam"
    out_bam_explicit_skip = tmp_path / "output.explicit_skip_strategy.bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_default_strategy,
        ignore_read_groups=True,
    )

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_explicit_skip,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.SKIP,
    )
    pairs = []
    a1: pysam.AlignedSegment
    a2: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_default_strategy) as source1, pysam.AlignmentFile(
        out_bam_explicit_skip
    ) as source2:
        pairs = list(zip(source1, source2))
    for a1, a2 in pairs:
        assert a1.query_name == a2.query_name
        if a1.is_supplementary:
            assert not a1.has_tag("HP")
            assert not a2.has_tag("HP")
            assert not a1.has_tag("PS")
            assert not a2.has_tag("PS")
        if a1.query_name == "R1" and not a1.is_supplementary:
            assert not a1.has_tag("HP")
            assert not a2.has_tag("HP")
            assert not a1.has_tag("PS")
            assert not a2.has_tag("PS")
        elif a1.query_name == "R2" and not a1.is_supplementary:
            assert a1.get_tag("HP") == a2.get_tag("HP")
            assert a1.get_tag("HP") == 1
            assert a1.get_tag("PS") == a2.get_tag("PS")
            assert a1.get_tag("PS") == 16849384


def test_run_haplotag_supplementary_copy_primary_no_strand_match_permissive_distance(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_copy_primary_strategy = tmp_path / "output.copy_primary.bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_copy_primary_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.COPY_PRIMARY,
        supplementary_strand_match=False,
        supplementary_distance_threshold=1_000_000,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_copy_primary_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                assert not a.has_tag("HP")
                assert not a.has_tag("PS")
            if a.query_name == "R2":
                if a.reference_name == "chr2":
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                else:
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384


def test_run_haplotag_supplementary_copy_primary_strand_match_permissive_distance(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_copy_primary_strategy = tmp_path / "output..bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_copy_primary_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.COPY_PRIMARY,
        supplementary_strand_match=True,
        supplementary_distance_threshold=1_000_000,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_copy_primary_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                assert not a.has_tag("HP")
                assert not a.has_tag("PS")
            if a.query_name == "R2":
                if a.reference_name == "chr2":
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                elif a.flag == 2064:
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                else:
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384


def test_run_haplotag_supplementary_copy_primary_strand_match_small_distance(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_copy_primary_strategy = tmp_path / "output.bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_copy_primary_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.COPY_PRIMARY,
        supplementary_strand_match=True,
        supplementary_distance_threshold=100,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_copy_primary_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                assert not a.has_tag("HP")
                assert not a.has_tag("PS")
            if a.query_name == "R2":
                if a.is_supplementary:
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                else:
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384


def test_run_haplotag_supplementary_copy_primary_no_strand_match_small_distance(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_copy_primary_strategy = tmp_path / "output.bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_copy_primary_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.COPY_PRIMARY,
        supplementary_strand_match=False,
        supplementary_distance_threshold=100,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_copy_primary_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                assert not a.has_tag("HP")
                assert not a.has_tag("PS")
            if a.query_name == "R2":
                if a.is_supplementary and a.cigarstring not in [
                    "8452S2673M2D375M15947S",
                    "2673M4D373M24401S",
                ]:
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                else:
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384


def test_run_haplotag_supplementary_independent_or_skip(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_independent_or_skip_strategy = tmp_path / "output..bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_independent_or_skip_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.INDEPENDENT_OR_SKIP,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_independent_or_skip_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                # chr1_PS1_H1
                if (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M49443S"
                    and a.flag == 2048
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # chr1_NPS_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "11433S6757M1D8551M6D3789M30346S"
                    and a.flag == 0
                ):
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                # chr1_PS2_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "30528S1424M4D655M28269S"
                    and a.flag == 2048
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 18103117
                # chr1_rcPS1_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "16839S3442M2D6222M5D1768M32605S"
                    and a.flag == 2064
                ):
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 16849384
                # chr2_PS1_H1
                elif (
                    a.reference_name == "chr2"
                    and a.cigarstring == "44037S5406M11433S"
                    and a.flag == 2048
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 26802880
                # chr1_rcPS1_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M49443S"
                    and a.flag == 2064
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # we should not get here to R1, so a failsafe
                else:
                    assert False
            if a.query_name == "R2":
                # chr1_PS1_H1
                if (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M16014S"
                    and a.flag == 0
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # chr1_NPS_sub_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "11432S2445M13570S"
                    and a.flag == 2048
                ):
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                # chr1_PS2_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "13868S1424M4D655M11500S"
                    and a.flag == 2048
                ):
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 18103117
                # chr1_rcPS1_sub_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "8452S2673M2D375M15947S"
                    and a.flag == 2064
                ):
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 16849384
                # chr2_PS1_H1
                elif (
                    a.reference_name == "chr2"
                    and a.cigarstring == "18995S5407M3045S"
                    and a.flag == 2048
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 26802880
                # chr1_rPS1_sub_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "2673M4D373M24401S"
                    and a.flag == 2064
                ):
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # we should not get here to R2, so a failsafe
                else:
                    assert False


def test_run_haplotag_supplementary_independent_or_copy_primary(tmp_path):
    var_file = "tests/data/supplementary_strategy_test.grch38.vcf.gz"
    alignment_file = "tests/data/supplementary_strategy_test.grch38.bam"

    out_bam_independent_or_copy_primary_strategy = tmp_path / "output..bam"

    run_haplotag(
        variant_file=var_file,
        alignment_file=alignment_file,
        output=out_bam_independent_or_copy_primary_strategy,
        ignore_read_groups=True,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.INDEPENDENT_OR_COPY_PRIMARY,
    )

    a: pysam.AlignedSegment
    with pysam.AlignmentFile(out_bam_independent_or_copy_primary_strategy) as source:
        for a in source:
            if a.query_name == "R1":
                # chr1_PS1_H1
                if (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M49443S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # chr1_NPS_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "11433S6757M1D8551M6D3789M30346S"
                    and a.flag == 0
                ):
                    assert not a.has_tag("HP")
                    assert not a.has_tag("PS")
                # chr1_PS2_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "30528S1424M4D655M28269S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 18103117
                # chr1_rcPS1_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "16839S3442M2D6222M5D1768M32605S"
                    and a.flag == 2064
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 16849384
                # chr2_PS1_H1
                elif (
                    a.reference_name == "chr2"
                    and a.cigarstring == "44037S5406M11433S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 26802880
                # chr1_rcPS1_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M49443S"
                    and a.flag == 2064
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # we should not get here to R1, so a failsafe
                else:
                    assert False
            if a.query_name == "R2":
                # chr1_PS1_H1
                if (
                    a.reference_name == "chr1"
                    and a.cigarstring == "3442M4D7991M16014S"
                    and a.flag == 0
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # chr1_NPS_sub_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "11432S2445M13570S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # chr1_PS2_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "13868S1424M4D655M11500S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 18103117
                # chr1_rcPS1_sub_H2
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "8452S2673M2D375M15947S"
                    and a.flag == 2064
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 2
                    assert a.get_tag("PS") == 16849384
                # chr2_PS1_H1
                elif (
                    a.reference_name == "chr2"
                    and a.cigarstring == "18995S5407M3045S"
                    and a.flag == 2048
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 26802880
                # chr1_rPS1_sub_H1
                elif (
                    a.reference_name == "chr1"
                    and a.cigarstring == "2673M4D373M24401S"
                    and a.flag == 2064
                ):
                    assert a.has_tag("HP")
                    assert a.has_tag("PS")
                    assert a.get_tag("HP") == 1
                    assert a.get_tag("PS") == 16849384
                # we should not get here to R2, so a failsafe
                else:
                    assert False


def test_haplotag_missing_SM_tag(tmp_path):
    outbam1 = tmp_path / "output1.bam"
    outbam2 = tmp_path / "output2.bam"
    outlist1 = tmp_path / "list1.tsv"
    outlist2 = tmp_path / "list2.tsv"

    # run haplotag
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist1,
        output=outbam1,
        ignore_read_groups=True,
    )
    # use copy of 'haplotag.bam' which lacks the 'SM' tag
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag_noSM.bam",
        haplotag_list=outlist2,
        output=outbam2,
        ignore_read_groups=True,
    )

    # results should be identical
    a1: pysam.AlignedSegment
    a2: pysam.AlignedSegment
    with pysam.AlignmentFile(outbam1) as source1, pysam.AlignmentFile(outbam2) as source2:
        for a1, a2 in zip(source1, source2):
            assert a1.query_name == a2.query_name
            if a1.has_tag("HP"):
                assert a1.get_tag("HP") == a2.get_tag("HP")
            for n, (line1, line2) in enumerate(zip(open(outlist1), open(outlist2))):
                fields1 = line1.split(sep="\t")
                fields2 = line2.split(sep="\t")
                assert len(fields1) == len(fields2) == 4
                if n == 0:
                    continue
                queryname1, haplotype1, phaseset1, chromosome1 = fields1
                queryname2, haplotype2, phaseset2, chromosome2 = fields2
                assert queryname1 == queryname2
                assert haplotype1 == haplotype2
                assert chromosome1 == chromosome2
            assert n == 20


def test_haplotag_missing_chromosome(tmp_path):
    outbam = tmp_path / "output.bam"

    # input BAM contains a chromosome for which there is no variant in the input VCF
    run_haplotag(
        variant_file="tests/data/haplotag.missing_chr.vcf.gz",
        alignment_file="tests/data/haplotag.large.bam",
        output=outbam,
    )
    ps_count = 0
    for alignment in pysam.AlignmentFile(outbam):
        if alignment.has_tag("PS"):
            ps_count += 1
    assert ps_count > 0


def test_contig_exists_in_bam_but_not_in_vcf_header(tmp_path):
    outbam = tmp_path / "output.bam"

    with pytest.raises(CommandLineError) as e:
        run_haplotag(
            variant_file="tests/data/haplotag.without_chr2.vcf.gz",
            alignment_file="tests/data/haplotag.large.bam",  # has reads mapped to chr2
            output=outbam,
        )
    assert "contig does not exist" in e.value.args[0]

    run_haplotag(
        variant_file="tests/data/haplotag.without_chr2.vcf.gz",
        alignment_file="tests/data/haplotag.large.bam",  # has reads mapped to chr2
        output=outbam,
        skip_missing_contigs=True,
    )


def test_haplotag_no_readgroups1(tmp_path):
    outbam1 = tmp_path / "output1.bam"
    outbam2 = tmp_path / "output2.bam"

    # run haplotag with/without --ignore-read-groups, results should be identical since files contain only data for one sample
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        output=outbam1,
    )
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag_noRG.bam",
        output=outbam2,
        ignore_read_groups=True,
    )
    count = 0
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP"):
            assert a1.get_tag("HP") == a2.get_tag("HP")
            count += 1
    assert count > 0


def test_haplotag_no_readgroups2():
    with pytest.raises((CommandLineError, ValueError)):
        # VCF contains multiple samples, there should be an error
        run_haplotag(
            variant_file="tests/data/haplotag_noRG.vcf.gz",
            alignment_file="tests/data/haplotag_noRG.bam",
            output="/dev/null",
            ignore_read_groups=True,
        )


def test_haplotag_sample_given(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag_sample.vcf.gz",
        alignment_file="tests/data/haplotag_sample.bam",
        given_samples=["mother"],
        output=outbam,
    )
    for alignment in pysam.AlignmentFile(outbam):
        if alignment.get_tag("RG") == "mother":
            assert alignment.has_tag("HP")
        else:
            assert not alignment.has_tag("HP")


def haplotag_different_sorting(tmp_path):
    outbam1 = tmp_path / "output1.bam"
    outbam2 = tmp_path / "output2.bam"

    # both VCFs contain the same positions, but chromosomes are sorted differently
    run_haplotag(
        variant_file="tests/data/haplotag.large.vcf.gz",
        alignment_file="tests/data/haplotag.large.bam",
        output=outbam1,
    )
    run_haplotag(
        variant_file="tests/data/haplotag.large.2.vcf.gz",
        alignment_file="tests/data/haplotag.large.bam",
        output=outbam2,
    )
    count = 0
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP"):
            assert a1.get_tag("HP") == a2.get_tag("HP")
            count += 1
    assert count > 0


def test_haplotag_10X(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag.10X.vcf.gz",
        alignment_file="tests/data/haplotag.10X.bam",
        output=outbam,
    )
    # map BX tag --> readlist
    bx_tag_to_readlist = defaultdict(list)
    for alignment in pysam.AlignmentFile(outbam):
        if alignment.has_tag("BX") and alignment.has_tag("HP"):
            bx_tag_to_readlist[alignment.get_tag("BX")].append(alignment)
    # reads having same BX tag need to be assigned to same haplotype
    for tag in bx_tag_to_readlist.keys():
        haplotype = bx_tag_to_readlist[tag][0].get_tag("HP")
        for read in bx_tag_to_readlist[tag]:
            assert haplotype == read.get_tag("HP")


def test_haplotag_10X_2(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag.10X_2.vcf.gz",
        alignment_file="tests/data/haplotag.10X.bam",
        output=outbam,
    )
    count = 0
    for a1, a2 in zip(
        pysam.AlignmentFile("tests/data/haplotag.10X.bam"), pysam.AlignmentFile(outbam)
    ):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP") and a2.has_tag("HP"):
            assert a1.get_tag("HP") == a2.get_tag("HP")
            count += 1
    assert count > 0


def test_haplotag_10X_ignore_linked_read(tmp_path):
    outbam_links = tmp_path / "with_links.bam"
    outbam_nolinks = tmp_path / "no_links.bam"
    run_haplotag(
        variant_file="tests/data/haplotag.10X.vcf.gz",
        alignment_file="tests/data/haplotag.10X_3.bam",
        output=outbam_links,
    )
    run_haplotag(
        variant_file="tests/data/haplotag.10X.vcf.gz",
        alignment_file="tests/data/haplotag.10X_3.bam",
        output=outbam_nolinks,
        ignore_linked_read=True,
    )
    expected_links = {"read1": [1, 4], "read2": [1, 4], "read3": [1, 11], "read4": [1, 11]}
    expected_no_links = {"read1": [2, 66], "read2": [1, 70], "read3": [2, 55], "read4": [1, 66]}
    for a1, a2 in zip(pysam.AlignmentFile(outbam_links), pysam.AlignmentFile(outbam_nolinks)):
        assert a1.query_name == a2.query_name
        name = a1.query_name
        if name == "read5":
            # read5 assigned according to other reads with same BX tag
            assert a1.has_tag("HP")
            assert a1.get_tag("HP") == 1
            # using --ignore-linked-read, read5 must be untagged
            assert not a2.has_tag("HP")
        else:
            assert a1.get_tag("HP") == expected_links[name][0]
            assert a1.get_tag("PC") == expected_links[name][1]
            assert a2.get_tag("HP") == expected_no_links[name][0]
            assert a2.get_tag("PC") == expected_no_links[name][1]


def test_haplotag_supplementary(tmp_path):
    # test --tag-supplementary option which assigns supplementary
    # reads to haplotypes based on the tag of their primary alignment.
    outbam1 = tmp_path / "supp-untagged.bam"
    outbam2 = tmp_path / "supp-tagged.bam"
    run_haplotag(
        variant_file="tests/data/haplotag.supplementary.vcf.gz",
        alignment_file="tests/data/haplotag.supplementary.bam",
        output=outbam1,
        ignore_read_groups=True,
    )
    run_haplotag(
        variant_file="tests/data/haplotag.supplementary.vcf.gz",
        alignment_file="tests/data/haplotag.supplementary.bam",
        output=outbam2,
        supplementary_strategy=SupplementaryHaplotaggingStrategy.COPY_PRIMARY,
        ignore_read_groups=True,
        supplementary_strand_match=False,
        supplementary_distance_threshold=1_000_000_000,
    )
    # map name->haplotype
    primary_to_tag = {}
    supplementary_to_tag = {}
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP") and a2.has_tag("HP"):
            assert a1.get_tag("HP") == a2.get_tag("HP")
            assert not a1.is_supplementary
        if a2.has_tag("HP"):
            tag = a2.get_tag("HP")
            if a2.is_supplementary:
                supplementary_to_tag[a2.query_name] = tag
            else:
                primary_to_tag[a2.query_name] = tag
    # check if supplementary and primary tags agree
    assert len(primary_to_tag.keys()) == len(supplementary_to_tag.keys()) == 3
    for r, t in supplementary_to_tag.items():
        assert r in primary_to_tag
        primary_tag = primary_to_tag[r]
        assert t == primary_tag


def test_haplotag_regions(tmp_path):
    outbam1 = tmp_path / "output1.bam"
    outbam2 = tmp_path / "output2.bam"
    outlist1 = tmp_path / "list1.tsv"
    outlist2 = tmp_path / "list2.tsv"

    # run haplotag with identical VCF, but once specifying regions
    # output must be identical
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist1,
        output=outbam1,
        regions=None,
    )
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist2,
        output=outbam2,
        regions=["chr1"],
    )
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP"):
            assert a2.has_tag("HP")
            assert a1.get_tag("HP") == a2.get_tag("HP")
    for n, (line1, line2) in enumerate(zip(open(outlist1), open(outlist2))):
        assert line1 == line2
    assert n == 20


def test_haplotag_nonexisting_region():
    with pytest.raises(ValueError):
        run_haplotag(
            variant_file="tests/data/haplotag_1.vcf.gz",
            alignment_file="tests/data/haplotag.bam",
            haplotag_list=None,
            output=None,
            regions=["chr2"],
        )


def test_haplotag_selected_regions(tmp_path):
    start1 = 1054025
    end1 = 1069500
    start2 = 1075700
    outbam = tmp_path / "output.bam"
    outlist = tmp_path / "haplolist.tsv"
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        haplotag_list=outlist,
        output=outbam,
        regions=["chr1:{}-{}".format(start1, end1), "chr1:{}".format(start2)],
    )

    var_region1 = set()
    var_region2 = set()
    unphased_variants = [1074910, 1075707, 1075715]
    with pysam.VariantFile("tests/data/haplotag_1.vcf.gz", "rb") as vcf:
        for variant in vcf:
            if variant.pos in unphased_variants:
                continue
            if start1 <= variant.start <= end1:
                var_region1.add(variant.start)
            elif start2 <= variant.start:
                var_region2.add(variant.start)
            else:
                pass
    # sanity check:
    # there are no variants in the VCF
    # overlapping region 1
    assert len(var_region1) == 0

    with pysam.AlignmentFile(outbam, "rb") as test_bam:
        # Since not all variants from the VCF are selected,
        # count how many variants are overlapping the read.
        # If more than 1 overlap, read must be phased / have HP tag
        for aln in test_bam:
            num_ovl = sum([int(aln.reference_start <= v <= aln.reference_end) for v in var_region2])
            if num_ovl > 1:
                assert aln.has_tag("HP")


def test_cram_output(tmp_path):
    outcram = tmp_path / "output.cram"
    run_haplotag(
        variant_file="tests/data/pacbio/phased.vcf.gz",
        alignment_file="tests/data/pacbio/pacbio.bam",
        reference="tests/data/pacbio/reference.fasta",
        output=outcram,
    )
    with pysam.AlignmentFile(outcram) as f:
        assert f.is_cram


def test_haplotag_unmapped_reads(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag.10X.vcf.gz",
        alignment_file="tests/data/unmapped.bam",
        output=outbam,
    )
    pysam.index(str(outbam))
    with pysam.AlignmentFile(outbam) as af:
        alignments = list(af.fetch(until_eof=True))
    assert len(alignments) == 6
    assert not alignments[4].is_unmapped
    assert alignments[5].is_unmapped


def test_haplotag_triploid(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag_triploid.vcf.gz",
        alignment_file="tests/data/haplotag_triploid.bam",
        ploidy=3,
        output=outbam,
    )

    # manually computed haplotag scores and haplotype assignments
    readname_to_score = {
        "S1_31286_NA19240_HAP2": 23,
        "S1_248595_HG00514_HAP1": 18,
        "S1_103518_HG00514_HAP2": 29,
    }
    readname_to_haplotype = {
        "S1_31286_NA19240_HAP2": 3,
        "S1_248595_HG00514_HAP1": 1,
        "S1_103518_HG00514_HAP2": 2,
    }
    count = 0
    with pysam.AlignmentFile(outbam) as af:
        for alignment in af:
            count += 1
            assert readname_to_score[alignment.query_name] == alignment.get_tag("PC")
            assert readname_to_haplotype[alignment.query_name] == alignment.get_tag("HP")
    assert count == 3


def test_haplotag_tetraploid(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag_poly.vcf.gz",
        alignment_file="tests/data/haplotag_poly.bam",
        ploidy=4,
        output=outbam,
    )

    # manually computed haplotag scores and haplotype assignments
    readname_to_score = {
        "S1_31286_NA19240_HAP2": 6,
        "S1_248595_HG00514_HAP1": 4,
        "S1_284251_NA19240_HAP1": 14,
        "S1_103518_HG00514_HAP2": 16,
        "chr1:2000000-2000099": None,
    }
    readname_to_haplotype = {
        "S1_31286_NA19240_HAP2": 4,
        "S1_248595_HG00514_HAP1": 1,
        "S1_284251_NA19240_HAP1": 3,
        "S1_103518_HG00514_HAP2": 2,
        "chr1:2000000-2000099": None,
    }
    count = 0
    with pysam.AlignmentFile(outbam) as af:
        for alignment in af:
            count += 1
            score = readname_to_score[alignment.query_name]
            if score is not None:
                assert score == alignment.get_tag("PC")
                assert readname_to_haplotype[alignment.query_name] == alignment.get_tag("HP")
    assert count == 5


def test_haplotag_duplicates_are_tagged(tmp_path):
    # Create a version of the BAM file where all reads are marked as duplicates
    inbam_dup = tmp_path / "haplotag-duplicates.bam"
    with pysam.AlignmentFile("tests/data/haplotag.bam") as infile:
        with pysam.AlignmentFile(inbam_dup, mode="wb", template=infile) as outfile:
            for record in infile:
                record.is_duplicate = True
                outfile.write(record)
    pysam.index(str(inbam_dup))
    outbam_dup = tmp_path / "output-nodup.bam"
    outbam_nodup = tmp_path / "output-dup.bam"

    # Run haplotag twice and compare results
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file="tests/data/haplotag.bam",
        output=outbam_nodup,
    )
    run_haplotag(
        variant_file="tests/data/haplotag_1.vcf.gz",
        alignment_file=inbam_dup,
        output=outbam_dup,
    )
    count = 0
    for r1, r2 in zip(pysam.AlignmentFile(outbam_nodup), pysam.AlignmentFile(outbam_dup)):
        assert r1.query_name == r2.query_name
        if r1.has_tag("PS"):
            assert r2.has_tag("PS")
            assert r1.get_tag("PS") == r2.get_tag("PS")
            count += 1
    assert count > 0


def test_haplotag_run_twice(tmp_path):
    outbam = tmp_path / "output.bam"
    run_haplotag(
        variant_file="tests/data/haplotag_sample.vcf.gz",
        alignment_file="tests/data/haplotag_sample.bam",
        given_samples=["mother"],
        output=outbam,
    )
    # Index bam file
    pysam.index(str(outbam))

    outbam2 = tmp_path / "output2.bam"
    run_haplotag(
        variant_file="tests/data/haplotag_sample.vcf.gz",
        alignment_file=outbam,
        given_samples=["mother"],
        output=outbam2,
    )

    # Check that there are two PG unique entries for whatshap
    with pysam.AlignmentFile(outbam2) as f:
        pg_entries = f.header.get("PG")
        whatshap_ids = [entry["ID"] for entry in pg_entries if entry["ID"].startswith("whatshap")]
        assert len(whatshap_ids) == 2
        assert len(set(whatshap_ids)) == 2
