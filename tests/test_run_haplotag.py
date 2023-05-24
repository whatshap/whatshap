from collections import defaultdict
import shutil
import pysam
import pytest

from whatshap.cli.haplotag import run_haplotag
from whatshap.cli import CommandLineError


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
    import argparse as argp

    outbam = tmp_path / "output.bam"
    parser = argp.ArgumentParser(description="haplotag_test_parser", prog="whatshap_pytest")
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
    for a1, a2 in zip(pysam.AlignmentFile(outbam1), pysam.AlignmentFile(outbam2)):
        assert a1.query_name == a2.query_name
        if a1.has_tag("HP"):
            assert a2.has_tag("HP")
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
            assert a2.has_tag("HP")
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
            assert a2.has_tag("HP")
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
        tag_supplementary=True,
        ignore_read_groups=True,
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
    }
    readname_to_haplotype = {
        "S1_31286_NA19240_HAP2": 4,
        "S1_248595_HG00514_HAP1": 1,
        "S1_284251_NA19240_HAP1": 3,
        "S1_103518_HG00514_HAP2": 2,
    }
    count = 0
    with pysam.AlignmentFile(outbam) as af:
        for alignment in af:
            count += 1
            assert readname_to_score[alignment.query_name] == alignment.get_tag("PC")
            assert readname_to_haplotype[alignment.query_name] == alignment.get_tag("HP")
    assert count == 4


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
