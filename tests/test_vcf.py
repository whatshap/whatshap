import math

from pytest import raises, approx, fixture
from whatshap.core import PhredGenotypeLikelihoods, Genotype
from whatshap.cli.phase import run_whatshap
from whatshap.vcf import (
    VcfReader,
    MixedPhasingError,
    PloidyError,
    VariantCallPhase,
    BiallelicVcfVariant,
    GenotypeLikelihoods,
    VcfIndexMissing,
)
from whatshap.testhelpers import (
    canonic_index_to_biallelic_gt,
    canonic_index_list_to_biallelic_gt_list,
)

import pysam


@fixture(params=["whatshap", "hapchat"])
def algorithm(request):
    return request.param


def test_read_phased():
    tables = list(VcfReader("tests/data/phasedinput.vcf", phases=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "ref"
    assert table.samples == ["sample"]
    assert len(table.variants) == 2
    assert table.variants[0].reference_allele == "A"
    assert table.variants[0].alternative_allele == "C"
    assert table.variants[1].reference_allele == "G"
    assert table.variants[1].alternative_allele == "T"
    assert table.genotypes[0][0] == table.genotypes[0][1] == canonic_index_to_biallelic_gt(1)


def test_read_multisample_vcf():
    tables = list(VcfReader("tests/data/multisample.vcf"))
    assert len(tables) == 2
    table, table_b = tables
    assert table_b.chromosome == "chrB"
    assert table_b.samples == ["sample1", "sample2"]

    assert table.chromosome == "chrA"
    assert len(table.variants) == 3
    assert table.samples == ["sample1", "sample2"]

    assert table.variants[0].reference_allele == "A"
    assert table.variants[0].alternative_allele == "T"
    assert table.variants[1].reference_allele == "C"
    assert table.variants[1].alternative_allele == "G"
    assert table.variants[2].reference_allele == "G"
    assert table.variants[2].alternative_allele == "T"

    assert len(table.genotypes) == 2
    assert list(table.genotypes[0]) == canonic_index_list_to_biallelic_gt_list([1, 1, 1])
    assert list(table.genotypes[1]) == canonic_index_list_to_biallelic_gt_list([1, 1, 0])

    assert list(table.genotypes_of("sample1")) == canonic_index_list_to_biallelic_gt_list([1, 1, 1])
    assert list(table.genotypes_of("sample2")) == canonic_index_list_to_biallelic_gt_list([1, 1, 0])


def test_read_phased_vcf():
    for filename in ["tests/data/phased-via-HP.vcf", "tests/data/phased-via-PS.vcf"]:
        print("Testing", filename)
        tables = list(VcfReader(filename, phases=True))
        assert len(tables) == 2
        table_a, table_b = tables

        assert table_a.chromosome == "chrA"
        assert len(table_a.variants) == 4
        assert table_a.samples == ["sample1", "sample2"]

        assert table_b.chromosome == "chrB"
        assert len(table_b.variants) == 2
        assert table_b.samples == ["sample1", "sample2"]

        assert len(table_a.genotypes) == 2
        assert list(table_a.genotypes[0]) == canonic_index_list_to_biallelic_gt_list([1, 2, 1, 1])
        assert list(table_a.genotypes[1]) == canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1])
        assert list(table_a.genotypes_of("sample1")) == canonic_index_list_to_biallelic_gt_list(
            [1, 2, 1, 1]
        )
        assert list(table_a.genotypes_of("sample2")) == canonic_index_list_to_biallelic_gt_list(
            [1, 1, 1, 1]
        )

        assert len(table_b.genotypes) == 2
        assert list(table_b.genotypes[0]) == canonic_index_list_to_biallelic_gt_list([0, 1])
        assert list(table_b.genotypes[1]) == canonic_index_list_to_biallelic_gt_list([1, 2])
        assert list(table_b.genotypes_of("sample1")) == canonic_index_list_to_biallelic_gt_list(
            [0, 1]
        )
        assert list(table_b.genotypes_of("sample2")) == canonic_index_list_to_biallelic_gt_list(
            [1, 2]
        )

        print(table_a.phases)
        assert len(table_a.phases) == 2
        expected_phase_sample1 = [
            None,
            None,
            VariantCallPhase(block_id=300, phase=(1, 0), quality=23),
            VariantCallPhase(block_id=300, phase=(0, 1), quality=42),
        ]
        expected_phase_sample2 = [
            VariantCallPhase(block_id=100, phase=(0, 1), quality=10),
            VariantCallPhase(block_id=100, phase=(1, 0), quality=20),
            VariantCallPhase(block_id=300, phase=(0, 1), quality=30),
            VariantCallPhase(block_id=300, phase=(0, 1), quality=None),
        ]
        assert list(table_a.phases[0]) == expected_phase_sample1
        assert list(table_a.phases[1]) == expected_phase_sample2
        assert list(table_a.phases_of("sample1")) == expected_phase_sample1
        assert list(table_a.phases_of("sample2")) == expected_phase_sample2

        assert len(table_b.phases) == 2
        assert list(table_b.phases[0]) == [None, None]
        assert list(table_b.phases[1]) == [None, None]
        assert list(table_b.phases_of("sample1")) == [None, None]
        assert list(table_b.phases_of("sample2")) == [None, None]


def test_mixed_phasing_vcf():
    with raises(MixedPhasingError):
        list(VcfReader("tests/data/phased-via-mixed-HP-PS.vcf", phases=True))


def test_vcf_variant_hashability():
    v = [
        BiallelicVcfVariant(10, "A", "TC"),
        BiallelicVcfVariant(10, "A", "TCA"),
        BiallelicVcfVariant(10, "C", "TC"),
        BiallelicVcfVariant(20, "A", "TC"),
        BiallelicVcfVariant(10, "A", "TCA"),
        BiallelicVcfVariant(20, "A", "TC"),
    ]
    assert len(set(v)) == 4


def test_phasing_to_reads():
    for filename in ["tests/data/phased-via-HP.vcf", "tests/data/phased-via-PS.vcf"]:
        tables = list(VcfReader(filename, phases=True))
        assert len(tables) == 2
        table_a, table_b = tables
        phase_reads_sample1 = list(
            table_a.phased_blocks_as_reads(
                "sample1", table_a.variants, 17, 18, default_quality=90, mapq=101
            )
        )
        print(phase_reads_sample1)
        assert len(phase_reads_sample1) == 2
        read1, read2 = phase_reads_sample1[0], phase_reads_sample1[1]
        assert len(read1) == len(read2) == 2
        assert read1.name == "sample1_phase_0_block_300"
        assert read2.name == "sample1_phase_1_block_300"
        assert read1.source_id == read2.source_id == 17
        assert read1.mapqs == read2.mapqs == (101,)
        assert read1[0].position == read2[0].position == 300 - 1
        assert read1[0].allele == 1 != read2[0].allele
        assert read1[0].quality == read2[0].quality == 23
        assert read1[1].position == read2[1].position == 350 - 1
        assert read1[1].allele == 0 != read2[1].allele
        assert read1[1].quality == read2[1].quality == 42

        phase_reads_sample2 = list(
            table_a.phased_blocks_as_reads(
                "sample2", table_a.variants, 11, 12, default_quality=91, mapq=102
            )
        )
        print(phase_reads_sample2)
        assert len(phase_reads_sample2) == 4
        read1, _, read2, _ = phase_reads_sample2
        assert len(read1) == len(read2) == 2
        if read1[0].position > read2[0].position:
            read1, read2 = read2, read1
        assert read1.name == "sample2_phase_0_block_100"
        assert read1.source_id == 11
        assert read1.mapqs == (102,)
        assert read1[0].position == 100 - 1
        assert read1[0].allele == 0
        assert read1[0].quality == 10
        assert read1[1].position == 150 - 1
        assert read1[1].allele == 1
        assert read1[1].quality == 20
        assert read2.name == "sample2_phase_0_block_300"
        assert read2.source_id == 11
        assert read2.mapqs == (102,)
        assert read2[0].position == 300 - 1
        assert read2[0].allele == 0
        assert read2[0].quality == 30
        assert read2[1].position == 350 - 1
        assert read2[1].allele == 0
        assert read2[1].quality == 91

        variants = [
            BiallelicVcfVariant(350 - 1, "G", "T"),
            BiallelicVcfVariant(300 - 1, "G", "T"),
            BiallelicVcfVariant(17, "A", "TTC"),
            BiallelicVcfVariant(1000, "C", "G"),
        ]
        phase_reads_sample2 = list(
            table_a.phased_blocks_as_reads(
                "sample2", variants, 11, 12, default_quality=91, mapq=102
            )
        )
        print(phase_reads_sample2)
        assert len(phase_reads_sample2) == 2
        read = phase_reads_sample2[0]
        assert len(read) == 2
        assert read.name == "sample2_phase_0_block_300"
        assert read.source_id == 11
        assert read.mapqs == (102,)
        assert read[0].position == 300 - 1
        assert read[0].allele == 0
        assert read[0].quality == 30
        assert read[1].position == 350 - 1
        assert read[1].allele == 0
        assert read[1].quality == 91


def test_phasing_to_reads_polyploid():
    for filename in [
        "tests/data/phased-via-HP-polyploid.vcf",
        "tests/data/phased-via-PS-polyploid.vcf",
    ]:
        tables = list(VcfReader(filename, phases=True, mav=True))
        assert len(tables) == 2
        table_a, table_b = tables
        reads = list(
            table_a.phased_blocks_as_reads(
                "sample1", table_a.variants, 17, 18, default_quality=90, mapq=101, target_ploidy=4
            )
        )
        print(reads)
        assert len(reads) == 4
        assert all(len(read) == 2 for read in reads)
        assert reads[0].name == "sample1_phase_0_block_300"
        assert reads[3].name == "sample1_phase_3_block_300"
        assert all(read.source_id == 17 for read in reads)
        assert all(read.mapqs == (101,) for read in reads)
        assert all(read[0].position == 300 - 1 for read in reads)
        assert all(read[0].quality == 23 for read in reads)
        assert reads[0][0].allele == 0
        assert reads[1][0].allele == 0
        assert reads[2][0].allele == 1
        assert reads[3][0].allele == 1
        assert all(read[1].position == 350 - 1 for read in reads)
        assert all(read[1].quality == 42 for read in reads)
        assert reads[0][1].allele == 0
        assert reads[1][1].allele == 0
        assert reads[2][1].allele == 1
        assert reads[3][1].allele == 0

        reads = list(
            table_a.phased_blocks_as_reads(
                "sample2", table_a.variants, 11, 12, default_quality=91, mapq=102, target_ploidy=4
            )
        )
        print(reads)
        assert len(reads) == 8
        assert all(len(read) == 2 for read in reads)
        assert reads[0].name == "sample2_phase_0_block_100"
        assert reads[3].name == "sample2_phase_3_block_100"
        assert all(read.source_id == 11 for read in reads)
        assert all(read.mapqs == (102,) for read in reads)
        assert all(read[0].position == 100 - 1 for read in reads[:4])
        assert all(read[0].quality == 10 for read in reads[:4])
        assert reads[0][0].allele == 0
        assert reads[1][0].allele == 0
        assert reads[2][0].allele == 1
        assert reads[3][0].allele == 1
        assert all(read[1].position == 150 - 1 for read in reads[:4])
        assert all(read[1].quality == 20 for read in reads[:4])
        assert reads[0][1].allele == 1
        assert reads[1][1].allele == 0
        assert reads[2][1].allele == 1
        assert reads[3][1].allele == 0
        assert all(read[0].position == 300 - 1 for read in reads[4:])
        assert all(read[0].quality == 30 for read in reads[4:])
        assert reads[4][0].allele == 0
        assert reads[5][0].allele == 0
        assert reads[6][0].allele == 0
        assert reads[7][0].allele == 1
        assert all(read[1].position == 350 - 1 for read in reads[4:])
        assert all(read[1].quality == 91 for read in reads[4:])
        assert reads[4][1].allele == 1
        assert reads[5][1].allele == 0
        assert reads[6][1].allele == 2
        assert reads[7][1].allele == 1


def test_unknown_genotype():
    """VCF with './.' genotype"""
    tables = list(VcfReader("tests/data/unknown-genotype.vcf"))
    assert tables[0].genotypes[1][0] == Genotype([])
    assert tables[0].genotypes[1][0].is_none()


def test_normalize():
    assert BiallelicVcfVariant(100, "A", "C").normalized() == BiallelicVcfVariant(100, "A", "C")
    assert BiallelicVcfVariant(100, "", "A").normalized() == BiallelicVcfVariant(100, "", "A")
    assert BiallelicVcfVariant(100, "A", "").normalized() == BiallelicVcfVariant(100, "A", "")
    assert BiallelicVcfVariant(100, "A", "AC").normalized() == BiallelicVcfVariant(101, "", "C")
    assert BiallelicVcfVariant(100, "AC", "A").normalized() == BiallelicVcfVariant(101, "C", "")
    assert BiallelicVcfVariant(100, "ACAGACC", "ACAGACT").normalized() == BiallelicVcfVariant(
        106, "C", "T"
    )
    assert BiallelicVcfVariant(100, "GCTG", "GCTAAA").normalized() == BiallelicVcfVariant(
        103, "G", "AAA"
    )
    assert BiallelicVcfVariant(100, "ATTA", "ATA").normalized() == BiallelicVcfVariant(101, "T", "")
    assert BiallelicVcfVariant(100, "ATTTC", "ATTTTTTC").normalized() == BiallelicVcfVariant(
        101, "", "TTT"
    )
    assert BiallelicVcfVariant(100, "GCTGTT", "GCTAAATT").normalized() == BiallelicVcfVariant(
        103, "G", "AAA"
    )


def test_read_duplicate_position():
    """Two rows with same position"""
    # As soon as we can actually work with multiple such rows, this test
    # needs to be updated since it currently just checks whether the second of
    # the positions is skipped.
    table = list(VcfReader("tests/data/duplicate-positions.vcf", only_snvs=False))[0]
    assert len(table.variants) == 2
    assert table.variants[0].position == 1
    assert table.variants[0].reference_allele == "A"
    assert table.variants[0].alternative_allele == "T"
    assert table.variants[1].position == 19
    assert table.variants[1].reference_allele == "G"
    assert table.variants[1].alternative_allele == "A"


def test_do_not_phase_duplicate_position(algorithm, tmpdir):
    """Ensure HP tag is added only to first of duplicate positions"""
    tmpvcf = str(tmpdir.join("duplicate-positions-phased.vcf"))
    run_whatshap(
        phase_input_files=["tests/data/oneread.bam"],
        variant_file="tests/data/duplicate-positions.vcf",
        output=tmpvcf,
        algorithm=algorithm,
    )
    seen_positions = set()
    records = list(pysam.VariantFile(tmpvcf))
    assert len(records) == 4
    for record in records:
        assert not (record.start in seen_positions and "HP" in record.format)
        seen_positions.add(record.start)


def test_multi_alt():
    """Skip multi-ALT in VCF"""
    table = list(VcfReader("tests/data/unknown-genotype.vcf"))[0]
    assert [variant.position for variant in table.variants] == [1, 4]


def assert_genotype_likelihoods(actual, expected):
    if expected is None:
        assert actual is None
        return
    for i in range(2):
        e = expected.log10_prob_of(i)
        a = actual.log10_prob_of(i)
        if e is None or a is None:
            assert a is None and e is None
        else:
            assert e == approx(a, rel=1e-6)


def test_read_genotype_likelihoods():
    tables = list(VcfReader("tests/data/genotype-likelihoods.vcf", genotype_likelihoods=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "chrA"
    assert table.samples == ["sample1", "sample2"]
    assert len(table.variants) == 4

    assert len(table.genotypes) == 2
    assert list(table.genotypes[0]) == canonic_index_list_to_biallelic_gt_list([2, 1, 1, 1])
    assert list(table.genotypes[1]) == canonic_index_list_to_biallelic_gt_list([1, 0, 0, 1])

    gl0 = GenotypeLikelihoods([-2.1206, -0.8195, -0.07525])
    gl1 = GenotypeLikelihoods([-10.3849, 0, -5.99143])
    gl2 = GenotypeLikelihoods([-2.1, None, -0.8])
    gl3 = GenotypeLikelihoods([0, -10.0, -0.6])

    assert len(table.genotype_likelihoods_of("sample1")) == 4
    assert len(table.genotype_likelihoods_of("sample2")) == 4

    expected1 = [gl0, gl2, None, gl0]
    expected2 = [gl1, gl3, None, gl1]
    for actual_gl, expected_gl in zip(table.genotype_likelihoods_of("sample1"), expected1):
        assert_genotype_likelihoods(actual_gl, expected_gl)
    for actual_gl, expected_gl in zip(table.genotype_likelihoods_of("sample2"), expected2):
        assert_genotype_likelihoods(actual_gl, expected_gl)


def test_genotype_likelihoods():
    assert list(PhredGenotypeLikelihoods([0, 0, 0])) == [0, 0, 0]
    assert list(PhredGenotypeLikelihoods([7, 1, 12])) == [7, 1, 12]
    gl = GenotypeLikelihoods([math.log10(x) for x in [1e-10, 0.5, 0.002]])
    assert list(gl.as_phred()) == [97, 0, 24]
    assert list(gl.as_phred(regularizer=0.01)) == [20, 0, 19]


def test_read_region():
    vcf_reader = VcfReader("tests/data/haplotag_1.vcf.gz")
    tableA = vcf_reader.fetch("chr1")
    tableB = vcf_reader.fetch("chr1", 1_069_570, 1_080_000)
    assert tableA.chromosome == tableB.chromosome
    assert len(tableA.variants) == len(tableB.variants)


def test_read_region_subsets():
    regions = [(1069570, 1070690), (1074910, 1076152)]
    vcf_reader = VcfReader("tests/data/haplotag_1.vcf.gz", only_snvs=False)
    table = vcf_reader.fetch_regions("chr1", regions)
    assert table.chromosome == "chr1"
    assert len(table.variants) == 8
    assert table.variants[5].reference_allele == "CG"
    assert table.variants[5].alternative_allele == "C"


def test_read_tetraploid_unphased():
    tables = list(VcfReader("tests/data/polyploid.chr22.unphased.vcf", phases=False))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "chr22"
    assert table.samples == ["HG00514_NA19240"]
    assert len(table.variants) == 8
    assert table.variants[0].reference_allele == "A"
    assert table.variants[0].alternative_allele == "C"
    assert table.variants[1].reference_allele == "G"
    assert table.variants[1].alternative_allele == "A"
    assert table.variants[2].reference_allele == "G"
    assert table.variants[2].alternative_allele == "T"
    assert table.variants[3].reference_allele == "G"
    assert table.variants[3].alternative_allele == "C"
    print("Got:")
    for genotype in table.genotypes[0]:
        print(genotype)
    print("Exp:")
    for genotype in canonic_index_list_to_biallelic_gt_list([3, 2, 0, 3, 3, 1, 1, 1]):
        print(genotype)
    assert table.genotypes[0] == canonic_index_list_to_biallelic_gt_list(
        [3, 2, 0, 3, 3, 1, 1, 1], 4
    )


def test_read_tetraploid_phased():
    tables = list(VcfReader("tests/data/polyploid.chr22.phased.vcf", phases=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "chr22"
    assert table.samples == ["HG00514_NA19240"]
    assert len(table.variants) == 8

    expected_phase = [
        VariantCallPhase(block_id=20000000, phase=(1, 0, 1, 1), quality=None),
        VariantCallPhase(block_id=20000000, phase=(1, 0, 1, 0), quality=None),
        None,
        VariantCallPhase(block_id=20000000, phase=(1, 0, 1, 1), quality=None),
        VariantCallPhase(block_id=20001000, phase=(1, 0, 1, 1), quality=None),
        VariantCallPhase(block_id=20001000, phase=(0, 0, 0, 1), quality=None),
        VariantCallPhase(block_id=20001000, phase=(0, 0, 0, 1), quality=None),
        VariantCallPhase(block_id=20001000, phase=(0, 0, 0, 1), quality=None),
    ]
    print("Got:")
    for variant in table.phases[0]:
        print(variant)
    print("Exp:")
    for variant in expected_phase:
        print(variant)
    assert list(table.phases[0]) == expected_phase


def test_read_tetraploid_genotype_likelihoods():
    tables = list(
        VcfReader(
            "tests/data/polyploid.chr22.unphased.vcf", phases=False, genotype_likelihoods=True
        )
    )
    assert len(tables) == 1
    table = tables[0]
    assert len(table.variants) == 8
    exp_gl = [
        GenotypeLikelihoods([-x / 10 for x in [19, 28, 29, 2, 10, 6]]),
        GenotypeLikelihoods([-x / 10 for x in [1, 8, 29, 24, 15, 23]]),
        GenotypeLikelihoods([-x / 10 for x in [25, 33, 35, 31, 0, 30]]),
        GenotypeLikelihoods([-x / 10 for x in [6, 27, 6, 3, 46, 42]]),
    ] * 2
    assert table.genotype_likelihoods_of(table.samples[0]) == exp_gl


def test_unsupported_ploidy():
    try:
        _ = list(VcfReader("tests/data/hexadecaploid.chr22.vcf", phases=False))
    except PloidyError:
        return
    assert False


def test_unsupported_ploidy_phased():
    try:
        _ = list(VcfReader("tests/data/hexadecaploid.chr22.vcf", phases=True))
    except PloidyError:
        return
    assert False


def test_inconsistent_ploidy():
    try:
        _ = list(VcfReader("tests/data/polyploid.chr22.inconsistent.vcf", phases=False))
    except PloidyError:
        return
    assert False


def test_inconsistent_ploidy_phased():
    try:
        _ = list(VcfReader("tests/data/polyploid.chr22.inconsistent.vcf", phases=True))
    except PloidyError:
        return
    assert False


def test_vcf_without_index(tmp_path):
    vcf_path = tmp_path / "file.vcf.gz"
    import shutil

    shutil.copy("tests/data/haplotag_1.vcf.gz", vcf_path)
    with raises(VcfIndexMissing):
        with VcfReader(vcf_path) as vr:
            list(vr.fetch("chr1"))
