import os
import pysam
import math

import pytest
from pysam import VariantFile
from whatshap.cli.genotype import run_genotype
from whatshap.cli import CommandLineError
from whatshap.vcf import VcfReader

trio_bamfile = "tests/data/trio.pacbio.bam"
trio_merged_bamfile = "tests/data/trio-merged-blocks.bam"
trio_paired_end_bamfile = "tests/data/paired_end.sorted.bam"
ped_samples_bamfile = "tests/data/ped_samples.bam"
recombination_breaks_bamfile = "tests/data/recombination_breaks.sorted.bam"
quartet2_bamfile = "tests/data/quartet2.bam"
short_bamfile = "tests/data/short-genome/short.bam"
indels_bamfile = "tests/data/indels.bam"

bam_files = [
    trio_bamfile,
    trio_merged_bamfile,
    trio_paired_end_bamfile,
    recombination_breaks_bamfile,
    quartet2_bamfile,
    short_bamfile,
    indels_bamfile,
]


def setup_module():
    # This function is run once for this module
    for bam_path in bam_files:
        assert bam_path.endswith(".bam")
        sam_path = bam_path[:-4] + ".sam"
        pysam.view(sam_path, "-b", "-o", bam_path, catch_stdout=False)
        pysam.index(bam_path, catch_stdout=False)


def test_one_variant():
    run_genotype(
        phase_input_files=["tests/data/oneread.bam"],
        variant_file="tests/data/onevariant.vcf",
        output="/dev/null",
    )


def test_default_output():
    """Output to stdout"""
    run_genotype(
        phase_input_files=["tests/data/oneread.bam"], variant_file="tests/data/onevariant.vcf"
    )


def test_bam_without_readgroup():
    run_genotype(
        phase_input_files=["tests/data/no-readgroup.bam"],
        variant_file="tests/data/onevariant.vcf",
        output="/dev/null",
        ignore_read_groups=True,
    )


def test_requested_sample_not_found():
    with pytest.raises(CommandLineError):
        run_genotype(
            phase_input_files=["tests/data/oneread.bam"],
            variant_file="tests/data/onevariant.vcf",
            output="/dev/null",
            samples=["DOES_NOT_EXIST"],
        )


def test_with_reference():
    run_genotype(
        phase_input_files=["tests/data/pacbio/pacbio.bam"],
        variant_file="tests/data/pacbio/variants.vcf",
        reference="tests/data/pacbio/reference.fasta",
    )


@pytest.mark.parametrize("priors", [True, False])
def test_only_snvs(tmpdir, priors):
    prioroutput = str(tmpdir.join("priors.vcf")) if priors else None
    outvcf = str(tmpdir.join("output_gl.vcf"))
    run_genotype(
        phase_input_files=["tests/data/pacbio/pacbio.bam"],
        variant_file="tests/data/pacbio/variants.vcf",
        reference="tests/data/pacbio/reference.fasta",
        output=outvcf,
        only_snvs=True,
        nopriors=not priors,
        prioroutput=prioroutput,
    )
    result_vcfs = [outvcf]
    if priors:
        result_vcfs.append(prioroutput)

    # make sure indels not genotyped (also in priors.vcf if computed)
    for o_vcf in result_vcfs:
        vcf_reader = VariantFile(o_vcf)
        default_l = math.log10(1 / 3.0)

        for record in vcf_reader:
            if record.alts is None:
                for call in record.samples.values():
                    assert set(call) == {"GT"}
            elif len(record.alts[0]) != len(record.ref):
                for call in record.samples.values():
                    for v in call["GL"]:
                        assert pytest.approx(default_l) == v


def test_multiallelic(tmpdir):
    outvcf = str(tmpdir.join("output_multi.vcf"))
    run_genotype(
        phase_input_files=["tests/data/pacbio/pacbio.bam"],
        variant_file="tests/data/multiallelic.vcf",
        reference="tests/data/pacbio/reference.fasta",
        output=outvcf,
        only_snvs=True,
    )
    vcf_reader = VariantFile(outvcf)
    for record in vcf_reader:
        n_alleles = len(record.alts) + 1
        if n_alleles > 1:
            for call in record.samples.values():
                assert len(call["GL"]) == ((n_alleles + 1) * n_alleles) / 2


def likeliest_genotype(a, b, c, thres):
    prob_a = 10**a
    prob_b = 10**b
    prob_c = 10**c

    prob = sorted([(prob_a, 0), (prob_b, 1), (prob_c, 2)])

    if prob[2][0] > prob[1][0] and prob[2][0] > thres:
        return prob[2][1]
    else:
        return None


@pytest.mark.parametrize("threshold", [0, 2, 3, 6, 13, 50])
def test_gt_quality_threshold(threshold, tmpdir):
    thres = 1 - 10 ** (-threshold / 10.0)

    out_vcf = str(tmpdir.join("out.vcf"))
    priors_vcf = str(tmpdir.join("priors.vcf"))
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio.vcf",
        output=out_vcf,
        gt_qual_threshold=threshold,
        only_snvs=True,
        prioroutput=priors_vcf,
    )

    for path in [out_vcf, priors_vcf]:
        for record in VariantFile(path):
            for call in record.samples.values():
                likelihoods = call["GL"]
                genotype = call["GT"]
                if genotype == (None,):
                    genotype = None
                else:
                    genotype = genotype[0] + genotype[1]
                gt = likeliest_genotype(likelihoods[0], likelihoods[1], likelihoods[2], thres)
                # print(likelihoods[0], likelihoods[1], likelihoods[2], gt, genotype, thres)
                assert gt == genotype


def test_genotyping_one_of_three_individuals(tmp_path):
    outvcf = tmp_path / "output.vcf"
    outpriors = tmp_path / "priors.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio.vcf",
        output=outvcf,
        samples=["HG003"],
        prioroutput=outpriors,
    )

    for outfile in [outvcf, outpriors]:
        assert os.path.isfile(outfile)

        tables = list(VcfReader(outfile, phases=True, genotype_likelihoods=True))
        assert len(tables) == 1
        table = tables[0]
        assert table.chromosome == "1"
        assert len(table.variants) == 5
        assert table.samples == ["HG004", "HG003", "HG002"]

        # there should be no genotype predictions for HG003/HG002
        default_l = math.log10(1 / 3.0)
        for l in [
            table.genotype_likelihoods_of("HG002"),
            table.genotype_likelihoods_of("HG004"),
        ]:
            for var in l:
                for v in var.log10_probs():
                    assert pytest.approx(default_l) == v


def test_use_ped_samples(tmp_path):
    outvcf = tmp_path / "output_ped_samples.vcf"
    run_genotype(
        phase_input_files=[ped_samples_bamfile],
        variant_file="tests/data/ped_samples.vcf",
        output=outvcf,
        ped="tests/data/trio.ped",
        genmap="tests/data/trio.map",
        use_ped_samples=True,
    )
    assert os.path.isfile(outvcf)

    tables = list(VcfReader(outvcf, phases=True, genotype_likelihoods=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "1"
    assert len(table.variants) == 5
    assert table.samples == ["HG004", "HG003", "HG002", "orphan"]

    default_l = math.log10(1 / 3.0)
    for var in table.genotype_likelihoods_of("orphan"):
        for v in var.log10_probs():
            assert pytest.approx(default_l) == v


@pytest.mark.parametrize(
    "sample_set",
    [["HG002"], ["HG003"], ["HG004"], ["HG002", "HG003"], ["HG002", "HG004"], ["HG003", "HG004"]],
)
def test_ped_sample(sample_set, tmp_path):
    # running with --ped and --sample on subset of trio,
    # should give same results as running with only --sample
    # the trio information should be ignored
    outvcf1 = tmp_path / "output1.vcf"
    outvcf2 = tmp_path / "output2.vcf"
    run_genotype(
        phase_input_files=[ped_samples_bamfile],
        variant_file="tests/data/ped_samples.vcf",
        output=outvcf1,
        ped="tests/data/trio.ped",
        samples=sample_set,
    )
    run_genotype(
        phase_input_files=[ped_samples_bamfile],
        variant_file="tests/data/ped_samples.vcf",
        output=outvcf2,
        samples=sample_set,
    )
    assert os.path.isfile(outvcf1)
    assert os.path.isfile(outvcf2)
    tables1 = list(VcfReader(outvcf1, phases=True, genotype_likelihoods=True))
    tables2 = list(VcfReader(outvcf2, phases=True, genotype_likelihoods=True))
    assert (len(tables1) == 1) and (len(tables2) == 1)
    table1, table2 = tables1[0], tables2[0]

    for individual in sample_set:
        for var1, var2 in zip(
            table1.genotype_likelihoods_of(individual), table2.genotype_likelihoods_of(individual)
        ):
            print(var1, var2)
            assert var1.log10_probs() == var2.log10_probs()


def test_genotyping_trio(tmp_path):
    outvcf = tmp_path / "output.vcf"
    outpriors = tmp_path / "priors.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio.vcf",
        output=outvcf,
        ped="tests/data/trio.ped",
        genmap="tests/data/trio.map",
        prioroutput=outpriors,
    )

    for outfile in [outvcf, outpriors]:
        assert os.path.isfile(outfile)

        tables = list(VcfReader(outfile, phases=True))
        assert len(tables) == 1
        table = tables[0]
        assert table.chromosome == "1"
        assert len(table.variants) == 5
        assert table.samples == ["HG004", "HG003", "HG002"]


@pytest.mark.parametrize("chromosome", ["1", "2"])
def test_genotyping_specific_chromosome(chromosome, tmp_path):
    outvcf = tmp_path / "output.vcf"
    outpriors = tmp_path / "priors.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio-two-chromosomes.vcf",
        output=outvcf,
        ped="tests/data/trio.ped",
        genmap="tests/data/trio.map",
        chromosomes=[chromosome],
        prioroutput=outpriors,
    )
    for outfile in [outvcf, outpriors]:
        assert os.path.isfile(outfile)
        tables = list(VcfReader(outfile, genotype_likelihoods=True))
        assert len(tables) == 2
        for table in tables:
            assert len(table.variants) == 5
            assert table.samples == ["HG004", "HG003", "HG002"]

        index = 0
        if chromosome == "1":
            index = 1

        # should be no genotype likelihoods for skipped chromosomes
        for s in tables[index].samples:
            assert tables[index].genotype_likelihoods_of(s) == [None] * 5
            assert tables[not index].genotype_likelihoods_of(s) != [None] * 5


def test_genotype_likelihoods_given(tmp_path):
    outvcf = tmp_path / "output_gl.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio_genotype_likelihoods.vcf",
        output=outvcf,
        ped="tests/data/trio.ped",
        genmap="tests/data/trio.map",
    )
    assert os.path.isfile(outvcf)
    tables = list(VcfReader(outvcf, phases=True, genotype_likelihoods=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "1"
    assert len(table.variants) == 5
    assert table.samples == ["HG004", "HG003", "HG002"]

    # check if PL likelihoods (that were present before) are deleted
    vcf_reader = VariantFile(outvcf)
    # print(list(vcf_reader.samples), outvcf)
    for record in vcf_reader:
        for call in record.samples.values():
            PL = call.get("PL", None)
            GL = call.get("GL", None)
            print("GL:", GL, "PL:", PL)
            assert PL == (None, None, None)
            assert GL is not None


# GL field was already present, make sure it is replaced by new likelihoods
def test_genotype_log_likelihoods_given(tmp_path):
    outvcf = tmp_path / "output_gl_log.vcf"
    outpriors = tmp_path / "priors.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio_genotype_log_likelihoods.vcf",
        output=outvcf,
        ped="tests/data/trio.ped",
        genmap="tests/data/trio.map",
        gt_qual_threshold=0,
        prioroutput=outpriors,
    )
    for outfile in [outvcf, outpriors]:
        assert os.path.isfile(outfile)
        tables = list(VcfReader(outfile, phases=True, genotype_likelihoods=True))
        assert len(tables) == 1
        table = tables[0]
        assert table.chromosome == "1"
        assert len(table.variants) == 5
        assert table.samples == ["HG004", "HG003", "HG002"]

        # check if GL likelihoods were replaced
        vcf_reader = VariantFile(outfile)
        print(list(vcf_reader.header.samples), outfile)
        for record in vcf_reader:
            for call in record.samples.values():
                GL = call.get("GL", None)
                GQ = call.get("GQ", None)
                print("GL:", GL, "GQ", GQ)
                assert GL != [-1, -1, -1]
                assert GQ != 100


def test_empty_format_field(tmp_path):
    outvcf = tmp_path / "output_empty_format.vcf"
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/empty_format.vcf",
        output=outvcf,
        gt_qual_threshold=0,
    )

    # check if sample fields now contain information
    assert os.path.isfile(outvcf)
    vcf_reader = VariantFile(outvcf)
    for record in vcf_reader:
        for sample, call in record.samples.items():
            assert set(call) == {"GT", "GL", "GQ"}


def test_phase_trio_paired_end_reads(tmp_path):
    outvcf = tmp_path / "output-paired_end.vcf"
    run_genotype(
        phase_input_files=[trio_paired_end_bamfile],
        variant_file="tests/data/paired_end.sorted.vcf",
        output=outvcf,
        ped="tests/data/trio_paired_end.ped",
        genmap="tests/data/trio.map",
    )
    assert os.path.isfile(outvcf)
    tables = list(VcfReader(outvcf, phases=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "1"
    assert len(table.variants) == 3
    assert table.samples == ["mother", "father", "child"]


def test_wrong_chromosome(tmp_path):
    outvcf = tmp_path / "output.vcf"
    with pytest.raises(CommandLineError):
        run_genotype(
            phase_input_files=[short_bamfile],
            ignore_read_groups=True,
            variant_file="tests/data/short-genome/wrongchromosome.vcf",
            output=outvcf,
        )


def extract_likelihoods(record):
    return [10.0**gl for gl in record.samples[0]["GL"]]


@pytest.mark.parametrize("constant", [0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 5, 10, 20, 100])
def test_adding_constant(constant, tmpdir):
    priors_raw_vcf = str(tmpdir.join("output.raw_priors.vcf"))
    outvcf_raw_vcf = str(tmpdir.join("output_raw.vcf"))
    priors_const_vcf = str(tmpdir.join("output.const_priors.vcf"))
    outvcf_const_vcf = str(tmpdir.join("output_raw.vcf"))

    # run genotyping without adding constant to priors
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio.vcf",
        prioroutput=priors_raw_vcf,
        output=outvcf_raw_vcf,
        only_snvs=True,
    )

    # run genotyping with modified priors
    run_genotype(
        phase_input_files=[trio_bamfile],
        variant_file="tests/data/trio.vcf",
        prioroutput=priors_const_vcf,
        output=outvcf_const_vcf,
        only_snvs=True,
        constant=constant,
    )

    # check if priors were modified properly
    with pysam.VariantFile(priors_raw_vcf) as f:
        records_raw = list(f)
    with pysam.VariantFile(priors_const_vcf) as f:
        records_const = list(f)

    assert len(records_raw) == len(records_const)

    for record_raw, record_const in zip(records_raw, records_const):
        likelihoods_raw = extract_likelihoods(record_raw)
        likelihoods_const = extract_likelihoods(record_const)
        norm_sum = likelihoods_raw[0] + likelihoods_raw[1] + likelihoods_raw[2] + 3.0 * constant
        # print(float(likelihoods_const[0]), float((likelihoods_raw[0] + constant)/norm_sum))
        # print(float(likelihoods_const[1]), float((likelihoods_raw[1] + constant)/norm_sum))
        # print(float(likelihoods_const[2]), float((likelihoods_raw[2] + constant)/norm_sum))

        for j in range(3):
            assert (
                pytest.approx(likelihoods_const[j], 1e-5)
                == (likelihoods_raw[j] + constant) / norm_sum
            )
