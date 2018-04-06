from nose.tools import raises, assert_almost_equals
import os
import math
from tempfile import TemporaryDirectory
from whatshap.vcf import VcfReader, MixedPhasingError, VariantCallPhase, VcfVariant, GenotypeLikelihoods
from whatshap.phase import run_whatshap
from whatshap.core import PhredGenotypeLikelihoods


def test_read_phased():
	tables = list(VcfReader('tests/data/phasedinput.vcf', phases=True))
	assert len(tables) == 1
	table = tables[0]
	assert table.chromosome == 'ref'
	assert table.samples == ['sample']
	assert len(table.variants) == 2
	assert table.variants[0].reference_allele == 'A'
	assert table.variants[0].alternative_allele == 'C'
	assert table.variants[1].reference_allele == 'G'
	assert table.variants[1].alternative_allele == 'T'
	assert table.genotypes[0][0] == table.genotypes[0][1] == 1


def test_read_multisample_vcf():
	tables = list(VcfReader('tests/data/multisample.vcf'))
	assert len(tables) == 2
	table, table_b = tables
	assert table_b.chromosome == 'chrB'
	assert table_b.samples == ['sample1', 'sample2']

	assert table.chromosome == 'chrA'
	assert len(table.variants) == 3
	assert table.samples == ['sample1', 'sample2']

	assert table.variants[0].reference_allele == 'A'
	assert table.variants[0].alternative_allele == 'T'
	assert table.variants[1].reference_allele == 'C'
	assert table.variants[1].alternative_allele == 'G'
	assert table.variants[2].reference_allele == 'G'
	assert table.variants[2].alternative_allele == 'T'

	assert len(table.genotypes) == 2
	assert list(table.genotypes[0]) == [1, 1, 1]
	assert list(table.genotypes[1]) == [1, 1, 0]

	assert list(table.genotypes_of('sample1')) == [1, 1, 1]
	assert list(table.genotypes_of('sample2')) == [1, 1, 0]


def test_read_phased_vcf():
	for filename in ['tests/data/phased-via-HP.vcf', 'tests/data/phased-via-PS.vcf']:
		print('Testing', filename)
		tables = list(VcfReader(filename, phases=True))
		assert len(tables) == 2
		table_a, table_b = tables

		assert table_a.chromosome == 'chrA'
		assert len(table_a.variants) == 4
		assert table_a.samples == ['sample1', 'sample2']

		assert table_b.chromosome == 'chrB'
		assert len(table_b.variants) == 2
		assert table_b.samples == ['sample1', 'sample2']

		assert len(table_a.genotypes) == 2
		assert list(table_a.genotypes[0]) == [1, 2, 1, 1]
		assert list(table_a.genotypes[1]) == [1, 1, 1, 1]
		assert list(table_a.genotypes_of('sample1')) == [1, 2, 1, 1]
		assert list(table_a.genotypes_of('sample2')) == [1, 1, 1, 1]

		assert len(table_b.genotypes) == 2
		assert list(table_b.genotypes[0]) == [0, 1]
		assert list(table_b.genotypes[1]) == [1, 2]
		assert list(table_b.genotypes_of('sample1')) == [0, 1]
		assert list(table_b.genotypes_of('sample2')) == [1, 2]

		print(table_a.phases)
		assert len(table_a.phases) == 2
		expected_phase_sample1 = [
			None, 
			None, 
			VariantCallPhase(block_id=300, phase=1, quality=23),
			VariantCallPhase(block_id=300, phase=0, quality=42)
		]
		expected_phase_sample2 = [
			VariantCallPhase(block_id=100, phase=0, quality=10),
			VariantCallPhase(block_id=100, phase=1, quality=20),
			VariantCallPhase(block_id=300, phase=0, quality=30),
			VariantCallPhase(block_id=300, phase=0, quality=None)
		]
		assert list(table_a.phases[0]) == expected_phase_sample1
		assert list(table_a.phases[1]) == expected_phase_sample2
		assert list(table_a.phases_of('sample1')) == expected_phase_sample1
		assert list(table_a.phases_of('sample2')) == expected_phase_sample2

		assert len(table_b.phases) == 2
		assert list(table_b.phases[0]) == [None, None]
		assert list(table_b.phases[1]) == [None, None]
		assert list(table_b.phases_of('sample1')) == [None, None]
		assert list(table_b.phases_of('sample2')) == [None, None]


@raises(MixedPhasingError)
def test_mixed_phasing_vcf():
	tables = list(VcfReader('tests/data/phased-via-mixed-HP-PS.vcf', phases=True))


def test_vcf_variant_hashability():
	v = [
		VcfVariant(10, 'A', 'TC'),
		VcfVariant(10, 'A', 'TCA'),
		VcfVariant(10, 'C', 'TC'),
		VcfVariant(20, 'A', 'TC'),
		VcfVariant(10, 'A', 'TCA'),
		VcfVariant(20, 'A', 'TC')
	]
	assert len(set(v)) == 4


def test_phasing_to_reads():
	for filename in ['tests/data/phased-via-HP.vcf', 'tests/data/phased-via-PS.vcf']:
		tables = list(VcfReader(filename, phases=True))
		assert len(tables) == 2
		table_a, table_b = tables
		phase_reads_sample1 = list(table_a.phased_blocks_as_reads('sample1', table_a.variants, 17, 18, default_quality=90, mapq=101))
		print(phase_reads_sample1)
		assert len(phase_reads_sample1) == 1
		read = phase_reads_sample1[0]
		assert len(read) == 2
		assert read.name == 'sample1_block_300'
		assert read.source_id == 17
		assert read.mapqs == (101,)
		assert read[0].position == 300 - 1
		assert read[0].allele == 1
		assert read[0].quality == 23
		assert read[1].position == 350 - 1
		assert read[1].allele == 0
		assert read[1].quality == 42

		phase_reads_sample2 = list(table_a.phased_blocks_as_reads('sample2', table_a.variants, 11, 12, default_quality=91, mapq=102))
		print(phase_reads_sample2)
		assert len(phase_reads_sample2) == 2
		read1, read2 = phase_reads_sample2
		assert len(read1) == len(read2) == 2
		if read1[0].position > read2[0].position:
			read1, read2 = read2, read1
		assert read1.name == 'sample2_block_100'
		assert read1.source_id == 11
		assert read1.mapqs == (102,)
		assert read1[0].position == 100 - 1
		assert read1[0].allele == 0
		assert read1[0].quality == 10
		assert read1[1].position == 150 - 1
		assert read1[1].allele == 1
		assert read1[1].quality == 20
		assert read2.name == 'sample2_block_300'
		assert read2.source_id == 11
		assert read2.mapqs == (102,)
		assert read2[0].position == 300 - 1
		assert read2[0].allele == 0
		assert read2[0].quality == 30
		assert read2[1].position == 350 - 1
		assert read2[1].allele == 0
		assert read2[1].quality == 91

		variants = [
			VcfVariant(350 - 1, 'G', 'T'),
			VcfVariant(300 - 1, 'G', 'T'),
			VcfVariant(17, 'A', 'TTC'),
			VcfVariant(1000, 'C', 'G')
		]
		phase_reads_sample2 = list(table_a.phased_blocks_as_reads('sample2', variants, 11, 12, default_quality=91, mapq=102))
		print(phase_reads_sample2)
		assert len(phase_reads_sample2) == 1
		read = phase_reads_sample2[0]
		assert len(read) == 2
		assert read.name == 'sample2_block_300'
		assert read.source_id == 11
		assert read.mapqs == (102,)
		assert read[0].position == 300 - 1
		assert read[0].allele == 0
		assert read[0].quality == 30
		assert read[1].position == 350 - 1
		assert read[1].allele == 0
		assert read[1].quality == 91


def test_unknown_genotype():
	"""VCF with './.' genotype"""
	tables = list(VcfReader('tests/data/unknown-genotype.vcf'))
	assert tables[0].genotypes[1][0] == -1


def test_normalize():
	assert VcfVariant(100, 'A', 'C').normalized() == VcfVariant(100, 'A', 'C')
	assert VcfVariant(100, '', 'A').normalized() == VcfVariant(100, '', 'A')
	assert VcfVariant(100, 'A', '').normalized() == VcfVariant(100, 'A', '')
	assert VcfVariant(100, 'A', 'AC').normalized() == VcfVariant(101, '', 'C')
	assert VcfVariant(100, 'AC', 'A').normalized() == VcfVariant(101, 'C', '')
	assert VcfVariant(100, 'ACAGACC', 'ACAGACT').normalized() == VcfVariant(106, 'C', 'T')
	assert VcfVariant(100, 'GCTG', 'GCTAAA').normalized() == VcfVariant(103, 'G', 'AAA')
	assert VcfVariant(100, 'ATTA', 'ATA').normalized() == VcfVariant(101, 'T', '')
	assert VcfVariant(100, 'ATTTC', 'ATTTTTTC').normalized() == VcfVariant(101, '', 'TTT')
	assert VcfVariant(100, 'GCTGTT', 'GCTAAATT').normalized() == VcfVariant(103, 'G', 'AAA')


def test_read_duplicate_position():
	"""Two rows with same position"""
	# As soon as we can actually work with multiple such rows, this test
	# needs to be updated since it currently just checks whether the second of
	# the positions is skipped.
	table = list(VcfReader('tests/data/duplicate-positions.vcf', indels=True))[0]
	assert len(table.variants) == 2
	assert table.variants[0].position == 1
	assert table.variants[0].reference_allele == 'A'
	assert table.variants[0].alternative_allele == 'T'
	assert table.variants[1].position == 19
	assert table.variants[1].reference_allele == 'G'
	assert table.variants[1].alternative_allele == 'A'


def test_do_not_phase_duplicate_position():
	"""Ensure HP tag is added only to first of duplicate positions"""
	with TemporaryDirectory() as tmpdir:
		tmpvcf = os.path.join(tmpdir, 'duplicate-positions-phased.vcf')
		run_whatshap(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/duplicate-positions.vcf',
			output=tmpvcf)
		import vcf
		seen_positions = set()
		records = list(vcf.Reader(filename=tmpvcf))
		assert len(records) == 4
		for record in records:
			assert not (record.start in seen_positions and hasattr(record.samples[0].data, 'HP'))
			seen_positions.add(record.start)


def test_multi_alt():
	"""Skip multi-ALT in VCF"""
	table = list(VcfReader('tests/data/unknown-genotype.vcf'))[0]
	assert [ variant.position for variant in table.variants ] == [1, 4]


def assert_genotype_likelihoods(actual, expected):
	if expected is None:
		assert actual is None
		return
	for i in range(2):
		assert_almost_equals(actual.log10_prob_of(i), expected.log10_prob_of(i), places=10)


def test_read_genotype_likelihoods():
	tables = list(VcfReader('tests/data/genotype-likelihoods.vcf', genotype_likelihoods=True))
	assert len(tables) == 1
	table = tables[0]
	assert table.chromosome == 'chrA'
	assert table.samples == ['sample1', 'sample2']
	assert len(table.variants) == 4

	assert len(table.genotypes) == 2
	assert list(table.genotypes[0]) == [2, 1, 1, 1]
	assert list(table.genotypes[1]) == [1, 0, 0, 1]

	gl0 = GenotypeLikelihoods(-2.1206, -0.8195, -0.07525)
	gl1 = GenotypeLikelihoods(-10.3849, 0, -5.99143)
	gl2 = GenotypeLikelihoods(-2.1, -0.8, -0.8)
	gl3 = GenotypeLikelihoods(0, -10.0, -0.6)

	assert len(table.genotype_likelihoods_of('sample1')) == 4
	assert len(table.genotype_likelihoods_of('sample2')) == 4

	expected1 = [gl0, gl2, None, gl0]
	expected2 = [gl1, gl3, None, gl1]
	for actual_gl, expected_gl in zip(table.genotype_likelihoods_of('sample1'), expected1):
		assert_genotype_likelihoods(actual_gl, expected_gl)
	for actual_gl, expected_gl in zip(table.genotype_likelihoods_of('sample2'), expected2):
		assert_genotype_likelihoods(actual_gl, expected_gl)


def test_genotype_likelihoods():
	assert list(PhredGenotypeLikelihoods()) == [0, 0, 0]
	assert list(PhredGenotypeLikelihoods(7, 1, 12)) == [7, 1, 12]
	gl = GenotypeLikelihoods( *(math.log10(x) for x in [1e-10, 0.5, 0.002]) )
	assert list(gl.as_phred()) == [97, 0, 24]
	assert list(gl.as_phred(regularizer=0.01)) == [20, 0, 19]
