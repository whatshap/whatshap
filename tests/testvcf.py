from nose.tools import raises
from whatshap.vcf import VcfReader, MixedPhasingError


def test_read_phased():
	tables = list(VcfReader('tests/data/phasedinput.vcf'))
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
		tables = list(VcfReader(filename))
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
		assert list(table_a.phases[0]) == [None, None, (300,1), (300,0)]
		assert list(table_a.phases[1]) == [(100,0), (100,1), (300,0), (300,0)]
		assert list(table_a.phases_of('sample1')) == [None, None, (300,1), (300,0)]
		assert list(table_a.phases_of('sample2')) == [(100,0), (100,1), (300,0), (300,0)]

		assert len(table_b.phases) == 2
		assert list(table_b.phases[0]) == [None, None]
		assert list(table_b.phases[1]) == [None, None]
		assert list(table_b.phases_of('sample1')) == [None, None]
		assert list(table_b.phases_of('sample2')) == [None, None]


@raises(MixedPhasingError)
def test_mixed_phasing_vcf():
	tables = list(VcfReader('tests/data/phased-via-mixed-HP-PS.vcf'))