from nose.tools import raises
from whatshap.vcf import parse_vcf

def test_read_phased():
	l = list(parse_vcf('tests/data/phasedinput.vcf'))
	assert len(l) == 1
	sample, chromosome, variants = l[0]
	assert sample == 'sample'
	assert chromosome == 'ref'
	assert len(variants) == 2
	assert variants[0].reference_allele == 'A'
	assert variants[0].alternative_allele == 'C'
	assert variants[1].reference_allele == 'G'
	assert variants[1].alternative_allele == 'T'

