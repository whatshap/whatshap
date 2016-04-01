from nose.tools import raises
from whatshap.vcf import VcfReader

def test_read_phased():
	l = list(VcfReader('tests/data/phasedinput.vcf'))
	assert len(l) == 1
	chromosome, samples = l[0]
	assert list(samples.keys()) == ['sample']
	assert chromosome == 'ref'
	calls = samples['sample']
	assert len(calls) == 2
	assert calls[0].reference_allele == 'A'
	assert calls[0].alternative_allele == 'C'
	assert calls[1].reference_allele == 'G'
	assert calls[1].alternative_allele == 'T'
	assert calls[0].genotype == calls[1].genotype == 1


def test_read_multisample_vcf():
	l = list(VcfReader('tests/data/multisample.vcf'))
	assert len(l) == 2
	assert l[1][0] == 'chrB'
	chromosome, samples = l[0]
	assert chromosome == 'chrA'
	assert len(samples) == 2
	assert set(samples.keys()) == set(['sample1', 'sample2'])
	calls1 = samples['sample1']
	calls2 = samples['sample2']
	assert len(calls1) == 3
	assert len(calls2) == 3
	assert calls1[0].reference_allele == 'A'
	assert calls1[0].alternative_allele == 'T'
	assert calls1[1].reference_allele == 'C'
	assert calls1[1].alternative_allele == 'G'
	assert calls1[2].reference_allele == 'G'
	assert calls1[2].alternative_allele == 'T'
