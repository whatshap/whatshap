from whatshap.unphase import remove_phasing
from io import StringIO


def test_unphase():
	out = StringIO()
	remove_phasing('tests/data/phased-via-mixed-HP-PS.vcf', out)
	with open('tests/data/unphased.vcf') as f:
		expected = f.read()
	actual = out.getvalue()
	assert expected == actual
