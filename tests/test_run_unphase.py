from whatshap.unphase import run_unphase


def test_unphase(tmpdir):
	out = tmpdir.join("out.vcf")
	run_unphase('tests/data/phased-via-mixed-HP-PS.vcf', str(out))
	with open('tests/data/unphased.vcf') as f:
		expected = f.read()
	assert expected == out.read_text(encoding='ascii')
