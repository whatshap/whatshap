from whatshap.bam import SampleBamReader

def test_read():
	sbr = SampleBamReader('tests/data/oneread.bam')
	reads = list(sbr.fetch('ref', 'sample'))
	assert len(reads) == 1
	read = reads[0]
	assert read.opt('RG') == '1'
