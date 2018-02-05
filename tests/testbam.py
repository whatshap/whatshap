from nose.tools import raises
from whatshap.bam import SampleBamReader, SampleNotFoundError, AlignmentFileNotIndexedError


def test_read():
	sbr = SampleBamReader('tests/data/oneread.bam')
	reads = list(sbr.fetch('ref', 'sample'))
	assert len(reads) == 1
	read = reads[0]
	assert read.bam_alignment.opt('RG') == '1'


@raises(SampleNotFoundError)
def test_read_sample_not_found():
	sbr = SampleBamReader('tests/data/oneread.bam')
	reads = list(sbr.fetch('ref', 'non-existing-sample'))


def test_read_cram():
	sbr = SampleBamReader('tests/data/oneread.cram', reference='tests/data/oneread-ref.fasta')
	reads = list(sbr.fetch('ref', 'sample'))
	assert len(reads) == 1
	assert reads[0].bam_alignment.opt('RG') == '1'


@raises(AlignmentFileNotIndexedError)
def test_no_index():
	SampleBamReader('tests/data/not-indexed.bam')
