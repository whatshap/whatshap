import os.path
from nose.tools import raises
from whatshap.utils import UnknownFileFormatError, detect_file_format, IndexedFasta, FastaNotIndexedError


def test_detect_alignment_file_format():
	assert detect_file_format('tests/data/oneread.bam') == 'BAM'
	assert detect_file_format('tests/data/oneread.cram') == 'CRAM'
	assert detect_file_format('tests/data/onevariant.vcf') == 'VCF'
	assert detect_file_format('tests/data/onevariant.vcf.gz') == 'VCF'


@raises(UnknownFileFormatError)
def test_detect_alignment_file_format():
	detect_file_format('tests/data/pedigree.ped')


@raises(FastaNotIndexedError)
def test_missing_fai():
	assert not os.path.exists('tests/data/not-indexed.fasta.fai')
	IndexedFasta('tests/data/not-indexed.fasta')
