import tempfile
import os
import shutil
import pysam
from whatshap.scripts.whatshap import run_whatshap

def test_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	assert LooseVersion(pysam_version) >= LooseVersion("0.8.1")


def test_one_variant():
	run_whatshap(bam=['tests/data/oneread.bam'], vcf='tests/data/onevariant.vcf',
		output='/dev/null')


def test_bam_without_readgroup():
	run_whatshap(bam=['tests/data/no-readgroup.bam'], vcf='tests/data/onevariant.vcf',
		output='/dev/null', ignore_read_groups=True)


def test_phase_three_individuals():
	tempdir = tempfile.mkdtemp()
	try:
		bamfile = tempdir + '/trio.pacbio.bam'
		outvcf = tempdir + '/output.vcf'
		pysam.view('tests/data/trio.pacbio.sam', '-Sb', '-o', bamfile, catch_stdout=False)
		pysam.index(bamfile, catch_stdout=False)
		run_whatshap(bam=[bamfile], vcf='tests/data/trio.vcf', output=outvcf)
		assert os.path.isfile(outvcf)
	finally:
		shutil.rmtree(tempdir)
