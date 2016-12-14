import os
from tempfile import TemporaryDirectory
from subprocess import call
from whatshap.haplofasta import run_haplofasta


def files_equal(path1, path2):
	return call(['diff', '-u', path1, path2]) == 0


def test_haplofasta():
	with TemporaryDirectory() as tempdir:
		run_haplofasta(
			vcf_path='tests/data/haplofasta.vcf',
			reference_path='tests/data/haplofasta.fasta',
			destination_folder=tempdir)
		base = os.path.join(tempdir, 'sample')
		assert files_equal('tests/data/haplofasta.1.fasta', base + '.1.fasta')
		assert files_equal('tests/data/haplofasta.2.fasta', base + '.2.fasta')
