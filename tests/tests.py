import tempfile
import os
import shutil
import pysam
from whatshap.__main__ import run_whatshap
from whatshap.vcf import VcfReader


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


def assert_phasing(phases, expected_phases):
	print('assert_phasing({}, {})'.format(phases, expected_phases))
	assert len(phases) == len(expected_phases)
	p_unchanged = []
	p_inverted = []
	p_expected = []
	for phase, expected_phase in zip(phases, expected_phases):
		if (phase is None) and (expected_phase is None):
			continue
		block_id, phase_bit = phase
		block_id_expected, phase_bit_expected = expected_phase
		assert block_id == block_id_expected
		p_unchanged.append(phase_bit)
		p_inverted.append(1-phase_bit)
		p_expected.append(phase_bit_expected)
	assert (p_unchanged == p_expected) or (p_inverted == p_expected)


def test_phase_three_individuals():
	tempdir = tempfile.mkdtemp()
	try:
		bamfile = tempdir + '/trio.pacbio.bam'
		outvcf = tempdir + '/output.vcf'
		pysam.view('tests/data/trio.pacbio.sam', '-Sb', '-o', bamfile, catch_stdout=False)
		pysam.index(bamfile, catch_stdout=False)
		run_whatshap(bam=[bamfile], vcf='tests/data/trio.vcf', output=outvcf)
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		assert_phasing(table.phases_of('HG004'), [None, (60907394,0), (60907394,0), (60907394,0), None])
		assert_phasing(table.phases_of('HG003'), [(60906167,0), None, (60906167,0), None, None])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None])

	finally:
		shutil.rmtree(tempdir)


def test_phase_trio():
	tempdir = tempfile.mkdtemp()
	try:
		bamfile = tempdir + '/trio.pacbio.bam'
		outvcf = tempdir + '/output.vcf'
		pysam.view('tests/data/trio.pacbio.sam', '-Sb', '-o', bamfile, catch_stdout=False)
		pysam.index(bamfile, catch_stdout=False)
		run_whatshap(bam=[bamfile], vcf='tests/data/trio.vcf', output=outvcf, ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		assert_phasing(table.phases_of('HG004'), [(60906167,0), (60906167,0), (60906167,0), (60906167,0), (60906167,0)])
		assert_phasing(table.phases_of('HG003'), [(60906167,0), None, (60906167,0), (60906167,0), (60906167,0)])
		assert_phasing(table.phases_of('HG002'), [None, (60906167,0), None, None, None])

	finally:
		shutil.rmtree(tempdir)
