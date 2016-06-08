from tempfile import TemporaryDirectory
import os
import shutil
import pysam
from nose.tools import raises

from whatshap.phase import run_whatshap
from whatshap.vcf import VcfReader, VariantCallPhase

trio_bamfile = 'tests/data/trio.pacbio.bam'
trio_merged_bamfile = 'tests/data/trio-merged-blocks.bam'
trio_paired_end_bamfile = 'tests/data/paired_end.sorted.bam'


def setup_module():
	# This function is run once for this module
	pysam.view('tests/data/trio.pacbio.sam', '-b', '-o', trio_bamfile, catch_stdout=False)
	pysam.index(trio_bamfile, catch_stdout=False)
	pysam.view('tests/data/trio-merged-blocks.sam', '-b', '-o', trio_merged_bamfile, catch_stdout=False)
	pysam.index(trio_merged_bamfile, catch_stdout=False)
	pysam.view('tests/data/paired_end.sorted.sam', '-b', '-o', trio_paired_end_bamfile, catch_stdout=False)
	pysam.index(trio_paired_end_bamfile, catch_stdout=False)


def teardown_module():
	os.remove(trio_bamfile)
	os.remove(trio_bamfile + '.bai')
	os.remove(trio_merged_bamfile)
	os.remove(trio_merged_bamfile + '.bai')
	os.remove(trio_paired_end_bamfile)
	os.remove(trio_paired_end_bamfile + '.bai')


def test_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	assert LooseVersion(pysam_version) >= LooseVersion("0.8.1")


def test_one_variant():
	run_whatshap(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null')


def test_default_output():
	"""Output to stdout"""
	run_whatshap(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf')


def test_bam_without_readgroup():
	run_whatshap(phase_input_files=['tests/data/no-readgroup.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', ignore_read_groups=True)


@raises(SystemExit)
def test_requested_sample_not_found():
	run_whatshap(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', samples=['DOES_NOT_EXIST'])


def test_with_reference():
	run_whatshap(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
		reference='tests/data/pacbio/reference.fasta')


def assert_phasing(phases, expected_phases):
	print('assert_phasing({}, {})'.format(phases, expected_phases))
	assert len(phases) == len(expected_phases)
	p_unchanged = []
	p_inverted = []
	p_expected = []
	for phase, expected_phase in zip(phases, expected_phases):
		if (phase is None) and (expected_phase is None):
			continue
		assert phase is not None and expected_phase is not None
		assert phase.block_id == expected_phase.block_id
		p_unchanged.append(phase.phase)
		p_inverted.append(1-phase.phase)
		p_expected.append(expected_phase.phase)
	assert (p_unchanged == p_expected) or (p_inverted == p_expected)


def test_phase_three_individuals():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf)
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		phase1 = VariantCallPhase(60906167, 0, None)
		phase3 = VariantCallPhase(60907394, 0, None)
		assert_phasing(table.phases_of('HG004'), [None, phase3, phase3, phase3, None])
		assert_phasing(table.phases_of('HG003'), [phase1, None, phase1, None, None])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None])


def test_phase_one_of_three_individuals():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf, samples=['HG003'])
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		phase0 = VariantCallPhase(60906167,0,None)
		assert_phasing(table.phases_of('HG004'), [None, None, None, None, None])
		assert_phasing(table.phases_of('HG003'), [phase0, None, phase0, None, None])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None])


def test_phase_trio():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		phase0 = VariantCallPhase(60906167, 0, None)
		assert_phasing(table.phases_of('HG004'), [phase0, phase0, phase0, phase0, phase0])
		assert_phasing(table.phases_of('HG003'), [phase0, None, phase0, phase0, phase0])
		assert_phasing(table.phases_of('HG002'), [None, phase0, None, None, None])


def test_phase_trio_merged_blocks():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output-merged-blocks.vcf'
		run_whatshap(phase_input_files=[trio_merged_bamfile], variant_file='tests/data/trio-merged-blocks.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 8
		assert table.samples == ['HG002', 'HG003', 'HG004']
		assert table.num_of_blocks_of('HG004') == 1
		assert table.num_of_blocks_of('HG003') == 1
		assert table.num_of_blocks_of('HG002') == 1

		phase0 = VariantCallPhase(752566, 0, None)
		phase1 = VariantCallPhase(752566, 1, None)
		assert_phasing(table.phases_of('HG004'), [phase1, phase1, phase1, None, phase1, phase1, phase1, phase1])
		assert_phasing(table.phases_of('HG003'), [None, None, None, None, phase0, phase0, phase0, phase1])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None, None, None, phase1])


def test_phase_trio_dont_merge_blocks():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output-merged-blocks.vcf'
		run_whatshap(phase_input_files=[trio_merged_bamfile], variant_file='tests/data/trio-merged-blocks.vcf', output=outvcf,
				ped='tests/data/trio.ped', genmap='tests/data/trio.map', genetic_haplotyping=False)
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 8
		assert table.samples == ['HG002', 'HG003', 'HG004']
		assert table.num_of_blocks_of('HG004') == 2
		assert table.num_of_blocks_of('HG003') == 1
		assert table.num_of_blocks_of('HG002') == 1

		phase1 = VariantCallPhase(752566, 1, None)
		phase2_0 = VariantCallPhase(853954, 0, None)
		phase2_1 = VariantCallPhase(853954, 1, None)
		assert_phasing(table.phases_of('HG004'), [phase1, phase1, phase1, None, phase2_1, phase2_1, phase2_1, phase2_1])
		assert_phasing(table.phases_of('HG003'), [None, None, None, None, phase2_0, phase2_0, phase2_0, phase2_1])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None, None, None, phase2_1])


def test_phase_mendelian_conflict():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio-mendelian-conflict.vcf', output=outvcf,
				ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		phase = VariantCallPhase(60906167, 0, None)
		assert_phasing(table.phases_of('HG004'), [phase, None, phase, phase, phase])
		assert_phasing(table.phases_of('HG003'), [phase, None, phase, phase, phase])
		assert_phasing(table.phases_of('HG002'), [None, None, None, None, None])


def test_phase_missing_genotypes():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio-missing-genotypes.vcf', output=outvcf,
				ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		phase = VariantCallPhase(60906167, 0, None)
		assert_phasing(table.phases_of('HG004'), [phase, phase, None, phase, None])
		assert_phasing(table.phases_of('HG003'), [phase, None, None, phase, None])
		assert_phasing(table.phases_of('HG002'), [None, phase, None, None, None])


def test_phase_specific_chromosome():
	for requested_chromosome in ['1','2']:
		with TemporaryDirectory() as tempdir:
			outvcf = tempdir + '/output.vcf'
			run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio-two-chromosomes.vcf', output=outvcf,
					ped='tests/data/trio.ped', genmap='tests/data/trio.map', chromosomes=[requested_chromosome])
			assert os.path.isfile(outvcf)

			tables = list(VcfReader(outvcf))
			assert len(tables) == 2
			for table in tables:
				assert len(table.variants) == 5
				assert table.samples == ['HG004', 'HG003', 'HG002']
				if table.chromosome == '1' == requested_chromosome:
					phase0 = VariantCallPhase(60906167, 0, None)
					assert_phasing(table.phases_of('HG004'), [phase0, phase0, phase0, phase0, phase0])
					assert_phasing(table.phases_of('HG003'), [phase0, None, phase0, phase0, phase0])
					assert_phasing(table.phases_of('HG002'), [None, phase0, None, None, None])
				else:
					assert_phasing(table.phases_of('HG004'), [None, None, None, None, None])
					assert_phasing(table.phases_of('HG003'), [None, None, None, None, None])
					assert_phasing(table.phases_of('HG002'), [None, None, None, None, None])


def test_phase_trio_paired_end_reads():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output-paired_end.vcf'
		run_whatshap(phase_input_files=[trio_paired_end_bamfile], variant_file='tests/data/paired_end.sorted.vcf', output=outvcf,
		        ped='tests/data/trio_paired_end.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 3
		assert table.samples == ['mother', 'father', 'child']
		assert table.num_of_blocks_of('mother') == 1
		assert table.num_of_blocks_of('father') == 0
		assert table.num_of_blocks_of('child') == 1

		phase0 = VariantCallPhase(80050, 0, None)
		phase1 = VariantCallPhase(80050, 1, None)

		assert_phasing(table.phases_of('mother'), [phase1, phase1, phase0])
		assert_phasing(table.phases_of('father'), [None, None, None])
		assert_phasing(table.phases_of('child'), [None, None, phase1])
