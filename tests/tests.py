from tempfile import TemporaryDirectory
import os
from io import StringIO
import pysam
from nose.tools import raises
from collections import namedtuple

from whatshap.phase import run_whatshap
from whatshap.haplotag import run_haplotag
from whatshap.compare import run_compare
from whatshap.vcf import VcfReader, VariantCallPhase

trio_bamfile = 'tests/data/trio.pacbio.bam'
trio_merged_bamfile = 'tests/data/trio-merged-blocks.bam'
trio_paired_end_bamfile = 'tests/data/paired_end.sorted.bam'
recombination_breaks_bamfile = 'tests/data/recombination_breaks.sorted.bam'

bam_files = [trio_bamfile, trio_merged_bamfile, trio_paired_end_bamfile, recombination_breaks_bamfile]


def setup_module():
	# This function is run once for this module
	for bam_path in bam_files:
		assert bam_path.endswith('.bam')
		sam_path = bam_path[:-4] + '.sam'
		pysam.view(sam_path, '-b', '-o', bam_path, catch_stdout=False)
		pysam.index(bam_path, catch_stdout=False)


def teardown_module():
	for path in bam_files:
		os.remove(path)
		os.remove(path + '.bai')


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


def test_with_reference_and_indels():
	run_whatshap(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
		reference='tests/data/pacbio/reference.fasta', indels=True)


def test_ps_tag():
	out = StringIO()
	run_whatshap(variant_file='tests/data/trio.vcf', phase_input_files=['tests/data/trio.pacbio.bam'],
	    output=out, tag='PS')
	out.seek(0)
	lines = [ line for line in out.readlines() if not line.startswith('#') ]

	# TODO This is quite an ugly way to test phased VCF writing
	assert lines[0] == "1\t60906167\t.\tG\tA\t.\tPASS\tAC=2;AN=6\tGT:PS\t0/1:.\t0|1:60906167\t0/0:.\n"
	assert lines[1]	== "1\t60907394\t.\tG\tA\t.\tPASS\tAC=4;AN=6\tGT:PS\t0|1:60907394\t1/1:.\t0/1:.\n"
	assert lines[2] == "1\t60907460\t.\tG\tT\t.\tPASS\tAC=2;AN=6\tGT:PS\t0|1:60907394\t0|1:60906167\t0/0:.\n"
	assert lines[3] == "1\t60907473\t.\tC\tA\t.\tPASS\tAC=2;AN=6\tGT:PS\t0|1:60907394\t0/1:.\t0/0:.\n"
	assert lines[4] == "1\t60909718\t.\tT\tC\t.\tPASS\tAC=2;AN=6\tGT\t0/1\t0/1\t0/0\n"


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
		
def test_phase_quartet_recombination_breakpoints():
	parameter_sets = [
		(False, {'genmap':'tests/data/recombination_breaks.map'}),
		(True, {'recombrate':1000000}),
		(False, {'recombrate':.0000001})
	]
	
	for expect_recombination, parameters in parameter_sets:
		with TemporaryDirectory() as tempdir:
			outvcf = tempdir + '/output-recombination_breaks.vcf'
			outlist = tempdir + '/output.recomb'
			run_whatshap(phase_input_files=[recombination_breaks_bamfile], variant_file='tests/data/quartet.vcf', output=outvcf,
					ped='tests/data/recombination_breaks.ped', recombination_list_filename = outlist, **parameters)
			assert os.path.isfile(outvcf)

			tables = list(VcfReader(outvcf))
			assert len(tables) == 1
			table = tables[0]
			assert table.chromosome == '1'
			assert len(table.variants) == 4
			assert table.samples == ['HG002', 'HG005', 'HG003', 'HG004']
			assert table.num_of_blocks_of('HG002') == 0
			assert table.num_of_blocks_of('HG005') == 0
			assert table.num_of_blocks_of('HG003') == 1
			assert table.num_of_blocks_of('HG004') == 0

			phase0 = VariantCallPhase(68735304, 0, None)
			phase1 = VariantCallPhase(68735304, 1, None)

			assert_phasing(table.phases_of('HG002'), [None, None, None, None])
			assert_phasing(table.phases_of('HG005'), [None, None, None, None])
			if expect_recombination:
				assert_phasing(table.phases_of('HG003'), [phase0, phase0, None, phase1])
			else:
				assert_phasing(table.phases_of('HG003'), [phase0, phase0, None, phase0])
			assert_phasing(table.phases_of('HG004'), [None, None, None, None])
			
			lines = open(outlist).readlines()
			if expect_recombination:
				assert len(lines) == 3
				assert lines[1]=='HG002 1 68735433 68738308 0 0 0 1 3\n'
				assert lines[2]=='HG005 1 68735433 68738308 0 0 0 1 3\n'
			else:
				assert len(lines) == 1


def test_phase_trio_zero_distance():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_whatshap(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/zero-genetic-distance.map')
		assert os.path.isfile(outvcf)


def test_phase_quartet_recombination_breakpoints():
	parameter_sets = [
		(False, {'genmap':'tests/data/recombination_breaks.map'}),
		(True, {'recombrate':1000000}),
		(False, {'recombrate':.0000001})
	]


def test_haplotag():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		outbam = tempdir + '/output.bam'
		run_whatshap(phase_input_files=[recombination_breaks_bamfile], variant_file='tests/data/quartet.vcf', output=outvcf, ped='tests/data/recombination_breaks.ped')
		run_haplotag(variant_file=outvcf, alignment_file=recombination_breaks_bamfile, output=outbam)
		ps_count = 0
		for alignment in pysam.AlignmentFile(outbam):
			if alignment.has_tag('PS'):
				ps_count += 1
		assert ps_count > 0


def test_compare1():
	with TemporaryDirectory() as tempdir:
		outtsv = tempdir + '/output.tsv'
		run_compare(vcf=['tests/data/phased1.vcf', 'tests/data/phased2.vcf'], names='p1,p2', tsv_pairwise=outtsv, sample='sample1')
		lines = [l.split('\t') for l in open(outtsv)]
		assert len(lines) == 3
		Fields = namedtuple('Fields', [ f.strip('#\n') for f in lines[0] ])
		entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

		assert entry_chrA.chromosome == 'chrA'
		assert entry_chrA.all_assessed_pairs == '4'
		assert entry_chrA.all_switches == '1'
		assert entry_chrA.all_switchflips == '1/0'
		assert entry_chrA.largestblock_assessed_pairs == '2'
		assert entry_chrA.largestblock_switches == '1'
		assert entry_chrA.largestblock_hamming == '1'

		assert entry_chrB.chromosome == 'chrB'
		assert entry_chrB.all_assessed_pairs == '1'
		assert entry_chrB.all_switches == '0'
		assert entry_chrB.all_switchflips == '0/0'
		assert entry_chrB.largestblock_assessed_pairs == '1'
		assert entry_chrB.largestblock_switches == '0'
		assert entry_chrB.largestblock_hamming == '0'


def test_compare2():
	with TemporaryDirectory() as tempdir:
		outtsv = tempdir + '/output.tsv'
		run_compare(vcf=['tests/data/phased1.vcf', 'tests/data/phased2.vcf'], names='p1,p2', tsv_pairwise=outtsv, sample='sample2')
		lines = [l.split('\t') for l in open(outtsv)]
		assert len(lines) == 3
		Fields = namedtuple('Fields', [ f.strip('#\n') for f in lines[0] ])
		entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

		assert entry_chrA.chromosome == 'chrA'
		assert entry_chrA.all_assessed_pairs == '6'
		assert entry_chrA.all_switches == '2'
		assert entry_chrA.all_switchflips == '0/1'
		assert entry_chrA.largestblock_assessed_pairs == '5'
		assert entry_chrA.largestblock_switches == '2'
		assert entry_chrA.largestblock_hamming == '1'

		assert entry_chrB.chromosome == 'chrB'
		assert entry_chrB.all_assessed_pairs == '1'
		assert entry_chrB.all_switches == '1'
		assert entry_chrB.all_switchflips == '1/0'
		assert entry_chrB.largestblock_assessed_pairs == '1'
		assert entry_chrB.largestblock_switches == '1'
		assert entry_chrB.largestblock_hamming == '1'


def test_compare_only_snps():
	with TemporaryDirectory() as tempdir:
		outtsv = tempdir + '/output.tsv'
		run_compare(vcf=['tests/data/phased1.vcf', 'tests/data/phased2.vcf'], names='p1,p2', tsv_pairwise=outtsv, sample='sample2', only_snps=True)
		lines = [l.split('\t') for l in open(outtsv)]
		assert len(lines) == 3
		Fields = namedtuple('Fields', [ f.strip('#\n') for f in lines[0] ])
		entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

		assert entry_chrA.chromosome == 'chrA'
		assert entry_chrA.all_assessed_pairs == '3'
		assert entry_chrA.all_switches == '2'
		assert entry_chrA.all_switchflips == '0/1'
		assert entry_chrA.largestblock_assessed_pairs == '3'
		assert entry_chrA.largestblock_switches == '2'
		assert entry_chrA.largestblock_hamming == '1'

		assert entry_chrB.chromosome == 'chrB'
		assert entry_chrB.all_assessed_pairs == '1'
		assert entry_chrB.all_switches == '1'
		assert entry_chrB.all_switchflips == '1/0'
		assert entry_chrB.largestblock_assessed_pairs == '1'
		assert entry_chrB.largestblock_switches == '1'
		assert entry_chrB.largestblock_hamming == '1'
