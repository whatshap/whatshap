from tempfile import TemporaryDirectory
import os
from io import StringIO
import pysam
from nose.tools import raises
from collections import namedtuple
import math
import vcf

from whatshap.genotyping import run_genotyping
from whatshap.phase import run_whatshap
from whatshap.haplotag import run_haplotag
from whatshap.hapcut2vcf import run_hapcut2vcf
from whatshap.compare import run_compare
from whatshap.vcf import VcfReader, VariantCallPhase, GenotypeLikelihoods

trio_bamfile = 'tests/data/trio.pacbio.bam'
trio_merged_bamfile = 'tests/data/trio-merged-blocks.bam'
trio_paired_end_bamfile = 'tests/data/paired_end.sorted.bam'
recombination_breaks_bamfile = 'tests/data/recombination_breaks.sorted.bam'
quartet2_bamfile = 'tests/data/quartet2.bam'
short_bamfile = 'tests/data/short-genome/short.bam'
indels_bamfile = 'tests/data/indels.bam'

bam_files = [trio_bamfile, trio_merged_bamfile, trio_paired_end_bamfile, recombination_breaks_bamfile, quartet2_bamfile, short_bamfile, indels_bamfile]


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
	run_genotyping(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null')


def test_default_output():
	"""Output to stdout"""
	run_genotyping(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf')


def test_bam_without_readgroup():
	run_genotyping(phase_input_files=['tests/data/no-readgroup.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', ignore_read_groups=True)


@raises(SystemExit)
def test_requested_sample_not_found():
	run_genotyping(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', samples=['DOES_NOT_EXIST'])


def test_with_reference():
	run_genotyping(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
		reference='tests/data/pacbio/reference.fasta')


def test_no_indels():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_gl.vcf'
		run_genotyping(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
			reference='tests/data/pacbio/reference.fasta', output=outvcf, indels=False)
		
		# make sure indels not genotyped
		vcf_reader = vcf.Reader(filename=outvcf)
		default_l = math.log10(1/3.0)

		for record in vcf_reader:
			if len(record.ALT[0]) > 1:
				for call in record.samples:
					GL = getattr(call.data, 'GL', None)
					print('GL:', GL, record.POS)
					assert(GL == [default_l, default_l, default_l])
				
		
		
def likeliest_genotype(a, b, c, thres):
	prob_a = 10 ** a
	prob_b = 10 ** b
	prob_c = 10 ** c
	
	to_sort = [(prob_a,0),(prob_b,1),(prob_c,2)]
	to_sort.sort(key=lambda x: x[0])

	if (to_sort[2][0] > to_sort[1][0]) and (to_sort[2][0] > thres):
		return to_sort[2][1]
	else:
		return '.'
	
def test_GtQualThreshold():
	for threshold in [0,2,3,6,13,50]:
		thres = 1-10**(-threshold/10.0)	
	
		out = StringIO()
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf',
			output=out, gt_qual_threshold=threshold, indels=False)
		out.seek(0)

		lines = [line for line in out.readlines() if not line.startswith('#')]
	
		for line in lines:
			entries = line.split()
			likelihood_str = entries[9].split(':')[-1:][0]
			likelihoods = [float(i) for i in likelihood_str.split(',')]
			genotype = entries[9].split(':')[0]
			if not genotype == '.':
				genotype = int(genotype[0]) + int(genotype[2])
			
			gt = likeliest_genotype(likelihoods[0], likelihoods[1], likelihoods[2], thres)
			print(10**likelihoods[0], 10**likelihoods[1], 10**likelihoods[2], gt, genotype, thres)
			assert(gt == genotype)

def test_genotyping_one_of_three_individuals():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf, samples=['HG003'])
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True,genotype_likelihoods=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']
		
		# there should be no genotype predicitons for HG003/HG002		
		default_l = math.log10(1/3.0)
		for l in [table.genotype_likelihoods_of('HG002'), table.genotype_likelihoods_of('HG004')]:
			for var in l:
				assert(var.log10_probs() == (default_l, default_l, default_l))

def test_genotyping_trio():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

def test_phase_specific_chromosome():
	for requested_chromosome in ['1','2']:
		with TemporaryDirectory() as tempdir:
			outvcf = tempdir + '/output.vcf'
			run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio-two-chromosomes.vcf', output=outvcf,
					ped='tests/data/trio.ped', genmap='tests/data/trio.map', chromosomes=[requested_chromosome])
			assert os.path.isfile(outvcf)

			tables = list(VcfReader(outvcf, genotype_likelihoods=True))

			assert len(tables) == 2
			for table in tables:
				assert len(table.variants) == 5
				assert table.samples == ['HG004', 'HG003', 'HG002']

			index = 0
			if requested_chromosome == '1':
				index = 1
				
			# should be no genotype likelihoods for skipped chromosomes
			for s in tables[index].samples:
				tables[index].genotype_likelihoods_of(s) == [None] * 5
				tables[not index].genotype_likelihoods_of(s) != [None] * 5
				
# given likelihoods should be deleted 
def test_genotype_likelihoods_given():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_gl.vcf'
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio_genotype_likelihoods.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True, genotype_likelihoods=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']

		# check if PL likelihoods (that were present before) are deleted
		vcf_reader = vcf.Reader(filename=outvcf)
		print(vcf_reader.samples, outvcf)
		for record in vcf_reader:
			for call in record.samples:
				PL = getattr(call.data, 'PL', None)
				GL = getattr(call.data, 'GL', None)
				print('GL:', GL, 'PL:', PL)
				assert(PL==None)
				assert(GL != None)
				
# GL field is already present, make sure it is replaced by new likelihoods
def test_genotype_log_likelihoods_given():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_gl_log.vcf'
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/trio_genotype_log_likelihoods.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map', gt_qual_threshold=0)
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True, genotype_likelihoods=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002']
		print("DONE")
		# check if GL likelihoods were replaced
		vcf_reader = vcf.Reader(filename=outvcf)
		print(vcf_reader.samples, outvcf)
		for record in vcf_reader:
			for call in record.samples:
				GL = getattr(call.data, 'GL', None)
				GQ = getattr(call.data, 'GQ', None)
				print('GL:', GL, 'GQ', GQ)
				assert(GL != [-1,-1,-1])
				assert(GQ != 100)
				
def test_no_gt_field():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_no_gt.vcf'
		run_genotyping(phase_input_files=[trio_bamfile], variant_file='tests/data/TEST_no_gt.vcf', output=outvcf, samples=['HG004'])

def test_phase_trio_paired_end_reads():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output-paired_end.vcf'
		run_genotyping(phase_input_files=[trio_paired_end_bamfile], variant_file='tests/data/paired_end.sorted.vcf', output=outvcf,
		        ped='tests/data/trio_paired_end.ped', genmap='tests/data/trio.map')
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 3
		assert table.samples == ['mother', 'father', 'child']

@raises(SystemExit)
def test_wrong_chromosome():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_genotyping(phase_input_files=[short_bamfile],
			ignore_read_groups=True,
			variant_file='tests/data/short-genome/wrongchromosome.vcf', output=outvcf)
