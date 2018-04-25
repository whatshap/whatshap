from tempfile import TemporaryDirectory
import os
import pysam
from nose.tools import raises
import math
import vcf

from whatshap.genotype import run_genotype
from whatshap.vcf import VcfReader

trio_bamfile = 'tests/data/trio.pacbio.bam'
trio_merged_bamfile = 'tests/data/trio-merged-blocks.bam'
trio_paired_end_bamfile = 'tests/data/paired_end.sorted.bam'
ped_samples_bamfile = 'tests/data/ped_samples.bam'
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


def test_one_variant():
	run_genotype(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null')


def test_default_output():
	"""Output to stdout"""
	run_genotype(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf')


def test_bam_without_readgroup():
	run_genotype(phase_input_files=['tests/data/no-readgroup.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', ignore_read_groups=True)


@raises(SystemExit)
def test_requested_sample_not_found():
	run_genotype(phase_input_files=['tests/data/oneread.bam'], variant_file='tests/data/onevariant.vcf',
		output='/dev/null', samples=['DOES_NOT_EXIST'])


def test_with_reference():
	run_genotype(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
		reference='tests/data/pacbio/reference.fasta')

def test_no_indels():
	with TemporaryDirectory() as tempdir:
		for priors in [[False, tempdir+'/priors.vcf'], [True,None]]:
			outvcf = tempdir + '/output_gl.vcf'
			run_genotype(phase_input_files=['tests/data/pacbio/pacbio.bam'], variant_file='tests/data/pacbio/variants.vcf',
				reference='tests/data/pacbio/reference.fasta', output=outvcf, indels=False, nopriors=priors[0], prioroutput=priors[1])

			result_vcfs = [outvcf]
			if priors[0]:
				result_vcfs.append(tempdir + '/priors.vcf')

			# make sure indels not genotyped (also in priors.vcf if computed)
			for o_vcf in result_vcfs:
				vcf_reader = vcf.Reader(filename=o_vcf)
				default_l = math.log10(1/3.0)

				for record in vcf_reader:
					if len(record.ALT[0]) > 1:
						for call in record.samples:
							GL = getattr(call.data, 'GL', None)
							print('GL:', GL, record.POS)
							assert(GL == [default_l, default_l, default_l])


def likeliest_genotype(a, b, c, thres):
	prob_a = 10**a
	prob_b = 10**b
	prob_c = 10**c

	to_sort = [(prob_a,0),(prob_b,1),(prob_c,2)]
	to_sort.sort(key=lambda x: x[0])

	if (to_sort[2][0] > to_sort[1][0]) and (to_sort[2][0] > thres):
		return to_sort[2][1]
	else:
		return '.'

def test_GtQualThreshold():
	for threshold in [0,2,3,6,13,50]:
		with TemporaryDirectory() as tempdir:
			thres = 1-10**(-threshold/10.0)

			out_vcf = tempdir + '/out.vcf'
			priors_vcf = tempdir + '/priors.vcf'
			run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf',
				output=out_vcf, gt_qual_threshold=threshold, indels=False, prioroutput=priors_vcf)

			for out in [open(out_vcf,'r'), open(priors_vcf,'r')]:
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
					print(likelihoods[0], likelihoods[1], likelihoods[2], gt, genotype, thres)
					assert(gt == genotype)

def test_genotyping_one_of_three_individuals():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		outpriors = tempdir + '/priors.vcf'
		run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf, samples=['HG003'],
		prioroutput=outpriors)

		for outfile in [outvcf, outpriors]:
			assert os.path.isfile(outfile)

			tables = list(VcfReader(outfile, phases=True,genotype_likelihoods=True))
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

def test_use_ped_samples():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_ped_samples.vcf'
		run_genotype(phase_input_files=[ped_samples_bamfile], variant_file='tests/data/ped_samples.vcf', output=outvcf,
			ped='tests/data/trio.ped', genmap='tests/data/trio.map', use_ped_samples=True)
		assert os.path.isfile(outvcf)

		tables = list(VcfReader(outvcf, phases=True, genotype_likelihoods=True))
		assert len(tables) == 1
		table = tables[0]
		assert table.chromosome == '1'
		assert len(table.variants) == 5
		assert table.samples == ['HG004', 'HG003', 'HG002', 'orphan']

		default_l = math.log10(1/3.0)
		for var in table.genotype_likelihoods_of('orphan'):
			assert(var.log10_probs() == (default_l, default_l, default_l))

def test_ped_sample():
	with TemporaryDirectory() as tempdir:
		# running with --ped and --sample on subset of trio, should give same results as running with only --sample
		# the trio information should be ignored
		outvcf1 = tempdir + '/output1.vcf'
		outvcf2 = tempdir + '/output2.vcf'
		for sample_set in [['HG002'], ['HG003'], ['HG004'], ['HG002','HG003'], ['HG002','HG004'], ['HG003','HG004']]:
			run_genotype(phase_input_files=[ped_samples_bamfile], variant_file='tests/data/ped_samples.vcf', output=outvcf1,
				ped='tests/data/trio.ped', samples=sample_set)
			run_genotype(phase_input_files=[ped_samples_bamfile], variant_file='tests/data/ped_samples.vcf', output=outvcf2,
				samples=sample_set)

			assert os.path.isfile(outvcf1)
			assert os.path.isfile(outvcf2)
			tables1 = list(VcfReader(outvcf1, phases=True, genotype_likelihoods=True))
			tables2 = list(VcfReader(outvcf2, phases=True, genotype_likelihoods=True))
			assert( (len(tables1) == 1) and (len(tables2) == 1) )
			table1, table2 = tables1[0], tables2[0]

			for individual in sample_set:
				for var1,var2 in zip(table1.genotype_likelihoods_of(individual), table2.genotype_likelihoods_of(individual)):
					print(var1, var2)
					assert(var1.log10_probs() == var2.log10_probs())

def test_genotyping_trio():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		outpriors = tempdir + 'priors.vcf'
		run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map', prioroutput=outpriors)

		for outfile in [outvcf,outpriors]:
			assert os.path.isfile(outfile)

			tables = list(VcfReader(outfile, phases=True))
			assert len(tables) == 1
			table = tables[0]
			assert table.chromosome == '1'
			assert len(table.variants) == 5
			assert table.samples == ['HG004', 'HG003', 'HG002']

def test_genotyping_specific_chromosome():
	for requested_chromosome in ['1','2']:
		with TemporaryDirectory() as tempdir:
			outvcf = tempdir + '/output.vcf'
			outpriors = tempdir + '/priors.vcf'
			run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio-two-chromosomes.vcf', output=outvcf,
					ped='tests/data/trio.ped', genmap='tests/data/trio.map', chromosomes=[requested_chromosome],
					prioroutput=outpriors)

			for outfile in [outvcf, outpriors]:
				assert os.path.isfile(outfile)

				tables = list(VcfReader(outfile, genotype_likelihoods=True))

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

def test_genotype_likelihoods_given():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_gl.vcf'
		run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio_genotype_likelihoods.vcf', output=outvcf,
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

# GL field was already present, make sure it is replaced by new likelihoods
def test_genotype_log_likelihoods_given():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_gl_log.vcf'
		outpriors = tempdir + '/priors.vcf'
		run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio_genotype_log_likelihoods.vcf', output=outvcf,
		        ped='tests/data/trio.ped', genmap='tests/data/trio.map', gt_qual_threshold=0, prioroutput=outpriors)

		for outfile in [outvcf, outpriors]:
			assert os.path.isfile(outfile)

			tables = list(VcfReader(outfile, phases=True, genotype_likelihoods=True))
			assert len(tables) == 1
			table = tables[0]
			assert table.chromosome == '1'
			assert len(table.variants) == 5
			assert table.samples == ['HG004', 'HG003', 'HG002']

			# check if GL likelihoods were replaced
			vcf_reader = vcf.Reader(filename=outfile)
			print(vcf_reader.samples, outfile)
			for record in vcf_reader:
				for call in record.samples:
					GL = getattr(call.data, 'GL', None)
					GQ = getattr(call.data, 'GQ', None)
					print('GL:', GL, 'GQ', GQ)
					assert(GL != [-1,-1,-1])
					assert(GQ != 100)

def test_empty_format_field():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output_empty_format.vcf'
		run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/empty_format.vcf', output=outvcf, gt_qual_threshold=0)

		# check if sample fields now contain information
		assert os.path.isfile(outvcf)
		vcf_reader = vcf.Reader(filename=outvcf)
		for record in vcf_reader:
			print(record.samples, outvcf)
			assert(len(record.samples) == 3)

def test_phase_trio_paired_end_reads():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output-paired_end.vcf'
		run_genotype(phase_input_files=[trio_paired_end_bamfile], variant_file='tests/data/paired_end.sorted.vcf', output=outvcf,
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
		run_genotype(phase_input_files=[short_bamfile],
			ignore_read_groups=True,
			variant_file='tests/data/short-genome/wrongchromosome.vcf', output=outvcf)

def extract_likelihoods(line):
	entries = line.split()
	likelihood_str = entries[9].split(':')[-1:][0]
	likelihoods = [10.0**float(i) for i in likelihood_str.split(',')]
	return likelihoods

def test_adding_constant():
	for const in [0.1,0.2,0.3,0.5,0.7,1,2,5,10,20,100]:
		with TemporaryDirectory() as tempdir:
			priors_raw_vcf = tempdir + '/output.raw_priors.vcf'
			outvcf_raw_vcf = tempdir + '/output_raw.vcf'
			priors_const_vcf = tempdir + '/output.const_priors.vcf'
			outvcf_const_vcf = tempdir + '/output_raw.vcf'

			# run genotyping without adding constant to priors
			run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf',
				prioroutput=priors_raw_vcf, output=outvcf_raw_vcf, indels=False)

			# run genotyping with modified priors
			run_genotype(phase_input_files=[trio_bamfile], variant_file='tests/data/trio.vcf',
				prioroutput=priors_const_vcf, output=outvcf_const_vcf, indels=False, constant=const)

			# check if priors were modified properly
			priors_raw = open(priors_raw_vcf, 'r')
			priors_const = open(priors_const_vcf, 'r')

			priors_raw.seek(0)
			priors_const.seek(0)

			lines_raw = [line for line in priors_raw.readlines() if not line.startswith('#')]
			lines_const = [line for line in priors_const.readlines() if not line.startswith('#')]

			assert(len(lines_raw) == len(lines_const))

			for i in range(len(lines_raw)):
				likelihoods_raw = extract_likelihoods(lines_raw[i])
				likelihoods_const = extract_likelihoods(lines_const[i])

				norm_sum = likelihoods_raw[0]+likelihoods_raw[1]+likelihoods_raw[2] + 3.0*const
				#print('raw likelihoods: ', likelihoods_raw, ' modified likelihoods: ', likelihoods_const)

				print(float(likelihoods_const[0]),float((likelihoods_raw[0] + const)/norm_sum))
				print(float(likelihoods_const[1]),float((likelihoods_raw[1] + const)/norm_sum))
				print(float(likelihoods_const[2]),float((likelihoods_raw[2] + const)/norm_sum))

				for j in range(3):
					assert(round(likelihoods_const[j],10) == round((likelihoods_raw[j] + const)/norm_sum,10))

			priors_raw.close()
			priors_const.close()
