#!/usr/bin/env python3
"""
Genotype variants

Runs only the genotyping algorithm. Genotype Likelihoods are computed using the
forward backward algorithm.
"""
import logging
import sys
import platform
import resource
from collections import defaultdict

from xopen import xopen

from contextlib import ExitStack
from .vcf import VcfReader, GenotypeVcfWriter
from . import __version__
from .core import ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods, GenotypeDPTable, compute_genotypes
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
			load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import AlignmentFileNotIndexedError, EmptyAlignmentFileError
from .timer import StageTimer
from .variants import ReadSetReader
from .utils import IndexedFasta, FastaNotIndexedError

from .phase import read_reads, select_reads, split_input_file_list, setup_pedigree


logger = logging.getLogger(__name__)


def determine_genotype(likelihoods, threshold_prob):
	"""given genotype likelihoods for 0/0,0/1,1/1, determines likeliest genotype"""

	to_sort = [(likelihoods[0],0),(likelihoods[1],1),(likelihoods[2],2)]
	to_sort.sort(key=lambda x: x[0])

	# make sure there is a unique maximum which is greater than the threshold
	if (to_sort[2][0] > to_sort[1][0]) and (to_sort[2][0] > threshold_prob):
		return to_sort[2][1]
	else:
		return -1


def run_genotype(phase_input_files, variant_file, reference=None,
		output=sys.stdout, samples=None, chromosomes=None,
		ignore_read_groups=False, indels=True, mapping_quality=20,
		max_coverage=15, nopriors=False,
		ped=None, recombrate=1.26, genmap=None, gt_qual_threshold=0,
		prioroutput=None, constant=0.0, overhang=10,affine_gap=False, gap_start=10, gap_extend=7, mismatch=15,
		write_command_line_header=True, use_ped_samples=False):
	"""
	For now: this function only runs the genotyping algorithm. Genotype likelihoods for
	all variants are computed using the forward backward algorithm
	"""
	print('running only genotyping algorithm')
	timers = StageTimer()
	timers.start('overall')
	logger.info("This is WhatsHap (genotyping) %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:

		# read the given input files (bams,vcfs,ref...)
		numeric_sample_ids = NumericSampleIds()
		phase_input_bam_filenames, phase_input_vcf_filenames = split_input_file_list(phase_input_files)
		try:
			readset_reader = stack.enter_context(ReadSetReader(
				phase_input_bam_filenames, reference, numeric_sample_ids, mapq_threshold=mapping_quality,
				overhang=overhang, affine=affine_gap, gap_start=gap_start, gap_extend=gap_extend, default_mismatch=mismatch))
		except OSError as e:
			logger.error(e)
			sys.exit(1)
		except AlignmentFileNotIndexedError as e:
			logger.error('The file %r is not indexed. Please create the appropriate BAM/CRAM '
				'index with "samtools index"', str(e))
			sys.exit(1)
		except EmptyAlignmentFileError as e:
			logger.error('No reads could be retrieved from %r. If this is a CRAM file, possibly the '
				'reference could not be found. Try to use --reference=... or check you '
				'$REF_PATH/$REF_CACHE settings', str(e))
			sys.exit(1)
		try:
			phase_input_vcf_readers = [VcfReader(f, indels=indels, phases=True) for f in phase_input_vcf_filenames]
		except OSError as e:
			logger.error(e)
			sys.exit(1)
		if reference:
			try:
				fasta = stack.enter_context(IndexedFasta(reference))
			except OSError as e:
				logger.error('%s', e)
				sys.exit(1)
			except FastaNotIndexedError as e:
				logger.error('An index file (.fai) for the reference %r could not be found. '
					'Please create one with "samtools faidx".', str(e))
				sys.exit(1)
		else:
			fasta = None
		del reference
		if isinstance(output, str):
			output = stack.enter_context(xopen(output, 'w'))
		if write_command_line_header:
			command_line = '(whatshap {}) {}'.format(__version__ , ' '.join(sys.argv[1:]))
		else:
			command_line=None

		# vcf writer for final genotype likelihoods
		vcf_writer = GenotypeVcfWriter(command_line=command_line, in_path=variant_file,
		        out_file=output)
		# vcf writer for only the prior likelihoods (if output is desired)
		prior_vcf_writer = None
		if prioroutput != None:
			prior_vcf_writer = GenotypeVcfWriter(command_line=command_line, in_path=variant_file, out_file=open(prioroutput,'w'))

		# parse vcf with input variants
		# remove all likelihoods that may already be present
		vcf_reader = VcfReader(variant_file, indels=indels, genotype_likelihoods=False, ignore_genotypes=True)

		if ignore_read_groups and not samples and len(vcf_reader.samples) > 1:
			logger.error('When using --ignore-read-groups on a VCF with '
				'multiple samples, --sample must also be used.')
			sys.exit(1)
		if not samples:
			samples = vcf_reader.samples

		# if --use-ped-samples is set, use only samples from PED file
		if ped and use_ped_samples:
			samples = set()
			for trio in PedReader(ped, numeric_sample_ids):
				if (trio.child is None or trio.mother is None or trio.father is None):
					continue
				samples.add(trio.mother)
				samples.add(trio.father)
				samples.add(trio.child)

		vcf_sample_set = set(vcf_reader.samples)
		for sample in samples:
			if sample not in vcf_sample_set:
				logger.error('Sample %r requested on command-line not found in VCF', sample)
				sys.exit(1)

		samples = frozenset(samples)
		# list of all trios across all families
		all_trios = dict()

		# Keep track of connected components (aka families) in the pedigree
		family_finder = ComponentFinder(samples)

		# if pedigree information present, parse it
		if ped:
			all_trios, pedigree_samples = setup_pedigree(ped, numeric_sample_ids, samples)
			if genmap:
				logger.info('Using region-specific recombination rates from genetic map %s.', genmap)
			else:
				logger.info('Using uniform recombination rate of %g cM/Mb.', recombrate)
			for trio in all_trios:
				family_finder.merge(trio.mother, trio.child)
				family_finder.merge(trio.father, trio.child)

		# map family representatives to lists of family members
		families = defaultdict(list)
		for sample in samples:
			families[family_finder.find(sample)].append(sample)
		# map family representatives to lists of trios for this family
		family_trios = defaultdict(list)
		for trio in all_trios:
			family_trios[family_finder.find(trio.child)].append(trio)
		largest_trio_count = max([0] + [len(trio_list) for trio_list in family_trios.values()])
		logger.info('Working on %d samples from %d famil%s', len(samples), len(families), 'y' if len(families)==1 else 'ies')

		if max_coverage + 2 * largest_trio_count > 25:
			logger.warning('The maximum coverage is too high! '
				'WhatsHap may take a long time to finish and require a huge amount of memory.')

		# Read phase information provided as VCF files, if provided.
		phase_input_vcfs = []
		timers.start('parse_phasing_vcfs')
		for reader, filename in zip(phase_input_vcf_readers, phase_input_vcf_filenames):
			# create dict mapping chromsome names to VariantTables
			m = dict()
			logger.info('Reading phased blocks from %r', filename)
			for variant_table in reader:
				m[variant_table.chromosome] = variant_table
			phase_input_vcfs.append(m)
		timers.stop('parse_phasing_vcfs')

		# compute genotype likelihood threshold
		gt_prob = 1.0-(10 ** (-gt_qual_threshold/10.0))

		timers.start('parse_vcf')
		for variant_table in vcf_reader:

			# create a mapping of genome positions to indices
			var_to_pos = dict()
			for i in range(len(variant_table.variants)):
				var_to_pos[variant_table.variants[i].position] = i

			chromosome = variant_table.chromosome
			timers.stop('parse_vcf')
			if (not chromosomes) or (chromosome in chromosomes):
				logger.info('======== Working on chromosome %r', chromosome)
			else:
				logger.info('Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)', chromosome)
				vcf_writer.write_genotypes(chromosome,variant_table,indels,leave_unchanged=True)
				if prioroutput != None:
					prior_vcf_writer.write_genotypes(chromosome,variant_table,indels,leave_unchanged=True)
				continue

			positions = [v.position for v in variant_table.variants]
			if not nopriors:
				# compute prior genotype likelihoods based on all reads
				for sample in samples:
					logger.info('---- Initial genotyping of %s', sample)
					with timers('read_bam'):
						bam_sample = None if ignore_read_groups else sample
						readset, vcf_source_ids = read_reads(readset_reader, chromosome, variant_table.variants, bam_sample, fasta, [], numeric_sample_ids, phase_input_bam_filenames)
						readset.sort()
						genotypes, genotype_likelihoods = compute_genotypes(readset, positions)
						# recompute genotypes based on given threshold
						reg_genotype_likelihoods = []
						for gl in range(len(genotype_likelihoods)):
							norm_sum = genotype_likelihoods[gl][0] + genotype_likelihoods[gl][1] + genotype_likelihoods[gl][2] + 3*constant
							regularized = ((genotype_likelihoods[gl][0]+constant)/norm_sum, (genotype_likelihoods[gl][1]+constant)/norm_sum, (genotype_likelihoods[gl][2]+constant)/norm_sum)
							genotypes[gl] = determine_genotype(regularized, gt_prob)
							reg_genotype_likelihoods.append(regularized)
						variant_table.set_genotype_likelihoods_of(sample, [PhredGenotypeLikelihoods(*gl) for gl in reg_genotype_likelihoods])
						variant_table.set_genotypes_of(sample, genotypes)
			else:

				# use uniform genotype likelihoods for all individuals
				for sample in samples:
					variant_table.set_genotype_likelihoods_of(sample, [PhredGenotypeLikelihoods(1.0/3.0,1.0/3.0,1.0/3.0)] * len(positions))

			# if desired, output the priors in separate vcf
			if prioroutput != None:
				prior_vcf_writer.write_genotypes(chromosome,variant_table,indels)

			# Iterate over all families to process, i.e. a separate DP table is created
			# for each family.
			for representative_sample, family in sorted(families.items()):
				if len(family) == 1:
					logger.info('---- Processing individual %s', representative_sample)
				else:
					logger.info('---- Processing family with individuals: %s', ','.join(family))
				max_coverage_per_sample = max(1, max_coverage // len(family))
				logger.info('Using maximum coverage per sample of %dX', max_coverage_per_sample)
				trios = family_trios[representative_sample]
				assert (len(family) == 1) or (len(trios) > 0)

				# Get the reads belonging to each sample
				readsets = dict()
				for sample in family:
					with timers('read_bam'):
						bam_sample = None if ignore_read_groups else sample
						readset, vcf_source_ids = read_reads(readset_reader, chromosome, variant_table.variants, bam_sample, fasta, phase_input_vcfs, numeric_sample_ids, phase_input_bam_filenames)

					with timers('select'):
						selected_reads = select_reads(readset, max_coverage_per_sample, preferred_source_ids = vcf_source_ids)
					readsets[sample] = selected_reads

				# Merge reads into one ReadSet (note that each Read object
				# knows the sample it originated from).
				all_reads = ReadSet()
				for sample, readset in readsets.items():
					for read in readset:
						assert read.is_sorted(), "Add a read.sort() here"
						all_reads.add(read)

				all_reads.sort()

				# Determine which variants can (in principle) be phased
				accessible_positions = sorted(all_reads.get_positions())
				logger.info('Variants covered by at least one phase-informative '
					'read in at least one individual after read selection: %d',
					len(accessible_positions))

				# Create Pedigree
				pedigree = Pedigree(numeric_sample_ids)
				for sample in family:
					# genotypes are assumed to be unknown, so ignore information that
					# might already be present in the input vcf
					all_genotype_likelihoods = variant_table.genotype_likelihoods_of(sample)
					genotype_l = [ all_genotype_likelihoods[var_to_pos[a_p]] for a_p in accessible_positions]
					pedigree.add_individual(sample, [3] * len(accessible_positions), genotype_l)
				for trio in trios:
					pedigree.add_relationship(
						mother_id=trio.mother,
						father_id=trio.father,
						child_id=trio.child)

				if genmap:
					# Load genetic map
					recombination_costs = recombination_cost_map(load_genetic_map(genmap), accessible_positions)
				else:
					recombination_costs = uniform_recombination_map(recombrate, accessible_positions)

				# Finally, run genotyping algorithm
				with timers('genotyping'):
					problem_name = 'genotyping'
					logger.info('Genotype %d sample%s by solving the %s problem ...',
						len(family), 's' if len(family) > 1 else '', problem_name)
					forward_backward_table = GenotypeDPTable(numeric_sample_ids, all_reads, recombination_costs, pedigree, accessible_positions)
					# store results
					for s in family:
						likelihood_list = variant_table.genotype_likelihoods_of(s)
						genotypes_list = variant_table.genotypes_of(s)

						for pos in range(len(accessible_positions)):
							likelihoods = forward_backward_table.get_genotype_likelihoods(s,pos)

							# compute genotypes from likelihoods and store information
							geno = determine_genotype(likelihoods, gt_prob)
							genotypes_list[var_to_pos[accessible_positions[pos]]] = geno
							likelihood_list[var_to_pos[accessible_positions[pos]]] = likelihoods

						variant_table.set_genotypes_of(s, genotypes_list)
						variant_table.set_genotype_likelihoods_of(s,likelihood_list)

			with timers('write_vcf'):
				logger.info('======== Writing VCF')
				vcf_writer.write_genotypes(chromosome,variant_table,indels)
				logger.info('Done writing VCF')

			logger.debug('Chromosome %r finished', chromosome)
			timers.start('parse_vcf')
		timers.stop('parse_vcf')

	logger.info('\n== SUMMARY ==')
	timers.stop('overall')
	if sys.platform == 'linux':
		memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		logger.info('Maximum memory usage: %.3f GB', memory_kb / 1E6)
	logger.info('Time spent reading BAM:                      %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF:                      %6.1f s', timers.elapsed('parse_vcf'))
	if len(phase_input_vcfs) > 0:
		logger.info('Time spent parsing input phasings from VCFs: %6.1f s', timers.elapsed('parse_phasing_vcfs'))
	logger.info('Time spent selecting reads:                  %6.1f s', timers.elapsed('select'))
	logger.info('Time spent genotyping:                          %6.1f s', timers.elapsed('genotyping'))
	logger.info('Time spent writing VCF:                      %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent on rest:                          %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:                          %6.1f s', timers.elapsed('overall'))


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be gzip-compressed)')
	arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
		help='BAM or VCF file(s) with phase information, either through sequencing reads (BAM) or through phased blocks (VCF)')

	arg('-o', '--output', default=sys.stdout,
		help='Output VCF file. Add .gz to the file name to get compressed output. '
			'If omitted, use standard output.')
	arg('--reference', '-r', metavar='FASTA',
		help='Reference file. Provide this to detect alleles through re-alignment. '
			'If no index (.fai) exists, it will be created')

	arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
	arg('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	arg('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	arg('--indels', dest='indels', default=False, action='store_true',
		help='Also genotype indels (default: genotype only SNPs)')
	arg('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM header and assume all reads come '
		'from the same sample.')
	arg('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
		help='Name of a sample to genotype. If not given, all samples in the '
		'input VCF are genotyped. Can be used multiple times.')
	arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
		help='Name of chromosome to genotyped. If not given, all chromosomes in the '
		'input VCF are genotyped. Can be used multiple times.')
	arg('--gt-qual-threshold', metavar='GTQUALTHRESHOLD', type=float, default=0,
		help='Phred scaled error probability threshold used for genotyping (default: %(default)s). Must be at least 0. '
		'If error probability of genotype is higher, genotype ./. is output.')
	arg('--no-priors', dest='nopriors', default=False, action='store_true',
		help='Skip initial prior genotyping and use uniform priors (default: %(default)s).')
	arg('-p', '--prioroutput', default=None,
		help='output prior genotype likelihoods to the given file.')
	arg('--overhang', metavar='OVERHANG', default=10, type=int,
		help='When --reference is used, extend alignment by this many bases to left and right when realigning (default: %(default)s).')
	arg('--constant', metavar='CONSTANT', default=0, type=float,
		help='This constant is used to regularize the priors (default: %(default)s).')
	arg('--affine-gap', default=False, action='store_true',
		help='When detecting alleles through re-alignment, use affine gap costs (EXPERIMENTAL).')
	arg('--gap-start', metavar='GAPSTART', default=10, type=float,
		help='gap starting penalty in case affine gap costs are used (default: %(default)s).')
	arg('--gap-extend', metavar='GAPEXTEND', default=7, type=float,
		help='gap extend penalty in case affine gap costs are used (default: %(default)s).')
	arg('--mismatch', metavar='MISMATCH', default=15, type=float,
		help='mismatch cost in case affine gap costs are used (default: %(default)s)')

	arg = parser.add_argument_group('Pedigree genotyping').add_argument
	arg('--ped', metavar='PED/FAM',
		help='Use pedigree information in PED file to improve genotyping '
		'(switches to PedMEC algorithm). Columns 2, 3, 4 must refer to child, '
		'mother, and father sample names as used in the VCF and BAM. Other '
		'columns are ignored (EXPERIMENTAL).')
	arg('--recombrate', metavar='RECOMBRATE', type=float, default=1.26,
		help='Recombination rate in cM/Mb (used with --ped). If given, a constant recombination '
		'rate is assumed (default: %(default)gcM/Mb).')
	arg('--genmap', metavar='FILE',
		help='File with genetic map (used with --ped) to be used instead of constant recombination '
		'rate, i.e. overrides option --recombrate.')
	arg('--use-ped-samples', dest='use_ped_samples',
		action='store_true', default=False,
		help='Only work on samples mentioned in the provided PED file.')


def validate(args, parser):
	if args.ignore_read_groups and args.ped:
		parser.error('Option --ignore-read-groups cannot be used together with --ped')
	if args.genmap and not args.ped:
		parser.error('Option --genmap can only be used together with --ped')
	if args.genmap and (len(args.chromosomes) != 1):
		parser.error('Option --genmap can only be used when working on exactly one chromosome (use --chromosome)')
	if len(args.phase_input_files) == 0:
		parser.error('Not providing any PHASEINPUT files not allowed for genotyping.')
	if args.gt_qual_threshold < 0:
		parser.error('Genotype quality threshold (gt-qual-threshold) must be at least 0.')
	if args.prioroutput != None and args.nopriors:
		parser.error('Genotype priors are only computed if --no-priors is NOT set.')
	if args.constant != 0  and args.nopriors:
		parser.error('--constant can only be used if --no-priors is NOT set..')
	if args.affine_gap and not args.reference:
		parser.error('Option --affine-gap can only be used together with --reference.')
	if args.use_ped_samples and not args.ped:
		parser.error('Option --use-ped-samples can only be used when PED file is provided (--ped).')
	if args.use_ped_samples and args.samples:
		parser.error('Option --use-ped-samples cannot be used together with --samples')


def main(args):
	run_genotype(**vars(args))
