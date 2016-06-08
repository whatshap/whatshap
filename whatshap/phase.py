#!/usr/bin/env python3
"""
Phase variants in a VCF with the WhatsHap algorithm

Read a VCF and one or more files with phase information (BAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
"""
import logging
import sys
import platform
from collections import defaultdict

import pyfaidx

try:
	from contextlib import ExitStack
except ImportError:
	from contextlib2 import ExitStack  # PY32
from .vcf import VcfReader, PhasedVcfWriter, VariantTable
from . import __version__
from .args import HelpfulArgumentParser as ArgumentParser
from .core import ReadSet, DPTable, readselection, Pedigree, PedigreeDPTable, NumericSampleIds
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map)
from .bam import BamIndexingError, SampleNotFoundError, HaplotypeBamWriter
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError

__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)


def find_components(phased_positions, reads, master_block=None):
	"""
	Return a dict that maps each variant position to the component it is in.
	Variants are considered to be in the same component if a read exists that
	covers both. A component is identified by the position of its leftmost
	variant.
	master_block -- List of positions in a "master block", i.e. all blocks containing
	                any of these positions are merged into one block.
	"""
	logger.debug('Finding connected components ...')
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [ variant.position for variant in read if variant.position in phased_positions ]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	if not master_block is None:
		for position in master_block[1:]:
			component_finder.merge(master_block[0], position)
	components = { position : component_finder.find(position) for position in phased_positions }
	return components


def find_largest_component(components):
	"""
	Determine the largest component and return a sorted list of positions
	contained in it.
	components -- dictionary mapping positin to block_id as returned by find_components.
	"""
	blocks = defaultdict(list)
	for position, block_id in components.items():
		blocks[block_id].append(position)
	largest = []
	for block in blocks.values():
		if len(block) > len(largest):
			largest = block
	largest.sort()
	return largest


def best_case_blocks(reads):
	"""
	Given a list of core reads, determine the number of phased blocks that
	would result if each variant were actually phased.

	Return the number of connected components and non-singleton components.
	"""
	positions = set()
	for read in reads:
		for variant in read:
			positions.add(variant.position)
	component_finder = ComponentFinder(positions)
	for read in reads:
		read_positions = [ variant.position for variant in read ]
		for position in read_positions[1:]:
			component_finder.merge(read_positions[0], position)
	# A dict that maps each component to the number of variants it contains
	component_sizes = defaultdict(int)
	for position in positions:
		component_sizes[component_finder.find(position)] += 1
	non_singletons = [ component for component, size in component_sizes.items() if size > 1]
	return len(component_sizes), len(non_singletons)


def read_reads(readset_reader, chromosome, variants, sample, fasta, phase_input_vcfs, numeric_sample_ids, phase_input_bam_filenames):
	"""Return a sorted ReadSet"""
	logger.info('Reading alignments for sample %r and detecting alleles ...', sample)
	reference = fasta[chromosome] if fasta else None
	try:
		readset = readset_reader.read(chromosome, variants, sample, reference)
	except SampleNotFoundError:
		logger.warning("Sample %r not found in any BAM file.", sample)
		readset = ReadSet()
	except ReadSetError as e:
		logger.error("%s", e)
		sys.exit(1)

	# Add phasing information from VCF files, if present
	for i, phase_input_vcf in enumerate(phase_input_vcfs):
		if chromosome in phase_input_vcf:
			vt = phase_input_vcf[chromosome]
			source_id = len(phase_input_bam_filenames) + i
			for read in vt.phased_blocks_as_reads(sample, variants, source_id, numeric_sample_ids[sample]):
				readset.add(read)

	# TODO is this necessary?
	for read in readset:
		read.sort()
	readset.sort()

	logger.info('Found %d reads covering %d variants', len(readset), len(readset.get_positions()))
	return readset


def select_reads(readset, max_coverage):
	readset = readset.subset([i for i, read in enumerate(readset) if len(read) >= 2])
	logger.info('Kept %d reads that cover at least two variants each', len(readset))
	logger.info('Reducing coverage to at most %dX by selecting most informative reads ...', max_coverage)
	selected_indices = readselection(readset, max_coverage)
	selected_reads = readset.subset(selected_indices)
	logger.info('Selected %d reads covering %d variants',
		len(selected_reads), len(selected_reads.get_positions()))

	return selected_reads


def phase_sample(sample, chromosome, reads, all_heterozygous, max_coverage, timers, stats, haplotype_bam_writer, numeric_sample_ids):
	"""
	Phase variants of a single sample on a single chromosome.
	"""
	with timers('select'):
		selected_reads = select_reads(reads, max_coverage)

	n_best_case_blocks, n_best_case_nonsingleton_blocks = best_case_blocks(reads)
	n_best_case_blocks_cov, n_best_case_nonsingleton_blocks_cov = best_case_blocks(selected_reads)
	stats.n_best_case_blocks += n_best_case_blocks
	stats.n_best_case_nonsingleton_blocks += n_best_case_nonsingleton_blocks
	stats.n_best_case_blocks_cov += n_best_case_blocks_cov
	stats.n_best_case_nonsingleton_blocks_cov += n_best_case_nonsingleton_blocks_cov
	logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
		n_best_case_nonsingleton_blocks, n_best_case_blocks)
	logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
		n_best_case_nonsingleton_blocks_cov, n_best_case_blocks_cov)

	with timers('phase'):
		logger.info('Phasing the variants (using %d reads)...', len(selected_reads))
		if all_heterozygous:
			# For the all heterozygous case we use a PedigreeDPTable, which is more memory efficient.
			# Once implemented, this should also be done for the "not all heterozygous" (="distrust genotypes")
			# case, see Issue #77.

			# all genotypes are heterozygous
			accessible_positions = selected_reads.get_positions()
			genotypes = [1] * len(accessible_positions)
			# create pedigree with only one sample
			pedigree = Pedigree(numeric_sample_ids)
			pedigree.add_individual(sample, genotypes)
			# recombination costs are zero
			recombination_costs = [0] * len(accessible_positions)
			# Run the core algorithm: construct DP table ...
			dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree)
			# ... and do the backtrace to get the solution
			superreads_list, transmission_vector = dp_table.get_super_reads()
			superreads = superreads_list[0]
		else:
			# Run the core algorithm: construct DP table ...
			dp_table = DPTable(selected_reads, all_heterozygous)
			# ... and do the backtrace to get the solution
			superreads = dp_table.get_super_reads()
		logger.info('MEC score of phasing: %d', dp_table.get_optimal_cost())

		n_homozygous = sum(1 for v1, v2 in zip(*superreads)
			if v1.allele == v2.allele and v1.allele in (0, 1))
		stats.n_homozygous += n_homozygous

	with timers('components'):
		# The variant.allele attribute can be either 0 (major allele), 1 (minor allele),
		# or 3 (equal scores). If all_heterozygous is on (default), we can get
		# the combinations 0/1, 1/0 and 3/3 (the latter means: unphased).
		# If all_heterozygous is off, we can also get all other combinations.
		# In both cases, we are interested only in 0/1 and 1/0.
		allowed = frozenset([(0, 1), (1, 0)])
		phased_positions = [ v1.position for v1, v2 in zip(*superreads)
			if (v1.allele, v2.allele) in allowed ]
		components = find_components(phased_positions, selected_reads)
		logger.info('No. of variants considered for phasing: %d', len(superreads[0]))
		logger.info('No. of variants that were phased: %d', len(phased_positions))

	n_phased_blocks = len(set(components.values()))
	stats.n_phased_blocks += n_phased_blocks
	logger.info('No. of phased blocks: %d', n_phased_blocks)
	if all_heterozygous:
		assert n_homozygous == 0
	else:
		logger.info('No. of heterozygous variants determined to be homozygous: %d', n_homozygous)

	if haplotype_bam_writer is not None:
		logger.info('Writing used reads to haplotype-specific BAM files')
		haplotype_bam_writer.write(selected_reads, dp_table.get_optimal_partitioning(), chromosome)

	return superreads, components


class UnknownInputFileError(Exception):
	pass


def split_input_file_list(input_files):
	bams = []
	vcfs = []
	#TODO: maybe take a peek at the content rather than determining file type based on filename ending.
	for filename in input_files:
		if filename.endswith('.bam'):
			bams.append(filename)
		elif filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
			vcfs.append(filename)
		else:
			raise UnknownInputFileError('Unable to determine type of input file '+filename)
	return bams, vcfs


def setup_pedigree(ped_path, numeric_sample_ids, samples):
	"""
	Read in PED file to set up list of relationships (individuals).

	Return a pair (individuals, pedigree_samples)

	ped_path -- path to PED file
	samples -- samples that exist in the VCF
	"""
	individuals = []
	pedigree_samples = set()
	for individual in PedReader(ped_path, numeric_sample_ids):
		if (individual.id is None or individual.mother_id is None or
				    individual.father_id is None):
			logger.warning('Relationship %s/%s/%s ignored '
			               'because at least one of the individuals is unknown',
			               individual.id, individual.mother_id, individual.father_id)
		else:
			individuals.append(individual)
			pedigree_samples.add(individual.id)
			pedigree_samples.add(individual.mother_id)
			pedigree_samples.add(individual.father_id)

	for sample in pedigree_samples:
		if sample not in samples:
			# TODO should that really be an error?
			logger.error('Sample %r not found in VCF', sample)
			sys.exit(1)
	for sample in samples:
		if sample not in pedigree_samples:
			# TODO should be single-individual-phased instead
			# or perhaps it does work with the PedMEC algorithm
			logger.warning('No relationship known for sample %r - '
			               'will not be phased', sample)
	return individuals, pedigree_samples


def run_whatshap(phase_input_files, variant_file, reference=None,
		output=None, samples=None, chromosomes=None, ignore_read_groups=False, indels=True,
		mapping_quality=20, max_coverage=15, all_heterozygous=True,
		haplotype_bams_prefix=None, ped=None, recombrate=1.26, genmap=None, genetic_haplotyping=True):
	"""
	Run WhatsHap.

	phase_input_files -- list of paths to BAM/VCF files
	variant_file -- path to input VCF
	reference -- path to reference FASTA
	output -- path to output VCF or use None for stdout
	samples -- names of samples to phase. an empty list means: phase all samples
	chromosomes -- names of chromosomes to phase. an empty list means: phase all chromosomes
	ignore_read_groups
	mapping_quality -- discard reads below this mapping quality
	max_coverage
	all_heterozygous
	genetic_haplotyping -- in ped mode, merge disconnected blocks based on genotype status
	"""
	class Statistics:
		pass
	stats = Statistics()
	timers = StageTimer()
	timers.start('overall')
	stats.n_homozygous = 0
	stats.n_phased_blocks = 0
	stats.n_best_case_blocks = 0
	stats.n_best_case_nonsingleton_blocks = 0
	stats.n_best_case_blocks_cov = 0
	stats.n_best_case_nonsingleton_blocks_cov = 0
	logger.info("This is WhatsHap %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		numeric_sample_ids = NumericSampleIds()
		phase_input_bam_filenames, phase_input_vcf_filenames = split_input_file_list(phase_input_files)
		try:
			readset_reader = stack.enter_context(ReadSetReader(phase_input_bam_filenames, numeric_sample_ids, mapq_threshold=mapping_quality))
		except (OSError, BamIndexingError) as e:
			logger.error(e)
			sys.exit(1)
		try:
			phase_input_vcf_readers = [VcfReader(f, indels=indels) for f in phase_input_vcf_filenames]
		except OSError as e:
			logger.error(e)
			sys.exit(1)
		if reference:
			try:
				fasta = stack.enter_context(pyfaidx.Fasta(reference, as_raw=True))
			except OSError as e:
				logger.error('%s', e)
				sys.exit(1)
		else:
			fasta = None
		del reference
		# no reference given -> do not re-align -> must normalize indels
		normalize = fasta is None
		if output is not None:
			output = stack.enter_context(open(output, 'w'))
		else:
			output = sys.stdout
		command_line = '(whatshap {}) {}'.format(__version__ , ' '.join(sys.argv[1:]))
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=variant_file, normalized=normalize, out_file=output)
		vcf_reader = VcfReader(variant_file, indels=indels)

		if ignore_read_groups and not samples and len(vcf_reader.samples) > 1:
			logger.error('When using --ignore-read-groups on a VCF with '
				'multiple samples, --sample must also be used.')
			sys.exit(1)
		if not samples:
			samples = vcf_reader.samples
		vcf_sample_set = set(vcf_reader.samples)
		for sample in samples:
			if sample not in vcf_sample_set:
				logger.error('Sample %r requested on command-line not found in VCF', sample)
				sys.exit(1)

		haplotype_bam_writer = None
		if haplotype_bams_prefix is not None:
			logger.warning('Writing haplotype BAMs only for the first sample')
			haplotype_bam_writer = HaplotypeBamWriter(phase_input_bam_filenames, haplotype_bams_prefix, samples[0])

		samples = frozenset(samples)

		if ped:
			individuals, pedigree_samples = setup_pedigree(ped, numeric_sample_ids, vcf_reader.samples)
			if genmap:
				logger.info('Using region-specific recombination rates from genetic map %s.', genmap)
			else:
				logger.info('Using uniform recombination rate of %g cM/Mb.', recombrate)

		# Read phase information provided as VCF files, if provided.
		# TODO: do this chromosome- and/or sample-wise on demand to save memory.
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

		timers.start('parse_vcf')
		for variant_table in vcf_reader:
			chromosome = variant_table.chromosome
			timers.stop('parse_vcf')
			if (not chromosomes) or (chromosome in chromosomes):
				logger.info('Working on chromosome %r', chromosome)
			else:
				logger.info('Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)', chromosome)
				with timers('write_vcf'):
					superreads, components = dict(), dict()
					vcf_writer.write(chromosome, superreads, components)
				continue
			# These two variables hold the phasing results for all samples
			superreads, components = dict(), dict()
			if ped:
				# variant indices with at least one missing genotype
				missing_genotypes = set()
				# variant indices with at least one Mendelian conflict
				mendelian_conflicts = set()
				# variant indices with at least one heterozygous genotype
				heterozygous = set()
				# variant indices with at least one homozygous genotype
				homozygous = set()
				for trio in individuals:
					# TODO fix attribute names of Individual class
					genotypes_mother = variant_table.genotypes_of(trio.mother_id)
					genotypes_father = variant_table.genotypes_of(trio.father_id)
					genotypes_child = variant_table.genotypes_of(trio.id)

					for index, (gt_mother, gt_father, gt_child) in enumerate(zip(
							genotypes_mother, genotypes_father, genotypes_child)):
						is_missing = False
						for gt in (gt_mother, gt_father, gt_child):
							if gt == -1:
								missing_genotypes.add(index)
								is_missing = True
							elif gt == 1:
								heterozygous.add(index)
							else:
								assert gt in [0,2]
								homozygous.add(index)
						if not is_missing:
							if mendelian_conflict(gt_mother, gt_father, gt_child):
								mendelian_conflicts.add(index)

				# retain variants that are heterozygous in at least one individual (anywhere in the pedigree)
				# and do not have neither missing genotypes nor Mendelian conflicts
				to_retain = heterozygous.difference(missing_genotypes).difference(mendelian_conflicts)
				# discard every variant that is not to be retained
				to_discard = set(range(len(variant_table))).difference(to_retain)

				# Determine positions of selected variants that are homozygous in at least one individual.
				# These are used later to merge blocks containing these variants into one block (since
				# the are conntected by "genetic haplotyping").
				homozygous_positions = [variant_table.variants[i].position for i in to_retain.intersection(homozygous)]

				# Remove calls where *any* trio has a mendelian conflict or
				# is homozygous in all three individuals
				variant_table.remove_rows_by_index(to_discard)

				logger.info('Number of variants skipped due to missing genotypes: %d', len(missing_genotypes))
				logger.info('Number of variants skipped due to Mendelian conflicts: %d', len(mendelian_conflicts))
				logger.info('Number of remaining variants heterozygous in at least one individual: %d', len(variant_table))

				# Get the reads belonging to each sample
				readsets = dict()  # TODO this could become a list
				for sample in variant_table.samples:
					with timers('read_bam'):
						readset = read_reads(readset_reader, chromosome, variant_table.variants, sample, fasta, phase_input_vcfs, numeric_sample_ids, phase_input_bam_filenames)

					# TODO: Read selection done w.r.t. all variants, where using heterozygous variants only
					# TODO: would probably give better results.
					with timers('select'):
						selected_reads = select_reads(readset, max_coverage)
					readsets[sample] = selected_reads

				accessible_positions = []
				for readset in readsets.values():
					accessible_positions.extend(readset.get_positions())
				accessible_positions = sorted(set(accessible_positions))
				logger.info('Variants covered by at least one phase-informative '
					'read in at least one individual after read selection: %d',
					len(accessible_positions))

				# Keep only accessible positions
				variant_table.subset_rows_by_position(accessible_positions)
				assert len(variant_table.variants) == len(accessible_positions)

				# Create Pedigree
				individual_ids = { sample: index for index, sample in enumerate(pedigree_samples) }
				pedigree = Pedigree(numeric_sample_ids)
				for sample in pedigree_samples:
					pedigree.add_individual(sample, variant_table.genotypes_of(sample))
				for individual in individuals:
					pedigree.add_relationship(
						mother_id=individual.mother_id,
						father_id=individual.father_id,
						child_id=individual.id)

				# Merge reads into one ReadSet (note that each Read object
				# knows the sample it originated from).
				all_reads = ReadSet()
				for sample, readset in readsets.items():
					for read in readset:
						assert read.is_sorted(), "Add a read.sort() here"
						all_reads.add(read)

				all_reads.sort()

				if genmap:
					# Load genetic map
					recombination_costs = recombination_cost_map(load_genetic_map(genmap), accessible_positions)
				else:
					recombination_costs = uniform_recombination_map(recombrate, accessible_positions)

				# Finally, run phasing algorithm
				with timers('phase'):
					logger.info('Phasing %d samples with the PedMEC algorithm ...',
						len(pedigree_samples))
					dp_table = PedigreeDPTable(all_reads, recombination_costs, pedigree)
					superreads_list, transmission_vector = dp_table.get_super_reads()
					logger.info('PedMEC cost: %d', dp_table.get_optimal_cost())

				with timers('components'):
					master_block = None
					if genetic_haplotyping:
						master_block = sorted(set(homozygous_positions).intersection(set(accessible_positions)))
					overall_components = find_components(accessible_positions, all_reads, master_block)
					n_phased_blocks = len(set(overall_components.values()))
					stats.n_phased_blocks += n_phased_blocks
					logger.info('No. of phased blocks: %d', n_phased_blocks)
					largest_component = find_largest_component(overall_components)
					if len(largest_component) > 0:
						logger.info('Largest component contains %d variants (%.1f%% of accessible variants) between position %d and %d', len(largest_component), len(largest_component)*100.0/len(accessible_positions), largest_component[0]+1, largest_component[-1]+1)

				if False:
					n_recombination = find_recombination(transmission_vector, overall_components, accessible_positions, recombcost, recombination_list_filename)
					logger.info('No. of detected recombination events: %d', n_recombination)

				# TODO Do superreads actually come out in the order in which the
				# individuals were added to the pedigree?
				for sample, sample_superreads in zip(pedigree_samples, superreads_list):
					superreads[sample] = sample_superreads
					# identical for all samples
					components[sample] = overall_components
			else:
				for sample, genotypes in zip(variant_table.samples, variant_table.genotypes):
					if sample not in samples:
						continue
					# pick variants heterozygous in this sample
					variants = [ v for v, gt in zip(variant_table.variants, genotypes) if gt == 1 ]
					logger.info('Found %d heterozygous variants on sample %r', len(variants), sample)
					bam_sample = None if ignore_read_groups else sample
					with timers('read_bam'):
						reads = read_reads(readset_reader, chromosome, variants, bam_sample, fasta, phase_input_vcfs, numeric_sample_ids, phase_input_bam_filenames)

					sample_superreads, sample_components = phase_sample(
						sample, chromosome, reads, all_heterozygous, max_coverage, timers, stats, haplotype_bam_writer, numeric_sample_ids)
					superreads[sample] = sample_superreads
					components[sample] = sample_components
			with timers('write_vcf'):
				vcf_writer.write(chromosome, superreads, components)
			logger.debug('Chromosome %r finished', chromosome)
			timers.start('parse_vcf')
		timers.stop('parse_vcf')

	logger.info('== SUMMARY ==')
	# TODO: Print more meaningful summary, including block sizes, mendelian conflicts, etc.
	if not ped:
		logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
			stats.n_best_case_nonsingleton_blocks, stats.n_best_case_blocks)
		logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
			stats.n_best_case_nonsingleton_blocks_cov, stats.n_best_case_blocks_cov)
	if all_heterozygous:
		assert stats.n_homozygous == 0
	else:
		logger.info('No. of heterozygous variants determined to be homozygous: %d', stats.n_homozygous)
	timers.stop('overall')
	logger.info('Time spent reading BAM:                      %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF:                      %6.1f s', timers.elapsed('parse_vcf'))
	if len(phase_input_vcfs) > 0:
		logger.info('Time spent parsing input phasings from VCFs: %6.1f s', timers.elapsed('parse_phasing_vcfs'))
	logger.info('Time spent selecting reads:                  %6.1f s', timers.elapsed('select'))
	logger.info('Time spent phasing:                          %6.1f s', timers.elapsed('phase'))
	logger.info('Time spent writing VCF:                      %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent finding components:               %6.1f s', timers.elapsed('components'))
	logger.info('Time spent on rest:                          %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:                          %6.1f s', timers.elapsed('overall'))


def add_arguments(parser):
	parser.add_argument('--version', action='version', version=__version__)
	parser.add_argument('-o', '--output', default=None,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('--reference', '-r', metavar='FASTA',
		help='Reference file. Provide this to detect alleles through re-alignment. '
			'If no index (.fai) exists, it will be created')
	parser.add_argument('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	parser.add_argument('--indels', dest='indels', default=False, action='store_true',
		help='Also phase indels (default: do not phase indels)')
	parser.add_argument('--distrust-genotypes', dest='all_heterozygous',
		action='store_false', default=True,
		help='Allow switching variants from hetero- to homozygous in an '
		'optimal solution (see documentation).')
	parser.add_argument('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM header and assume all reads come '
		'from the same sample.')
	parser.add_argument('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
		help='Name of a sample to phase. If not given, all samples in the '
		'input VCF are phased. Can be used multiple times.')
	parser.add_argument('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
		help='Name of chromosome to phase. If not given, all chromosomes in the '
		'input VCF are phased. Can be used multiple times.')
	parser.add_argument('--haplotype-bams', metavar='PREFIX', dest='haplotype_bams_prefix', default=None,
		help='Write reads that have been used for phasing to haplotype-specific BAM files. '
		'Creates PREFIX.1.bam and PREFIX.2.bam')
	parser.add_argument('--ped', metavar='PED/FAM',
		help='Use pedigree information in PED file to improve phasing '
		'(switches to PedMEC algorithm). Columns 2, 3, 4 must refer to child, '
		'mother, and father sample names as used in the VCF and BAM. Other '
		'columns are ignored.')
	parser.add_argument('--recombrate', metavar='RECOMBRATE', type=float, default=1.26,
		help='Recombination rate in cM/Mb (used with --ped). If given, a constant recombination '
		'rate is assumed (default: %(default)gcM/Mb).')
	parser.add_argument('--genmap', metavar='GENMAP',
		help='File with genetic map (used with --ped) to be used instead of constant recombination '
		'rate, i.e. overrides option --recombrate.')  # TODO describe what the file format is
	parser.add_argument('--no-genetic-haplotyping', dest='genetic_haplotyping',
		action='store_false', default=True,
		help='Do not merge blocks that are not connected by reads (i.e. solely based on genotype '
		'status). Default: when in --ped mode, merge all blocks that contain at least one '
		'homozygous genotype in at least one individual into one block.')
	parser.add_argument('variant_file', metavar='VCF', help='VCF file with variants to be phased (can be gzip-compressed)')
	parser.add_argument('phase_input_files', nargs='+', metavar='PHASEINPUT',
	    help='BAM or VCF file(s) with phase information, either through sequencing reads (BAM) or through phased blocks (VCF)')


def validate(args, parser):
	if args.ped and not args.all_heterozygous:
		parser.error('Option --distrust-genotypes cannot be used together with --ped')
	if args.ignore_read_groups and args.ped:
		parser.error('Option --ignore-read-groups cannot be used together with --ped')
	if args.genmap and not args.ped:
		parser.error('Option --genmap can only be used together with --ped')
	if args.genmap and (len(args.chromosomes) != 1):
		parser.error('Option --genmap can only be used when working on exactly one chromosome (use --chromosome)')
	if args.ped and args.samples:
		parser.error('Option --sample cannot be used together with --ped')


def main(args):
	run_whatshap(**vars(args))
