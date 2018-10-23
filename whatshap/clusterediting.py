"""
Compute transformation columnwise using cluster editing and compute consensus clustering solving MEC.

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
For each column, cluster all reads covering in (cluster editing), and compute
a consensus clustering using MEC.

"""
import sys
import logging
import platform
import resource

from collections import defaultdict
from copy import deepcopy
from math import log

from xopen import xopen
from networkx import Graph, number_of_nodes, number_of_edges, connected_components, node_connected_component, shortest_path

from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter, VcfGenotypeLikelihoods
from . import __version__
from .core import Read, ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, GenotypeLikelihoods, Genotype, compute_genotypes
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import AlignmentFileNotIndexedError, SampleNotFoundError, ReferenceNotFoundError, EmptyAlignmentFileError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError
from .utils import detect_file_format, IndexedFasta, FastaNotIndexedError
from .matrixtransformation import MatrixTransformation
from .phase import read_reads, select_reads, split_input_file_list, setup_pedigree, find_components, find_largest_component, write_read_list
from .dissimilarityplots import draw_plots

__author__ = "Jana Ebler" 

logger = logging.getLogger(__name__)

def print_readset(readset):
	result = ""
	positions = readset.get_positions()
	for read in readset:
		result += read.name + '\t' + '\t' + '\t'
		for pos in positions:
			if pos in read:
				# get corresponding variant
				for var in read:
					if var.position == pos:
						result += str(var.allele)
			else:
				result += ' '
		result += '\n'
	print(result)

def run_clusterediting(
	phase_input_files,
	variant_file,
	ploidy,
	reference=None,
	output=sys.stdout,
	samples=None,
	chromosomes=None,
	ignore_read_groups=False,
	indels=True,
	mapping_quality=20,
	tag='PS',
	write_command_line_header=True,
	read_list_filename=None,
	errorrate = 0.1, 
	min_overlap = 5
	):
	"""
	Run Polyploid Phasing.
	
	phase_input_files -- list of paths to BAM/CRAM/VCF files
	variant-file -- path to input VCF
	reference -- path to reference FASTA
	output -- path to output VCF or a file like object
	samples -- names of samples to phase. An empty list means: phase all samples
	chromosomes -- names of chromosomes to phase. An empty list means: phase all chromosomes
	ignore_read_groups
	mapping_quality -- discard reads below this mapping quality
	tag -- How to store phasing info in the VCF, can be 'PS' or 'HP'
	write_command_line_header -- whether to add a ##commandline header to the output VCF
	"""
	timers = StageTimer()
	timers.start('overall')
	logger.info("This is WhatsHap (polyploid) %s running under Python %s", __version__, platform.python_version())
	with ExitStack() as stack:
		numeric_sample_ids = NumericSampleIds()
		phase_input_bam_filenames, phase_input_vcf_filenames = split_input_file_list(phase_input_files)
		assert len(phase_input_bam_filenames) > 0
		try:
			readset_reader = stack.enter_context(ReadSetReader(phase_input_bam_filenames, reference,
				numeric_sample_ids, mapq_threshold=mapping_quality))
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
			command_line = '(whatshap {}) {}'.format(__version__, ' '.join(sys.argv[1:]))
		else:
			command_line = None
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=variant_file,
			out_file=output, tag=tag, ploidy=ploidy)
		# TODO for now, assume we always trust the genotypes
		vcf_reader = VcfReader(variant_file, indels=indels, genotype_likelihoods=False, ploidy=ploidy)

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

		samples = frozenset(samples)

		read_list_file = None
		if read_list_filename:
			read_list_file = create_read_list_file(read_list_filename)
		
		timers.start('parse_vcf')
		for variant_table in vcf_reader:
			chromosome = variant_table.chromosome
			timers.stop('parse_vcf')
			if (not chromosomes) or (chromosome in chromosomes):
				logger.info('======== Working on chromosome %r', chromosome)
			else:
				logger.info('Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)', chromosome)
				with timers('write_vcf'):
					superreads, components = dict(), dict()
					vcf_writer.write(chromosome, superreads, components)
				continue
			# These two variables hold the phasing results for all samples
			superreads, components = dict(), dict()

			# Iterate over all samples to process
			for sample in samples:
				logger.info('---- Processing individual %s', sample)
				missing_genotypes = set()
				heterozygous = set()
				homozygous = set()

				genotypes = variant_table.genotypes_of(sample)
				for index, gt in enumerate(genotypes):
					if gt.is_none():
						missing_genotypes.add(index)
					elif not gt.is_homozygous():
						heterozygous.add(index)
					else:
						assert gt.is_homozygous()

				to_discard = set(range(len(variant_table))).difference(heterozygous)
				phasable_variant_table = deepcopy(variant_table)

				# Remove calls to be discarded from variant table
				phasable_variant_table.remove_rows_by_index(to_discard)

				logger.info('Number of variants skipped due to missing genotypes: %d', len(missing_genotypes))
				logger.info('Number of remaining heterozygous variants: %d', len(phasable_variant_table))

				# Get the reads belonging to this sample
				bam_sample = None if ignore_read_groups else sample
				readset, vcf_source_ids = read_reads(readset_reader, chromosome, phasable_variant_table.variants, bam_sample, fasta, [], numeric_sample_ids, phase_input_bam_filenames)
				readset.sort()
				# TODO: len == min_overlap ?
				readset = readset.subset([i for i, read in enumerate(readset) if len(read) >= min_overlap])
				logger.info('Kept %d reads that cover at least two variants each', len(readset))

				# transform the allele matrix
				selected_reads = select_reads(readset, 5*ploidy, preferred_source_ids = vcf_source_ids)
				readset = selected_reads
				print_readset(readset)
				transformation = MatrixTransformation(readset, find_components(readset.get_positions(), readset), ploidy, errorrate, min_overlap)
				transformed_matrix = transformation.get_transformed_matrix()
				cluster_counts = transformation.get_cluster_counts()

				# TODO: remove print
				print_readset(transformed_matrix)
				# TODO: remove print
#				print('cluster counts:', cluster_counts)

				# TODO include this readselection step?
				selected_reads = select_reads(transformed_matrix, 7, preferred_source_ids = vcf_source_ids)
				transformed_matrix = selected_reads

				print_readset(transformed_matrix)
				# solve MEC to get overall partitioning, prepare input objects for this
				cluster_pedigree = Pedigree(numeric_sample_ids, ploidy)
				positions = transformed_matrix.get_positions()
				alleles = [i for i in range(0,ploidy)]
				cluster_pedigree.add_individual(sample, [Genotype(alleles)]*len(positions), [GenotypeLikelihoods(ploidy, ploidy,[0])]*len(positions))
				recombination_costs = uniform_recombination_map(1.26, positions)
				partitioning_dp_table = PedigreeDPTable(transformed_matrix, recombination_costs, cluster_pedigree, ploidy, False, cluster_counts, positions)
				optimal_partitioning = partitioning_dp_table.get_optimal_partitioning()
				print('Consensus MEC cost: ', partitioning_dp_table.get_optimal_cost())

				print('read partitioning: ', optimal_partitioning, partitioning_dp_table.get_optimal_cost())

				# printing
				clu_to_r = defaultdict(list)
				for read, partition in zip(transformed_matrix,optimal_partitioning):
					clu_to_r[partition].append(read.name)
				
#				for c,l in clu_to_r.items():
#					print(c,l)

				# TODO the order of the reads in clusters_per_window and readset can differ. Therefore, reorder read_partitioning
#				read_to_partition = {}
#				for i,read in enumerate(transformed_matrix):
#					read_to_partition[read.name] = read_partitioning[i]
#				optimal_partitioning =  [ read_to_partition[r.name] for r in transformed_matrix]	

				selected_names = [r.name for r in transformed_matrix]
				# keep only those reads in original readset that have been selected in transformed matrix
				to_keep = [ i for i,read in enumerate(readset) if read.name in selected_names]
				readset = readset.subset(to_keep)
				print_readset(readset)

				# prepare input for determining the best allele configuration
				# Determine which variants can (in principle) be phased
				accessible_positions = sorted(readset.get_positions())
				logger.info('Variants covered by at least one phase-informative read: %d', len(accessible_positions))

				# Keep only accessible positions
				phasable_variant_table.subset_rows_by_position(accessible_positions)
				assert len(phasable_variant_table.variants) == len(accessible_positions)

				pedigree = Pedigree(numeric_sample_ids, ploidy)
				ind_genotypes = [gt.get_genotype() for gt in phasable_variant_table.genotypes_of(sample)]
				pedigree.add_individual(sample, ind_genotypes, None)

				# For the given partitioning of the reads, determine best allele configurations
				with timers('phase'):
					logger.info('Phasing %s by determining best allele assignment ... ', sample)
					recombination_costs = uniform_recombination_map(1.26, accessible_positions)
					allele_assignment = PedigreeDPTable(readset, recombination_costs, pedigree, ploidy, False, None, accessible_positions, optimal_partitioning)
					superreads_list = allele_assignment.get_super_reads()
					logger.info('MEC cost: %d', allele_assignment.get_optimal_cost())
				with timers('components'):
					overall_components = find_components(accessible_positions, readset, None, None)
					n_phased_blocks = len(set(overall_components.values()))
					logger.info('No. of phased blocks: %d', n_phased_blocks)
					largest_component = find_largest_component(overall_components)
					if len(largest_component) > 0:
							logger.info('Largest component contains %d variants (%.1f%% of accessible variants) between position %d and %d', len(largest_component), len(largest_component)*100.0/len(accessible_positions), largest_component[0]+1, largest_component[-1]+1)

				assert(len(superreads_list) == 2)
				sample_superreads = superreads_list[0]
				superreads[sample] = sample_superreads[0]
				assert len(sample_superreads[0]) == ploidy
				sr_sample_id = sample_superreads[0][0].sample_id
				for sr in sample_superreads[0]:
					assert sr.sample_id == sr_sample_id == numeric_sample_ids[sample]
				components[sample] = overall_components

#				if read_list_file:
#					write_read_list(all_reads, allele_assignment.get_optimal_partitioning(), components, numeric_sample_ids, read_list_file)
			with timers('write_vcf'):
				logger.info('======== Writing VCF')
				changed_genotypes = vcf_writer.write(chromosome, superreads, components)
				logger.info('Done writing VCF')
				assert len(changed_genotypes) == 0
			logger.debug('Chromosome %r finished', chromosome)
			timers.start('parse_vcf')
		timers.stop('parse_vcf')
	
	if read_list_file:
		read_list_file.close()

	logger.info('\n== SUMMARY ==')
	timers.stop('overall')
	if sys.platform == 'linux':
		memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		logger.info('Maximum memory usage: %.3f GB', memory_kb / 1E6)
	logger.info('Time spent reading BAM/CRAM:                 %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF:                      %6.1f s', timers.elapsed('parse_vcf'))
	logger.info('Time spent selecting reads:                  %6.1f s', timers.elapsed('select'))
	logger.info('Time spent pruning readset:                  %6.1f s', timers.elapsed('prune'))
	logger.info('Time spent phasing:                          %6.1f s', timers.elapsed('phase'))
	logger.info('Time spent writing VCF:                      %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent finding components:               %6.1f s', timers.elapsed('components'))
	logger.info('Time spent on rest:                          %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:                          %6.1f s', timers.elapsed('overall'))


def add_arguments(parser):
	arg = parser.add_argument
	# Positional argument
	arg('variant_file', metavar='VCF',
		help='VCF file with variants to be phased (can be gzip-compressed)')
	arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
		help='BAM or CRAM with sequencing reads.')
	arg('ploidy', metavar='PLOIDY', type=int,
		help='The ploidy of the sample(s).')
	
	arg('-o', '--output', default=sys.stdout,
		help='Output VCF file. Add .gz to the file name to get compressed output. '
			'If omitted, use standard output.')
	arg('--reference', '-r', metavar='FASTA',
		help='Reference file. Provide this to detect alleles through re-alignment. '
			'If no index (.fai) exists, it will be created')
	arg('--tag', choices=('PS','HP'), default='PS',
		help='Store phasing information with PS tag (standardized) or '
			'HP tag (used by GATK ReadBackedPhasing) (default: %(default)s)')
	arg('--output-read-list', metavar='FILE', default=None, dest='read_list_filename',
		help='Write reads that have been used for phasing to FILE.')

	arg = parser.add_argument_group('Input pre-processing, selection, and filtering').add_argument
	arg('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	arg('--indels', dest='indels', default=False, action='store_true',
		help='Also phase indels (default: do not phase indels)')
	arg('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM/CRAM header and assume all reads come '
			'from the same sample.')
	arg('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
		help='Name of a sample to phase. If not given, all samples in the '
		'input VCF are phased. Can be used multiple times.')
	arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
		help='Name of chromosome to phase. If not given, all chromosomes in the '
		'input VCF are phased. Can be used multiple times.')

	arg = parser.add_argument_group('Parameters for cluster editing').add_argument
	arg('--errorrate', metavar='ERROR', type=float, default=0.1, help='Read error rate (default: %(default)s).')
	arg('--min-overlap', metavar='OVERLAP', type=int, default=5, help='Minimum required read overlap (default: %(default)s).')

def validate(args, parser):
	pass

def main(args):
	run_clusterediting(**vars(args))