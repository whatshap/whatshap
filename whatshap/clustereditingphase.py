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
from .core import Read, ReadSet, CoreAlgorithm, LightCompleteGraph, readselection, NumericSampleIds, GenotypeLikelihoods, Genotype, compute_genotypes
from .graph import ComponentFinder
from .bam import AlignmentFileNotIndexedError, SampleNotFoundError, ReferenceNotFoundError, EmptyAlignmentFileError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError
from .utils import detect_file_format, IndexedFasta, FastaNotIndexedError
from .matrixtransformation import MatrixTransformation
from .phase import read_reads, select_reads, split_input_file_list, setup_pedigree, find_components, find_largest_component, write_read_list
from .clustereditingplots import draw_plots_dissimilarity, draw_plots_scoring, draw_column_dissimilarity, draw_heatmaps, draw_superheatmap, draw_cluster_coverage, draw_cluster_blocks
from .readscoring import score, locality_sensitive_score
from .kclustifier import clusters_to_haps, clusters_to_blocks, avg_readlength, calc_consensus_blocks, subset_clusters
#from .core import clusters_to_haps, clusters_to_blocks, avg_readlength, calc_consensus_blocks, subset_clusters
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

def run_clustereditingphase(
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
	min_overlap = 5,
	transform = False,
	plot_heatmap = False,
	plot_haploblocks = False
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
		output_str = output
		if isinstance(output, str):
			output = stack.enter_context(xopen(output, 'w'))
		if write_command_line_header:
			command_line = '(whatshap {}) {}'.format(__version__, ' '.join(sys.argv[1:]))
		else:
			command_line = None
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=variant_file,
			out_file=output, tag=tag, ploidy=ploidy)
		# TODO for now, assume we always trust the genotypes
		vcf_reader = VcfReader(variant_file, indels=indels, phases=True, genotype_likelihoods=False, ploidy=ploidy)

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
				readset = readset.subset([i for i, read in enumerate(readset) if len(read) >= max(2,min_overlap)])
				logger.info('Kept %d reads that cover at least two variants each', len(readset))

				# sample allele matrix
				selected_reads = select_reads(readset, 5*ploidy, preferred_source_ids = vcf_source_ids)
				readset = selected_reads

				# Transform allele matrix, if option selected
				timers.start('transform_matrix')
				if transform:
					logger.info("Transforming allele matrix...")
					transformation = MatrixTransformation(readset, find_components(readset.get_positions(), readset), ploidy, errorrate, min_overlap)
					readset = transformation.get_transformed_matrix()
					cluster_counts = transformation.get_cluster_counts()
				timers.stop('transform_matrix')

				# Compute similarity values for all read pairs
				logger.info("Computing similarities for read pairs ...")
				similarities = score(readset, ploidy, errorrate, min_overlap)
				#similarities = locality_sensitive_score(readset, ploidy, min_overlap)
				
				# Create read graph object
				logger.info("Constructing graph ...")
				graph = LightCompleteGraph(len(readset),True)

				# Insert edges into read graph
				n_reads = len(readset)
				for id1 in range(n_reads):
					for id2 in range(id1+1, n_reads):
						if similarities.get(id1, id2) != 0:
							graph.setWeight(id1, id2, similarities.get(id1, id2))

				# Run cluster editing
				logger.info("Solving cluster editing ...")
				timers.start('solve_clusterediting')
				clusterediting = CoreAlgorithm(graph)	
				readpartitioning = clusterediting.run()
				timers.stop('solve_clusterediting')

				# Assemble clusters to haplotypes
				logger.info("Assembling haplotypes from read clusters ...")
				timers.start('assemble_haplotypes')

				#consensus_blocks = clusters_to_haps(readset, readpartitioning, ploidy, coverage_padding = 7, copynumber_max_artifact_len = 0.5, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = True)
			#	coverage, copynumbers, cluster_blocks, cut_positions = clusters_to_blocks(readset, readpartitioning, ploidy, coverage_padding = 7, copynumber_max_artifact_len = 0.5, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = True)

####################################### new ##################################
				#add dynamic programming for finding the most likely subset of clusters
				coverage, cut_positions, cluster_blocks = subset_clusters(readset, readpartitioning, ploidy)
##############################################################################

				#truth = get_phase(readset, phasable_variant_table)
				#truth_cons = construct_ground_truth(readset)

				#truth_blocks = construct_ground_truth_blockwise(readset, cut_positions)
				#for hap in truth:
				#	print("".join(list(map(lambda x: str(x) if x >= 0 else ".", hap))))
				#for hap in truth_cons:
				#	print("".join(list(map(lambda x: str(x) if x >= 0 else ".", hap))))
				#vec_error = vector_error_blockwise(consensus_blocks, truth, 1, 100000)
				#print("Vector error = "+str(vec_error))

				#vec_error = vector_error_blockwise(consensus_blocks, truth_cons, 1, 100000)
				#print("Vector error = "+str(vec_error))
				timers.stop('assemble_haplotypes')

				logger.info("Generating plots ...")
				if plot_heatmap:
					draw_superheatmap(readset, readpartitioning, phasable_variant_table, output_str+ ("_D" if transform else "_N") + ".superheatmapg.png", genome_space = False)
				if plot_haploblocks:
					draw_cluster_blocks(readset, readpartitioning, cluster_blocks, cut_positions, phasable_variant_table, output_str+ ("_D" if transform else "_N") + ".haploblocks.png", genome_space = False)
				
				#draw_cluster_coverage(readset, readpartitioning, output_str+ ("_D" if transform else "_N") + ".coverage.png")
				
				print("... finished")
				
			with timers('write_vcf'):
				logger.info('======== Writing VCF')
				
				# !!!
				# TODO: Write results into VCF: Create superreads and components for VCF-Write call below
				#changed_genotypes = vcf_writer.write(chromosome, superreads, components)
				#assert len(changed_genotypes) == 0
				# !!!

				logger.info('Done writing VCF')
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
	logger.info('Time spent transforming allele matrix:       %6.1f s', timers.elapsed('transform_matrix'))
	logger.info('Time spent solving cluster editing:          %6.1f s', timers.elapsed('solve_clusterediting'))
	logger.info('Time spent assembling haplotypes:            %6.1f s', timers.elapsed('assemble_haplotypes'))
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
	arg('--transform', dest='transform', default=False, action='store_true',
		help='Use transformed matrix for read similarity scoring (default: %(default)s).')
	arg('--plot-heatmap', dest='plot_heatmap', default=False, action='store_true',
		help='Plot a super heatmap for the computed clustering (default: %(default)s).')
	arg('--plot-haploblocks', dest='plot_haploblocks', default=False, action='store_true',
		help='Plot the haplotype blocks with contained reads (default: %(default)s).')

def validate(args, parser):
	pass

def main(args):
	run_clustereditingphase(**vars(args))
