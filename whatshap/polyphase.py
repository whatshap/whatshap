"""
Cluster all reads based on their pairwise similarity using cluster editing and assemble multiple haplotypes from the clusters.

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
Each read is represented as a node in a graph with their pairwise score (based
on similarity) as edge weights. Reads originating from the same haplotype with
high confidence are clusterd togehter using cluster editing. The haplotype
sequences are computed variant-wise by choosing a set of paths through the read
clusters and by respecting coverage and genotype constraints.

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
from .core import Read, ReadSet, CoreAlgorithm, DynamicSparseGraph, readselection, NumericSampleIds, GenotypeLikelihoods, Genotype, compute_genotypes
from .graph import ComponentFinder
from .bam import AlignmentFileNotIndexedError, SampleNotFoundError, ReferenceNotFoundError, EmptyAlignmentFileError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError
from .utils import detect_file_format, IndexedFasta, FastaNotIndexedError
from .matrixtransformation import MatrixTransformation
from .phase import read_reads, select_reads, split_input_file_list, setup_pedigree, find_components, find_largest_component, write_read_list
from .clustereditingplots import draw_plots_dissimilarity, draw_plots_scoring, draw_column_dissimilarity, draw_heatmaps, draw_superheatmap, draw_cluster_coverage, draw_cluster_blocks, draw_dp_threading
from .readscoring import score_global, score_local, score_local_patternbased
from .kclustifier import clusters_to_haps, clusters_to_blocks, avg_readlength, calc_consensus_blocks
from .threading import subset_clusters, get_cluster_consensus_local
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

def run_polyphase(
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
	ce_bundle_edges = False,
	ce_score_global = False,
	ce_score_with_patterns = False,
	min_overlap = 2,
	transform = False,
	plot_clusters = False,
	plot_haploblocks = False,
	plot_threading = False,
	single_block = False,
	cpp_threading = False,
	ce_refine = False
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
				timers.start('read_bam')
				bam_sample = None if ignore_read_groups else sample
				readset, vcf_source_ids = read_reads(readset_reader, chromosome, phasable_variant_table.variants, bam_sample, fasta, [], numeric_sample_ids, phase_input_bam_filenames)
				readset.sort()
				# TODO: len == min_overlap ?
				readset = readset.subset([i for i, read in enumerate(readset) if len(read) >= max(2,min_overlap)])
				logger.info('Kept %d reads that cover at least two variants each', len(readset))

				#adapt the variant table to the subset of reads
				variant_table.subset_rows_by_position(readset.get_positions())
				
				#compute the genotypes that belong to the variant table and create a list of all genotypes				
				all_genotypes = variant_table.genotypes_of(sample)
				genotype_list = []
				for pos in range(len(all_genotypes)):
					gen = 0
					for allele in all_genotypes[pos].get_genotype().as_vector():
						gen += allele
					genotype_list.append(gen)

				# sample allele matrix
				#selected_reads = select_reads(readset, 5*ploidy, preferred_source_ids = vcf_source_ids)
				#readset = selected_reads
				timers.stop('read_bam')

				# Transform allele matrix, if option selected
				timers.start('transform_matrix')
				if transform:
					logger.info("Transforming allele matrix...")
					transformation = MatrixTransformation(readset, find_components(readset.get_positions(), readset), ploidy, min_overlap)
					readset = transformation.get_transformed_matrix()
					cluster_counts = transformation.get_cluster_counts()
				timers.stop('transform_matrix')

				# Compute similarity values for all read pairs
				timers.start('compute_graph')
				logger.info("Computing similarities for read pairs ...")
				if ce_score_with_patterns:
					similarities = score_local_patternbased(readset, ploidy, errorrate, min_overlap, 4)
				elif ce_score_global:
					similarities = score_global(readset, ploidy, min_overlap)
				else:
					similarities = score_local(readset, ploidy, min_overlap)
				
				# Create read graph object
				logger.info("Constructing graph ...")
				graph = DynamicSparseGraph()

				# Insert edges into read graph
				for (read1, read2) in similarities:
					graph.addEdge(read1, read2, similarities.get(read1, read2))
				timers.stop('compute_graph')

				# Run cluster editing
				logger.info("Solving cluster editing instance with {} nodes and {} edges ...".format(len(readset), len(similarities)))
				timers.start('solve_clusterediting')
				clusterediting = CoreAlgorithm(graph, ce_bundle_edges)
				clustering = clusterediting.run()
				
				if ce_refine:
					runs_remaining = 3
					refine = True
					while refine and runs_remaining > 0:
						refine = False
						runs_remaining -= 1
						#for every cluster and every variant position, compute the relative coverage and consensus
						#logger.info("Detecting single variant inconsistencies in clusters ...")

						# Map genome positions to [0,l)
						num_clusters = len(clustering)
						index = {}
						rev_index = []
						num_vars = 0

						for position in readset.get_positions():
							index[position] = num_vars
							rev_index.append(position)
							num_vars += 1

						#cov_map maps each position to the list of cluster IDS that appear at this position
						cov_positions = defaultdict()
						cov_map = defaultdict(list)
						for c_id in range(num_clusters):
							covered_positions = []
							for read in clustering[c_id]:
								for pos in [index[var.position] for var in readset[read] ]:		
									if pos not in covered_positions:
										covered_positions.append(pos)
							for p in covered_positions:
								cov_map[p].append(c_id)
						#	assert(len(cov_map.keys()) == num_vars)

						#compute for each position the amount of clusters that 'cover' this position
						for key in cov_map.keys():
							cov_positions[key] = len(cov_map[key])
						cov_positions_sorted = sorted(cov_positions.items(), key=lambda x: x[1], reverse=True)

						#restrict the number of clusters at each position to a maximum of 2*ploidy clusters with the largest number of reads	
						for key in cov_map.keys():
							largest_clusters = sorted(cov_map[key], key=lambda x: len(clustering[x]), reverse=True)[:min(len(cov_map[key]), 2*ploidy)]
							cov_map[key] = largest_clusters

						#create dictionary mapping a clusterID to a pair (starting position, ending position)
						positions = {}
						counter = 0
						for c_id in range(num_clusters):
							read = clustering[c_id][0]
							start = index[readset[read][0].position]
							end = index[readset[read][-1].position]
							for read in clustering[c_id]:
								readstart = index[readset[read][0].position]
								readend = index[readset[read][-1].position]
								if (readstart < start):
									start = readstart
								if (readend > end):
									end = readend
							positions[c_id] = (start,end)
						assert(len(positions) == num_clusters)

						coverage = [dict() for pos in range(num_vars)]
						coverage_sum = [0 for pos in range(num_vars)]
						for c_id in range(num_clusters):
							for read in clustering[c_id]:
								for pos in [index[var.position] for var in readset[read]]:
									if c_id not in coverage[pos]:
										coverage[pos][c_id] = 0
									coverage[pos][c_id] += 1
									coverage_sum[pos] += 1

						for pos in range(num_vars):
							for c_id in coverage[pos]:
								coverage[pos][c_id] = coverage[pos][c_id]/coverage_sum[pos]

						consensus = get_cluster_consensus_local(readset, clustering, cov_map, positions, fractional = True)

						inconsistency_count = 0
						for pos in range(num_vars):
							for c_id in coverage[pos]:
								if coverage[pos][c_id] > 0.333 and 0.30 <= consensus[pos][c_id] <= 0.70:
									#print("   inconsistency in cluster "+str(c_id)+" at position"+str(pos))
									refine = True
									inconsistency_count += 1
									zero_reads = []
									one_reads = []
									for read in clustering[c_id]:
										for var in readset[read]:
											if index[var.position] == pos:
												if var.allele == 0:
													zero_reads.append(read)
												else:
													one_reads.append(read)
									for r0 in zero_reads:
										for r1 in one_reads:
											similarities.set(r0, r1, -10000)

						refined_graph = DynamicSparseGraph()
						for (read1, read2) in similarities:
							refined_graph.addEdge(read1, read2, similarities.get(read1, read2))

						if refine:
							logger.info(str(inconsistency_count)+" inconsistent variants found. Refining clusters ...")
							refined_clusterediting = CoreAlgorithm(refined_graph, ce_bundle_edges)
							clustering = refined_clusterediting.run()
				
				timers.stop('solve_clusterediting')

				# Assemble clusters to haplotypes
				logger.info("Assembling {} haplotypes from {} clusters over {} positions ...".format(ploidy, len(clustering), len(readset.get_positions())))
				timers.start('assemble_haplotypes')

				# Add dynamic programming for finding the most likely subset of clusters
				cut_positions, cluster_blocks, components, superreads, coverage, paths, haplotypes = subset_clusters(readset, clustering, ploidy, sample,genotype_list, single_block, cpp_threading)

				timers.stop('assemble_haplotypes')

				# Plot options
				timers.start('create_plots')
				if plot_clusters or plot_haploblocks or plot_threading:
					logger.info("Generating plots ...")
					if plot_clusters:
						if ce_score_with_patterns:
							draw_plots_scoring(readset, similarities, output_str+".scoringplot.pdf", ploidy, errorrate, min_overlap)
						draw_superheatmap(readset, clustering, phasable_variant_table, output_str+".clusters.pdf", genome_space = False)
					if plot_haploblocks:
						draw_cluster_blocks(readset, clustering, cluster_blocks, cut_positions, phasable_variant_table, output_str+".haploblocks.pdf", genome_space = False)
					if plot_threading:
						draw_dp_threading(coverage, paths, cut_positions, haplotypes, readset, phasable_variant_table, output_str+".threading.pdf")
				timers.stop('create_plots')

			with timers('write_vcf'):
				logger.info('======== Writing VCF')
				changed_genotypes = vcf_writer.write(chromosome, superreads, components)
				# TODO: Use genotype information to polish results
				#assert len(changed_genotypes) == 0
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
	logger.info('Time spent computing read graph:             %6.1f s', timers.elapsed('compute_graph'))
	logger.info('Time spent solving cluster editing:          %6.1f s', timers.elapsed('solve_clusterediting'))
	logger.info('Time spent assembling haplotypes:            %6.1f s', timers.elapsed('assemble_haplotypes'))
	if plot_clusters or plot_haploblocks or plot_threading:
		logger.info('Time spent creating plots:                   %6.1f s', timers.elapsed('create_plots'))
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
	arg('--ce-score-global', dest='ce_score_global', default=False, action='store_true',
		help='Reads are scored with respect to their location inside the chromosome. (default: %(default)s).')
	arg('--ce-score-with-patterns', dest='ce_score_with_patterns', default=False, action='store_true',
		help='Uses a scoring method for reads, which is based on local haplotype inference through local patterns of the reads (default: %(default)s).')
	arg('--ce-bundle-edges', dest='ce_bundle_edges', default=False, action='store_true',
		help='Influences the cluster editing heuristic. Only for debug/developing purpose (default: %(default)s).')
	arg('--min-overlap', metavar='OVERLAP', type=int, default=2, help='Minimum required read overlap (default: %(default)s).')
	arg('--transform', dest='transform', default=False, action='store_true',
		help='Use transformed matrix for read similarity scoring (default: %(default)s).')
	arg('--plot-clusters', dest='plot_clusters', default=False, action='store_true',
		help='Plot a super heatmap for the computed clustering (default: %(default)s).')
	arg('--plot-haploblocks', dest='plot_haploblocks', default=False, action='store_true',
		help='Plot the haplotype blocks with contained reads (default: %(default)s).')
	arg('--plot-threading', dest='plot_threading', default=False, action='store_true',
		help='Plot the haplotypes\' threading through the read clusters (default: %(default)s).')
	arg('--single-block', dest='single_block', default=False, action='store_true',
		help='Output only one single block.')
	arg('--cpp-threading', dest='cpp_threading', default=False, action='store_true',
		help='Uses a C++ implementation to perform the haplotype threading.')
	arg('--ce-refine', dest='ce_refine', default=False, action='store_true',
		help='Refines the output of cluster editing by detecting inconsistencies within the clusters and rerunning the clustering.')

def validate(args, parser):
	pass

def main(args):
	run_polyphase(**vars(args))

