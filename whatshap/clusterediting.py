"""
Run polyploid phasing based on cluster editing.

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.

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
from .core import Read, ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, GenotypeLikelihoods, Genotype, compute_genotypes, LightCompleteGraph, CoreAlgorithm
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import AlignmentFileNotIndexedError, SampleNotFoundError, ReferenceNotFoundError, EmptyAlignmentFileError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError
from .utils import detect_file_format, IndexedFasta, FastaNotIndexedError
#from .readsetpruning import ReadSetPruning
from .readscoring import score
from .phase import read_reads, select_reads, split_input_file_list, setup_pedigree, find_components, find_largest_component, write_read_list

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
#	read_list_filename=None
	errorrate = 0.1,
	min_overlap = 3
	):
	"""
	Run Polyploid Phasing based on cluster editing.
	
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
	reads_per_window -- Maximum number of reads to consider in an window
	variants_per_window -- Minimum number of variants that need to be covered by all reads in a window
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
#		if read_list_filename:
#			read_list_file = create_read_list_file(read_list_filename)
		
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
				readset = readset.subset([i for i, read in enumerate(readset) if len(read) >= 2])

				# TODO: Keep readselection step?
				selected_reads = select_reads(readset, 2*ploidy, preferred_source_ids = vcf_source_ids)
				readset = selected_reads

				print('readset: length: ' + str(len(readset)))
				print_readset(readset)
				logger.info('Kept %d reads that cover at least two variants each', len(readset))

				# create read graph object
				graph = LightCompleteGraph(len(readset),True)
				# compute pairwise read similarities
				similarities = score(readset, ploidy, errorrate, min_overlap)
				# compute pairwise read similarities
				print(similarities)
				# insert edges into read graph
				n_reads = len(readset)
				for id1 in range(n_reads):
					for id2 in range(id1+1, n_reads):
						graph.setWeight(id1, id2, similarities[id1][id2 - id1 - 1])
				# TODO: keep track of mapping of read ids to readnames
				# run cluster editing
				clusterediting = CoreAlgorithm(graph)	
				readpartitioning = clusterediting.run()
				print(readpartitioning)

				# TODO: how to proceed??

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
#	arg('--tag', choices=('PS','HP'), default='PS',
#		help='Store phasing information with PS tag (standardized) or '
#			'HP tag (used by GATK ReadBackedPhasing) (default: %(default)s)')
#	arg('--output-read-list', metavar='FILE', default=None, dest='read_list_filename',
#		help='Write reads that have been used for phasing to FILE.')

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

	arg = parser.add_argument_group('Cluster editing').add_argument
	arg('--errorrate', metavar='ERROR', default=0.1, type=float,
		help='Read error rate (default: %(default)s)')
	arg('--min-overlap', metavar='OVERLAP', default=3, type=int,
		help='Minimum read overlap (default: %(default)s)')

def validate(args, parser):
	pass

def main(args):
	run_clusterediting(**vars(args))
