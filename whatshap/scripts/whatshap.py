#!/usr/bin/env python3
"""
Read a VCF and a BAM file and phase the variants. The phased VCF is written to
standard output.
"""
"""
 0: ref allele
 1: alt allele
 -: unphasable: no coverage of read that covers at least 2 SNPs
 X: unphasable: there is coverage, but still not phasable (tie)

TODO
* Perhaps simplify slice_reads() such that it only creates and returns one slice
* it would be cleaner to not open the input VCF twice
* convert parse_vcf to a class so that we can access VCF header info before
  starting to iterate (sample names)
"""
import os
import logging
import sys
import random
import gzip
import time
import itertools
from collections import defaultdict
from contextlib import ExitStack, closing
import pysam

from ..vcf import parse_vcf, PhasedVcfWriter
from .. import __version__
from ..args import HelpfulArgumentParser as ArgumentParser
from ..core import Read, ReadSet, DPTable, IndexSet
from ..graph import ComponentFinder


__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)


def find_alleles(variants, start, bam_read, core_read):
	"""
	Check whether the bam_read covers some of the variants in variants[start:] and, if so,
	add this information to the given core_read object.
	"""
	pos = bam_read.pos
	cigar = bam_read.cigar

	j = start  # index into variants list
	p = pos
	s = 0  # absolute index into the read string [0..len(read)]
	errors = 0
	for cigar_op, length in cigar:
		# The mapping of CIGAR operators to numbers is:
		# MIDNSHPX= => 012345678
		if cigar_op in (0, 7, 8):  # we're in a matching subregion
			s_next = s + length
			p_next = p + length
			r = p + length  # size of this subregion
			# skip over all SNPs that come before this region
			while j < len(variants) and variants[j].position < p:
				j += 1
			# iterate over all positions in this subregion and
			# check whether any of them coincide with one of the SNPs ('hit')
			while j < len(variants) and p < r:
				if variants[j].position == p:  # we have a hit
					base = bam_read.seq[s:s+1]
					al = None
					if base == variants[j].reference_allele:
						al = 0  # REF allele
					elif base == variants[j].alternative_allele:
						al = 1  # ALT allele
					else:
						errors += 1
					if al is not None:
						# Just ignore duplicate variants encountered when two reads in a pair overlap
						# TODO: Handle this case properly
						if not p in core_read:
							core_read.add_variant(p, base, al, ord(bam_read.qual[s:s+1])-33)
					j += 1
				s += 1 # advance both read and reference
				p += 1
			s = s_next
			p = p_next
		elif cigar_op == 1:  # an insertion
			s += length
		elif cigar_op == 2 or cigar_op == 3:  # a deletion or a reference skip
			p += length
		elif cigar_op == 4:  # soft clipping
			s += length
		elif cigar_op == 5 or cigar_op == 6:  # hard clipping or padding
			pass
		else:
			logger.error("Unsupported CIGAR operation: %d", cigar_op)
			sys.exit(1)
	return errors


class BamReader:
	"""
	Associate variants with reads.
	"""
	def __init__(self, path, mapq_threshold=20):
		"""
		path -- path to BAM file
		"""
		bai1 = path + '.bai'
		bai2 = os.path.splitext(path)[0] + '.bai'
		if not os.path.exists(bai1) and not os.path.exists(bai2):
			logger.info('BAM index not found, creating it now.')
			pysam.index(path)
		self._samfile = pysam.Samfile(path)
		self._mapq_threshold = mapq_threshold
		self._initialize_sample_to_group_ids()

	def _initialize_sample_to_group_ids(self):
		"""
		Return a dictionary that maps a sample name to a set of read group ids.
		"""
		read_groups = self._samfile.header['RG']  # a list of dicts
		logger.debug('Read groups in SAM header: %s', read_groups)
		samples = defaultdict(list)
		for read_group in read_groups:
			samples[read_group['SM']].append(read_group['ID'])
		self._sample_to_group_ids = {
			id: frozenset(values) for id, values in samples.items() }

	def read(self, chromosome, variants, sample):
		"""
		chromosome -- name of chromosome to work on
		variants -- list of Variant objects (obtained from VCF with parse_vcf)
		sample -- name of sample to work on. If None, read group information is
			ignored and all reads in the file are used.

		Return a ReadSet object.
		"""
		if sample is not None:
			read_groups = self._sample_to_group_ids[sample]

		# resulting set of reads
		result = ReadSet()

		i = 0  # keep track of position in variants array (which is in order)
		for bam_read in self._samfile.fetch(chromosome):
			if sample is not None and not bam_read.opt('RG') in read_groups:
				continue
			# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
			if bam_read.flag & 2048 != 0:
				# print('Skipping additional alignment for read ', bam_read.qname)
				continue
			if bam_read.is_secondary:
				continue
			if bam_read.is_unmapped:
				continue
			if bam_read.mapq < self._mapq_threshold:
				continue
			if not bam_read.cigar:
				continue

			# since reads are ordered by position, we do not need to consider
			# positions that are too small
			while i < len(variants) and variants[i].position < bam_read.pos:
				i += 1
			try:
				core_read = result[bam_read.qname]
				former_length = len(core_read)
				find_alleles(variants, i, bam_read, core_read)
				# If variants on the current (part of the) read have been added,
				# then also record its MAPQ
				if len(core_read) > former_length:
					core_read.addMapq(bam_read.mapq)
			except KeyError:
				core_read = Read(bam_read.qname, bam_read.mapq)
				find_alleles(variants, i, bam_read, core_read)
				# only add new read if it contained at least one variant
				if len(core_read) > 0:
					result.add(core_read)
		return result

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()

	def close(self):
		self._samfile.close()


class CoverageMonitor:
	'''TODO: This is a most simple, naive implementation. Could do this smarter.'''
	def __init__(self, length):
		self.coverage = [0] * length

	def max_coverage_in_range(self, begin, end):
		return max(self.coverage[begin:end])

	def add_read(self, begin, end):
		for i in range(begin, end):
			self.coverage[i] += 1


def slice_reads(reads, max_coverage):
	"""
	Iterate over all read in random order and greedily retain those reads whose
	addition does not lead to a local physical coverage exceeding the given threshold.
	Return a ReadSet containing the retained reads.

	max_coverage -- Slicing ensures that the (physical) coverage does not exceed max_coverage anywhere along the chromosome.
	reads -- a ReadSet
	"""
	shuffled_indices = list(range(len(reads)))
	random.shuffle(shuffled_indices)

	position_list = reads.get_positions()
	logger.info('Found %d SNP positions', len(position_list))

	# dictionary to map SNP position to its index
	position_to_index = { position: index for index, position in enumerate(position_list) }

	# List of slices, start with one empty slice ...
	slices = [IndexSet()]
	# ... and the corresponding coverages along each slice
	slice_coverages = [CoverageMonitor(len(position_list))]
	skipped_reads = 0
	accessible_positions = set()
	for index in shuffled_indices:
		read = reads[index]
		# Skip reads that cover only one SNP
		if len(read) < 2:
			skipped_reads += 1
			continue
		for position, base, allele, quality in read:
			accessible_positions.add(position)
		first_position, first_base, first_allele, first_quality = read[0]
		last_position, last_base, last_allele, last_quality = read[len(read)-1]
		begin = position_to_index[first_position]
		end = position_to_index[last_position] + 1
		slice_id = 0
		while True:
			# Does current read fit into this slice?
			if slice_coverages[slice_id].max_coverage_in_range(begin, end) < max_coverage:
				slice_coverages[slice_id].add_read(begin, end)
				slices[slice_id].add(index)
				break
			else:
				slice_id += 1
				# do we have to create a new slice?
				if slice_id == len(slices):
					slices.append(IndexSet())
					slice_coverages.append(CoverageMonitor(len(position_list)))
	logger.info('Skipped %d reads that only cover one SNP', skipped_reads)

	unphasable_snps = len(position_list) - len(accessible_positions)
	if position_list:
		logger.info('%d out of %d variant positions (%.1d%%) do not have a read '
			'connecting them to another variant and are thus unphasable',
			unphasable_snps, len(position_list),
			100. * unphasable_snps / len(position_list))

	# Print stats
	for slice_id, index_set in enumerate(slices):
		logger.info('Slice %d contains %d reads', slice_id, len(index_set))

	return reads.subset(slices[0])


def find_components(superreads, reads):
	"""
	Return a dict that maps each position to the component it is in. A
	component is identified by the position of its leftmost variant.
	"""
	logger.info('Finding connected components ...')
	assert len(superreads) == 2
	assert len(superreads[0]) == len(superreads[1])

	phased_positions = [ position for position, base, allele, quality in superreads[0] if allele in [0, 1] ]  # TODO set()
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [ position for position, base, allele, quality in read if position in phased_positions ]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	components = { position : component_finder.find(position) for position in phased_positions }
	logger.info('No. of variants considered for phasing: %d', len(superreads[0]))
	logger.info('No. of variants that were phased: %d', len(phased_positions))
	return components


def best_case_blocks(reads):
	"""
	Given a list of core reads, determine the number of phased blocks that
	would result if each variant were actually phased.

	Return the number of connected components.
	"""
	positions = set()
	for read in reads:
		for position, _, _, _ in read:
			positions.add(position)
	component_finder = ComponentFinder(positions)
	for read in reads:
		read_positions = [ position for position, _, _, _ in read ]
		for position in read_positions[1:]:
			component_finder.merge(read_positions[0], position)
	# A dict that maps each position to the component it is in.
	components = { component_finder.find(position) for position in positions }
	return len(components)


def ensure_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
		sys.exit("WhatsHap requires pysam >= 0.8.1")


def main():
	ensure_pysam_version()
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = ArgumentParser(prog='whatshap', description=__doc__)
	parser.add_argument('--version', action='version', version=__version__)
	parser.add_argument('-o', '--output', default=None,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality')
	parser.add_argument('--seed', default=123, type=int, help='Random seed (default: %(default)s)')
	parser.add_argument('--all-het', action='store_true', default=False,
		help='Assume all positions to be heterozygous (that is, fully trust SNP calls).')
	parser.add_argument('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM header and assume all reads come '
		'from the same sample.')
	parser.add_argument('--sample', metavar='SAMPLE', default=None,
		help='Name of a sample to phase. If not given, only the first sample '
			'in the input VCF is phased.')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('bam', metavar='BAM', help='BAM file')
	args = parser.parse_args()
	random.seed(args.seed)

	start_time = time.time()
	with ExitStack() as stack:
		try:
			bam_reader = stack.enter_context(closing(BamReader(args.bam, mapq_threshold=args.mapping_quality)))
		except OSError as e:
			logging.error(e)
			sys.exit(1)
		if args.output is not None:
			out_file = stack.enter_context(open(args.output, 'w'))
		else:
			out_file = sys.stdout
		command_line = ' '.join(sys.argv[1:])
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=args.vcf, out_file=out_file)
		for sample, chromosome, variants in parse_vcf(args.vcf, args.sample):
			logger.info('Read %d variants on chromosome %s', len(variants), chromosome)
			if args.ignore_read_groups:
				sample = None
			logger.info('Reading the BAM file ...')

			reads = bam_reader.read(chromosome, variants, sample)
			
			# Sort the variants stored in each read
			# TODO: Check whether this is already ensured by construction
			for read in reads:
				read.sort()
			# Sort reads in read set by position
			reads.sort()

			sliced_reads = slice_reads(reads, args.max_coverage)
			logger.info('Best-case phasing would result in %d phased blocks (%d with slicing)',
				best_case_blocks(reads), best_case_blocks(sliced_reads))
			logger.info('Phasing the variants (using %d reads)...', len(sliced_reads))

			# Run the core algorithm: construct DP table ...
			dp_table = DPTable(sliced_reads, args.all_het)
			# ... and do the backtrace to get the solution
			superreads = dp_table.get_super_reads()

			components = find_components(superreads, sliced_reads)
			logger.info('No. of phased blocks: %d', len(set(components.values())))
			logger.info('Writing phased variants on chromosome %s ...', chromosome)
			vcf_writer.write(chromosome, sample, superreads, components)

	logger.info('Elapsed time: %.1fs', time.time() - start_time)
