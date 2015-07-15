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
import platform
from collections import defaultdict, Counter
from contextlib import contextmanager
try:
	from contextlib import ExitStack
except ImportError:
	from contextlib2 import ExitStack  # PY32
from ..vcf import parse_vcf, PhasedVcfWriter, remove_overlapping_variants
from .. import __version__
from ..args import HelpfulArgumentParser as ArgumentParser
from ..core import Read, ReadSet, DPTable, readselection
from ..graph import ComponentFinder
from ..coverage import CovMonitor
from ..bam import MultiBamReader, SampleBamReader, BamIndexingError, SampleNotFoundError, HaplotypeBamWriter


__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)


def covered_variants(variants, start, bam_read, source_id):
	"""
	Find the variants that are covered by the given bam_read and return a
	core.Read instance that represents those variants. The instance may be
	empty.

	start -- index of the first variant (in the variants list) to check
	"""
	core_read = Read(bam_read.qname, bam_read.mapq, source_id)
	i = 0  # index into CIGAR
	j = start  # index into variants
	ref_pos = bam_read.pos  # position relative to reference
	query_pos = 0  # position relative to read

	for cigar_op, length in bam_read.cigar:
		# The mapping of CIGAR operators to numbers is:
		# MIDNSHPX= => 012345678
		if cigar_op in (0, 7, 8):  # we are in a matching region
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1

			# Iterate over all variants that are in this region
			while j < len(variants) and variants[j].position < ref_pos + length:
				if len(variants[j].reference_allele) == len(variants[j].alternative_allele) == 1:
					# Variant is a SNP
					offset = variants[j].position - ref_pos
					base = bam_read.seq[query_pos + offset]
					allele = None
					if base == variants[j].reference_allele:
						allele = 0
					elif base == variants[j].alternative_allele:
						allele = 1
					if allele is not None:
						# TODO
						# Fix this: we can actually have indel and SNP
						# calls at identical positions. For now, ignore the
						# second variant.
						if variants[j].position in core_read:
							logger.debug("Found two variants at identical positions. Ignoring the second one: %s", variants[j])
						else:
							# Do not use bam_read.qual here as it is extremely slow.
							# If we ever decide to be compatible with older pysam
							# versions, cache bam_read.qual somewhere - do not
							# access it within this loop (3x slower otherwise).
							core_read.add_variant(variants[j].position, allele, bam_read.query_qualities[query_pos + offset])
				elif len(variants[j].reference_allele) == 0:
					assert len(variants[j].alternative_allele) > 0
					# This variant is an insertion. Since we are in a region of
					# matches, the insertion was *not* observed (reference allele).
					qual = 30  #  TODO average qualities of "not inserted" bases?
					core_read.add_variant(variants[j].position, allele=0, quality=qual)
				elif len(variants[j].alternative_allele) == 0:
					assert len(variants[j].reference_allele) > 0
					# This variant is a deletion that was not observed.
					# Add it only if the next variant is not located 'within'
					# the deletion.
					deletion_end = variants[j].position + len(variants[j].reference_allele)
					if not (j + 1 < len(variants) and variants[j+1].position < deletion_end):
						qual = 30  # TODO
						core_read.add_variant(variants[j].position, allele=0, quality=qual)
					else:
						logger.info('Skipped a deletion overlapping another variant at pos. %d', variants[j].position)
						# Also skip all variants that this deletion overlaps
						while j + 1 < len(variants) and variants[j+1].position < deletion_end:
							j += 1
						# One additional j += 1 is done below
				else:
					assert False, "Strange variant: {}".format(variants[j])
				j += 1
			query_pos += length
			ref_pos += length
		elif cigar_op == 1:  # an insertion
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1
			if j < len(variants) and variants[j].position == ref_pos and \
					len(variants[j].reference_allele) == 0 and \
					variants[j].alternative_allele == bam_read.seq[query_pos:query_pos + length]:
				qual = 30  # TODO
				assert variants[j].position not in core_read
				core_read.add_variant(variants[j].position, allele=1, quality=qual)
				j += 1
			query_pos += length
		elif cigar_op == 2:  # a deletion
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1
			# We only check the length of the deletion, not the sequence
			# that gets deleted since we don’t have the reference available.
			# (We could parse the MD tag if it exists.)
			if j < len(variants) and variants[j].position == ref_pos and \
					len(variants[j].alternative_allele) == 0 and \
					len(variants[j].reference_allele) == length:
				qual = 30  # TODO
				deletion_end = variants[j].position + len(variants[j].reference_allele)
				if not (j + 1 < len(variants) and variants[j+1].position < deletion_end):
					qual = 30  # TODO
					assert variants[j].position not in core_read
					core_read.add_variant(variants[j].position, allele=1, quality=qual)
				else:
					logger.info('Skipped a deletion overlapping another variant at pos. %d', variants[j].position)
					# Also skip all variants that this deletion overlaps
					while j + 1 < len(variants) and variants[j+1].position < deletion_end:
						j += 1
					# One additional j += 1 is done below
				j += 1
			ref_pos += length
		elif cigar_op == 3:  # a reference skip
			ref_pos += length
		elif cigar_op == 4:  # soft clipping
			query_pos += length
		elif cigar_op == 5 or cigar_op == 6:  # hard clipping or padding
			pass
		else:
			logger.error("Unsupported CIGAR operation: %d", cigar_op)
			sys.exit(1)
	return core_read


class ReadSetError(Exception):
	pass


class ReadSetReader:
	"""
	Associate VCF variants with BAM reads.
	"""
	def __init__(self, paths, mapq_threshold=20):
		self._mapq_threshold = mapq_threshold
		if len(paths) == 1:
			self._reader = SampleBamReader(paths[0])
		else:
			self._reader = MultiBamReader(paths)

	def read(self, chromosome, variants, sample):
		"""
		chromosome -- name of chromosome to work on
		variants -- list of Variant objects (obtained from VCF with parse_vcf)
		sample -- name of sample to work on. If None, read group information is
			ignored and all reads in the file are used.

		Return a ReadSet object.
		"""
		# Since variants are identified by position, positions must be unique.
		if __debug__ and variants:
			varposc = Counter(variant.position for variant in variants)
			pos, count = varposc.most_common()[0]
			assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

		# Map read name to a list of Read objects. The list has two entries
		# if it is a paired-end read, one entry if the read is single-end.
		reads = defaultdict(list)

		i = 0  # keep track of position in variants array (which is in order)
		for alignment in self._reader.fetch(reference=chromosome, sample=sample):
			# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
			if alignment.bam_alignment.flag & 2048 != 0:
				# print('Skipping additional alignment for read ', alignment.bam_alignment.qname)
				continue
			if alignment.bam_alignment.mapq < self._mapq_threshold:
				continue
			if alignment.bam_alignment.is_secondary:
				continue
			if alignment.bam_alignment.is_unmapped:
				continue
			if not alignment.bam_alignment.cigar:
				continue

			# Skip variants that are to the left of this read.
			while i < len(variants) and variants[i].position < alignment.bam_alignment.pos:
				i += 1

			core_read = covered_variants(variants, i, alignment.bam_alignment, alignment.source_id)
			# Only add new read if it covers at least one variant.
			if core_read:
				reads[(alignment.source_id, alignment.bam_alignment.qname)].append(core_read)

		# Prepare resulting set of reads.
		read_set = ReadSet()

		for readlist in reads.values():
			assert len(readlist) > 0
			if len(readlist) > 2:
				raise ReadSetError("Read name {!r} occurs more than twice in the input file".format(readlist[0].name))
			if len(readlist) == 1:
				read_set.add(readlist[0])
			else:
				read_set.add(self._merge_pair(*readlist))
		return read_set

	def _merge_pair(self, read1, read2):
		"""
		Merge the two ends of a paired-end read into a single core.Read. Also
		takes care of self-overlapping read pairs.

		TODO this can be simplified as soon as a variant in a read can be
		modified.
		"""
		if read2:
			result = Read(read1.name, read1.mapqs[0], read1.source_id)
			result.add_mapq(read2.mapqs[0])
		else:
			return read1

		i1 = 0
		i2 = 0

		def add1():
			result.add_variant(read1[i1].position, read1[i1].allele, read1[i1].quality)

		def add2():
			result.add_variant(read2[i2].position, read2[i2].allele, read2[i2].quality)

		while i1 < len(read1) or i2 < len(read2):
			if i1 == len(read1):
				add2()
				i2 += 1
				continue
			if i2 == len(read2):
				add1()
				i1 += 1
				continue
			variant1 = read1[i1]
			variant2 = read2[i2]
			if variant2.position < variant1.position:
				add2()
				i2 += 1
			elif variant2.position > variant1.position:
				add1()
				i1 += 1
			else:
				# Variant on self-overlapping read pair
				assert read1[i1].position == read2[i2].position
				# If both alleles agree, merge into single variant and add up qualities
				if read1[i1].allele == read2[i2].allele:
					quality = read1[i1].quality + read2[i2].quality
					result.add_variant(read1[i1].position, read1[i1].allele, quality)
				else:
					# Otherwise, take variant with highest base quality and discard the other.
					if read1[i1].quality >= read2[i2].quality:
						add1()
					else:
						add2()
				i1 += 1
				i2 += 1
		return result

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()

	def close(self):
		self._reader.close()


def find_components(superreads, reads):
	"""
	Return a dict that maps each position to the component it is in. A
	component is identified by the position of its leftmost variant.
	"""
	logger.debug('Finding connected components ...')

	# The variant.allele attribute can be either 0 (major allele), 1 (minor allele),
	# or 3 (equal scores). If all_heterozygous is on (default), we can get
	# the combinations 0/1, 1/0 and 3/3 (the latter means: unphased).
	# If all_heterozygous is off, we can also get all other combinations.
	# In both cases, we are interested only in 0/1 and 1/0.
	phased_positions = [ v1.position for v1, v2 in zip(*superreads)
		if (v1.allele, v2.allele) in ((0, 1), (1, 0))
	]
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [ variant.position for variant in read if variant.position in phased_positions ]
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
	# A dict that maps each component to the number of SNPs it contains
	component_sizes = defaultdict(int)
	for position in positions:
		component_sizes[component_finder.find(position)] += 1
	non_singletons = [ component for component, size in component_sizes.items() if size > 1]
	return len(component_sizes), len(non_singletons)


class StageTimer:
	"""Measure run times of different stages of the program"""
	def __init__(self):
		self._start = dict()
		self._elapsed = defaultdict(float)

	def start(self, stage):
		"""Start measuring elapsed time for a stage"""
		self._start[stage] = time.time()

	def stop(self, stage):
		"""Stop measuring elapsed time for a stage."""
		t = time.time() - self._start[stage]
		self._elapsed[stage] += t
		return t

	def elapsed(self, stage):
		"""
		Return total time spent in a stage, which is the sum of the time spans
		between calls to start() and stop(). If the timer is currently running,
		its current invocation is not counted.
		"""
		return self._elapsed[stage]

	def total(self):
		"""Return sum of all times"""
		return sum(self._elapsed.values())

	@contextmanager
	def __call__(self, stage):
		self.start(stage)
		yield
		self.stop(stage)


def ensure_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
		sys.exit("WhatsHap requires pysam >= 0.8.1")


def run_whatshap(bam, vcf,
		output=sys.stdout, sample=None, ignore_read_groups=False, indels=True,
		mapping_quality=20, max_coverage=15, all_heterozygous=True, seed=123, haplotype_bams_prefix=None):
	"""
	Run WhatsHap.

	bam -- list of paths to BAM files
	vcf -- path to input VCF
	output -- path to output VCF or sys.stdout
	sample -- name of sample to phase. None means: phase first found sample.
	ignore_read_groups
	mapping_quality -- discard reads below this mapping quality
	max_coverage
	all_heterozygous
	seed -- seed for random numbers
	"""
	random.seed(seed)

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
		try:
			bam_reader = stack.enter_context(ReadSetReader(bam, mapq_threshold=mapping_quality))
		except (OSError, BamIndexingError) as e:
			logger.error(e)
			sys.exit(1)
		if output is not sys.stdout:
			output = stack.enter_context(open(output, 'w'))
		command_line = ' '.join(sys.argv[1:])
		vcf_writer = PhasedVcfWriter(command_line=command_line, in_path=vcf, out_file=output)
		vcf_reader = parse_vcf(vcf, sample=sample, indels=indels)
		haplotype_bam_writer = None
		if haplotype_bams_prefix is not None:
			haplotype_bam_writer = HaplotypeBamWriter(bam, haplotype_bams_prefix, sample)
		timers.start('parse_vcf')
		for sample, chromosome, variants in vcf_reader:
			variants = remove_overlapping_variants(variants)
			timers.stop('parse_vcf')
			logger.info('Working on chromosome %s', chromosome)
			logger.info('Read %d variants', len(variants))
			bam_sample = None if ignore_read_groups else sample
			logger.info('Reading the BAM file ...')
			try:
				timers.start('read_bam')
				reads = bam_reader.read(chromosome, variants, bam_sample)
				bam_time = timers.stop('read_bam')
			except SampleNotFoundError:
				logger.error("Sample %r is not among the read groups (RG tags) "
					"in the BAM header.", bam_sample)
				sys.exit(1)
			except ReadSetError as e:
				logger.error("%s", e)
				sys.exit(1)
			logger.info('Read %d reads in %.1f s', len(reads), bam_time)

			with timers('slice'):
				# Sort the variants stored in each read
				# TODO: Check whether this is already ensured by construction
				for read in reads:
					read.sort()
				# Sort reads in read set by position
				reads.sort()

				selected_reads, uninformative_read_count = readselection(reads, max_coverage)
				sliced_reads = reads.subset(selected_reads)

				position_list = reads.get_positions()
				accessible_positions = sliced_reads.get_positions()
				informative_read_count = len(reads) - uninformative_read_count
				unphasable_snps = len(position_list) - len(accessible_positions)
				logger.info('%d variants are covered by at least one read', len(position_list))
				logger.info('Skipped %d reads that only cover one variant', uninformative_read_count)
				if position_list:
					logger.info('%d out of %d variant positions (%.1d%%) do not have a read '
						'connecting them to another variant and are thus unphasable',
						unphasable_snps, len(position_list),
						100. * unphasable_snps / len(position_list)
					)
				if reads:
					logger.info('After read selection: Using %d of %d (%.1f%%) reads that cover two or more variants',
						len(selected_reads), informative_read_count, (100. * len(selected_reads) / informative_read_count if informative_read_count > 0 else float('nan'))
					)


			n_best_case_blocks, n_best_case_nonsingleton_blocks = best_case_blocks(reads)
			n_best_case_blocks_cov, n_best_case_nonsingleton_blocks_cov = best_case_blocks(sliced_reads)
			stats.n_best_case_blocks += n_best_case_blocks
			stats.n_best_case_nonsingleton_blocks += n_best_case_nonsingleton_blocks
			stats.n_best_case_blocks_cov += n_best_case_blocks_cov
			stats.n_best_case_nonsingleton_blocks_cov += n_best_case_nonsingleton_blocks_cov
			logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
				n_best_case_nonsingleton_blocks, n_best_case_blocks)
			logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
				n_best_case_nonsingleton_blocks_cov, n_best_case_blocks_cov)

			with timers('phase'):
				logger.info('Phasing the variants (using %d reads)...', len(sliced_reads))
				# Run the core algorithm: construct DP table ...
				dp_table = DPTable(sliced_reads, all_heterozygous)
				# get the mec score
				mec_score = dp_table.get_optimal_cost()
				# ... and do the backtrace to get the solution
				superreads = dp_table.get_super_reads()

				# output the MEC score of phasing
				logger.info('MEC score of phasing: %d', mec_score)

				n_homozygous = sum(1 for v1, v2 in zip(*superreads)
					if v1.allele == v2.allele and v1.allele in (0, 1))
				stats.n_homozygous += n_homozygous

			components = find_components(superreads, sliced_reads)
			n_phased_blocks = len(set(components.values()))
			stats.n_phased_blocks += n_phased_blocks
			logger.info('No. of phased blocks: %d', n_phased_blocks)
			if all_heterozygous:
				assert n_homozygous == 0
			else:
				logger.info('No. of heterozygous variants determined to be homozygous: %d', n_homozygous)

			with timers('write_vcf'):
				vcf_writer.write(chromosome, sample, superreads, components)

			if haplotype_bam_writer is not None:
				logger.info('Writing used reads to haplotype-specific BAM files')
				haplotype_bam_writer.write(sliced_reads, dp_table.get_optimal_partitioning(), chromosome)
			logger.info('Chromosome %s finished', chromosome)
			timers.start('parse_vcf')
		timers.stop('parse_vcf')
	logger.info('== SUMMARY ==')
	logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
		stats.n_best_case_nonsingleton_blocks, stats.n_best_case_blocks)
	logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
		stats.n_best_case_nonsingleton_blocks_cov, stats.n_best_case_blocks_cov)
	if all_heterozygous:
		assert stats.n_homozygous == 0
	else:
		logger.info('No. of heterozygous variants determined to be homozygous: %d', stats.n_homozygous)
	timers.stop('overall')
	logger.info('Time spent reading BAM: %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF: %6.1f s', timers.elapsed('parse_vcf'))
	logger.info('Time spent slicing:     %6.1f s', timers.elapsed('slice'))
	logger.info('Time spent phasing:     %6.1f s', timers.elapsed('phase'))
	logger.info('Time spent writing VCF: %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent on rest:     %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:     %6.1f s', timers.elapsed('overall'))


class NiceFormatter(logging.Formatter):
	"""
	Do not prefix "INFO:" to info-level log messages (but do it for all other
	levels).

	Based on http://stackoverflow.com/a/9218261/715090 .
	"""
	def format(self, record):
		if record.levelno != logging.INFO:
			record.msg = '{}: {}'.format(record.levelname, record.msg)
		return super().format(record)


def main():
	ensure_pysam_version()
	logger.setLevel(logging.INFO)
	handler = logging.StreamHandler()
	handler.setFormatter(NiceFormatter())
	logger.addHandler(handler)
	parser = ArgumentParser(prog='whatshap', description=__doc__)
	parser.add_argument('--version', action='version', version=__version__)
	parser.add_argument('--debug', action='store_true', default=False,
		help='Show more verbose output')
	parser.add_argument('-o', '--output', default=sys.stdout,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	parser.add_argument('--seed', default=123, type=int, help='Random seed (default: %(default)s)')
	parser.add_argument('--indels', dest='indels', default=False, action='store_true',
		help='Also phase indels (default: do not phase indels)')
	parser.add_argument('--distrust-genotypes', dest='all_heterozygous',
		action='store_false', default=True,
		help='Allow switching variants from hetero- to homozygous in an '
		'optimal solution (see documentation).')
	parser.add_argument('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM header and assume all reads come '
		'from the same sample.')
	parser.add_argument('--sample', metavar='SAMPLE', default=None,
		help='Name of a sample to phase. If not given, the first sample in the '
		'input VCF is phased.')
	parser.add_argument('--haplotype-bams', metavar='PREFIX', dest='haplotype_bams_prefix', default=None,
		help='Write reads that have been used for phasing to haplotype-specific BAM files. '
		'Creates PREFIX.1.bam and PREFIX.2.bam')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('bam', nargs='+', metavar='BAM', help='BAM file')

	args = parser.parse_args()
	logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
	del args.debug
	run_whatshap(**vars(args))
