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
import math
from collections import defaultdict, Counter, namedtuple
from contextlib import contextmanager
try:
	from contextlib import ExitStack
except ImportError:
	from contextlib2 import ExitStack  # PY32
from ..vcf import parse_vcf, PhasedVcfWriter, remove_overlapping_variants, VcfVariant
from .. import __version__
from ..args import HelpfulArgumentParser as ArgumentParser
from ..core import Read, ReadSet, DPTable, readselection
from ..graph import ComponentFinder
from ..coverage import CovMonitor
from ..pedigree import mendelian_conflict, load_genetic_map, recombination_cost_map
from ..bam import MultiBamReader, SampleBamReader, BamIndexingError, SampleNotFoundError, HaplotypeBamWriter


__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"


logger = logging.getLogger(__name__)

large_value = 20000
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
	def __init__(self, paths, source_id, mapq_threshold=20):
		self._mapq_threshold = mapq_threshold
		assert len(paths) == 1
		self._reader = SampleBamReader(paths[0], source_id=source_id)

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
			assert 0 < len(readlist) <= 2
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


def find_components(phased_positions, reads):
	"""
	Return a dict that maps each position to the component it is in. A
	component is identified by the position of its leftmost variant.
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
	components = { position : component_finder.find(position) for position in phased_positions }
	return components


def find_recombination(transmission_vector, components, positions, recombcost, recombination_list_filename=None):
	assert len(transmission_vector) == len(positions) == len(recombcost)
	assert set(components.keys()).issubset(set(positions))
	position_to_index = { pos: i for i, pos in enumerate(positions) }
	blocks = defaultdict(list)
	for position, block_id in components.items():
		blocks[block_id].append(position)

	RecombinationEvent = namedtuple('RecombinationEvent', ['position1', 'position2', 'inheritance_value1', 'inheritance_value2' , 'recombination_cost'])
	event_count = 0
	event_list = []
	cum_recomb_cost = 0
	for block_id, block in blocks.items():
		block.sort()
		block_transmission_vector = [ transmission_vector[position_to_index[i]] for i in block ]
		block_recomb_cost = [ recombcost[position_to_index[i]] for i in block ]
		if len(block) <= 2:
			continue
		for i in range(2, len(block)):
			if block_transmission_vector[i-1] != block_transmission_vector[i]:
				event_list.append(RecombinationEvent(block[i-1], block[i], block_transmission_vector[i-1], block_transmission_vector[i], block_recomb_cost[i]))
				cum_recomb_cost += block_recomb_cost[i]
				event_count += 1

	if recombination_list_filename is not None:
		event_list.sort()
		f = open(recombination_list_filename, 'w')
		print('#position1', 'position2', 'inheritance_value1', 'inheritance_value2' , 'recombination_cost', file=f)
		for e in event_list:
			print(e.position1, e.position2, e.inheritance_value1, e.inheritance_value2, e.recombination_cost, file=f)
		f.close()

	logger.info('Cost accounted for by recombination events: %d', cum_recomb_cost)
	return event_count


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


def run_whatshap(chromosome, genmap, bamm, vcfm, bamf, vcff, bamc, vcfc,
		outputm=sys.stdout,outputf=sys.stdout,outputc=sys.stdout, samplem=None, samplef=None, samplec=None, ignore_read_groups=False, indels=True,
		mapping_quality=20, max_coveragem=5, max_coveragef=5, max_coveragec=5, seed=123, haplotype_bams_prefix=None, recombination_list_filename=None):
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

	assert len(chromosome) == 1
	chromosome = chromosome[0]

	assert len(genmap) == 1
	genmap = genmap[0]

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
			bam_readerm = stack.enter_context(ReadSetReader(bamm, source_id=1, mapq_threshold=mapping_quality))
			bam_readerf = stack.enter_context(ReadSetReader(bamf, source_id=2, mapq_threshold=mapping_quality))
			bam_readerc = stack.enter_context(ReadSetReader(bamc, source_id=0, mapq_threshold=mapping_quality))
		except (OSError, BamIndexingError) as e:
			logger.error(e)
			sys.exit(1)
		if outputm is not sys.stdout:
			outputm = stack.enter_context(open(outputm, 'w'))
		if outputf is not sys.stdout:
			outputf = stack.enter_context(open(outputf, 'w'))
		if outputc is not sys.stdout:
			outputc = stack.enter_context(open(outputc, 'w'))
		command_line = '(whatshap {}) {}'.format(__version__ , ' '.join(sys.argv[1:]))

		vcf_writerm = PhasedVcfWriter(command_line=command_line, in_path=' '.join(vcfm), out_file=outputm)
		vcf_writerf = PhasedVcfWriter(command_line=command_line, in_path=' '.join(vcff), out_file=outputf)
		vcf_writerc = PhasedVcfWriter(command_line=command_line, in_path=' '.join(vcfc), out_file=outputc)
		with timers('parse_vcf'):
			vcf_readerm = parse_vcf(' '.join(vcfm), sample=samplem, indels=indels)
			vcf_readerf = parse_vcf(' '.join(vcff), sample=samplef, indels=indels)
			vcf_readerc = parse_vcf(' '.join(vcfc), sample=samplec, indels=indels)
		#TODO: Reactive haplotype BAM writing.
		#haplotype_bam_writer = None
		if haplotype_bams_prefix is not None:
			logger.error('Option --haplotype-bams not implemented right now. Sorry.')
			return 1
		#	haplotype_bam_writer = HaplotypeBamWriter(bamm, haplotype_bams_prefix, samplem)
		#if haplotype_bams_prefix is not None:
			#haplotype_bam_writerm = HaplotypeBamWriter(bamm, haplotype_bams_prefix, sample)
			#haplotype_bam_writerf = HaplotypeBamWriter(bamf, haplotype_bams_prefix, sample)
			#haplotype_bam_writerc = HaplotypeBamWriter(bamc, haplotype_bams_prefix, sample)

		for m,f,c in zip(vcf_readerm,vcf_readerf,vcf_readerc):
			samplem, chromosomem, variantsm = m
			samplef, chromosomef, variantsf = f
			samplec, chromosomec, variantsc = c
			assert chromosomem == chromosomef == chromosomec, 'Input VCFs out of sync. They have to contain the same variants.'
			if chromosome != chromosomem:
				logger.info('Skipping chromosome %s found in VCF files.', chromosomem)
				continue
			logger.info('Processing chromosome %s of trio mother=%s, father=%s, child=%s', chromosome, samplem, samplef, samplec)

			with timers('variant_filtering'):
				variantsm = remove_overlapping_variants(variantsm)
				variantsf = remove_overlapping_variants(variantsf)
				variantsc = remove_overlapping_variants(variantsc)
				assert len(variantsm) == len(variantsf) == len(variantsc), 'We are assuming that all input VCFs contain the same variants'
				logger.info('Number of variants: %d', len(variantsm))
				
				variants = []
				mendelian_conflicts = 0
				for variantm, variantf, variantc in zip(variantsm, variantsf, variantsc):
					assert variantm.position == variantf.position == variantc.position, 'Input VCFs out of sync'
					assert variantm.reference_allele == variantf.reference_allele == variantc.reference_allele, 'Input VCFs out of sync'
					assert variantm.alternative_allele == variantf.alternative_allele == variantc.alternative_allele, 'Input VCFs out of sync'
					if (variantm.genotype == 1) or (variantf.genotype == 1) or (variantc.genotype == 1):
						if mendelian_conflict(variantm.genotype, variantf.genotype, variantc.genotype):
							mendelian_conflicts += 1
						else:
							variants.append(VcfVariant(
								position = variantm.position,
								reference_allele = variantm.reference_allele,
								alternative_allele = variantm.alternative_allele,
								genotype = [variantm.genotype, variantf.genotype, variantc.genotype] )
							)
				logger.info('Number of variants skipped due to Mendelian conflicts: %d', mendelian_conflicts)
				logger.info('Number of remaining variants hetorzygous in at least one individual: %d', len(variants))
			

			with timers('read_bam'):
				bam_samplem = None if ignore_read_groups else samplem
				bam_samplef = None if ignore_read_groups else samplef
				bam_samplec = None if ignore_read_groups else samplec
				logger.info('Reading the BAM file ...')
				try:
					readsm = bam_readerm.read(chromosomem, variants, bam_samplem)
					readsf = bam_readerf.read(chromosomef, variants, bam_samplef)
					readsc = bam_readerc.read(chromosomec, variants, bam_samplec)
					logger.info('Variant informative reads: mother=%d, father=%d, child=%d', len(readsm), len(readsf),len(readsc))
				except SampleNotFoundError:
					#logger.error("Sample %r is not among the read groups (RG tags) "
						#"in the BAM header.", bam_sample)
					sys.exit(1)
				except ReadSetError as e:
					logger.error("%s", e)
					sys.exit(1)

			with timers('slice'):
				# Sort the variants stored in each read
				# TODO: Check whether this is already ensured by construction
				for read in readsm:
					read.sort()
				for read in readsf:
					read.sort()
				for read in readsc:
					read.sort()
				# Sort reads in read set by position
				readsm.sort()
				readsf.sort()
				readsc.sort()
				# TODO: Read selection done w.r.t. all variants, where using heterozygous variants only 
				# TODO: would probably give better results.
				selected_readsm, uninformative_read_countm = readselection(readsm, max_coveragem)
				logger.info('Done selecting reads in mother, selected %d of %d reads', len(selected_readsm), len(readsm))
				selected_readsf, uninformative_read_countf = readselection(readsf, max_coveragef)
				logger.info('Done selecting reads in father, selected %d of %d reads', len(selected_readsf), len(readsf))
				selected_readsc, uninformative_read_countc = readselection(readsc, max_coveragec)
				logger.info('Done selecting reads in child, selected %d of %d reads', len(selected_readsc), len(readsc))
				sliced_readsm = readsm.subset(selected_readsm)
				sliced_readsf = readsf.subset(selected_readsf)
				sliced_readsc = readsc.subset(selected_readsc)
				accessible_positions = set()
				for sliced_reads in [sliced_readsm, sliced_readsf, sliced_readsc]:
					accessible_positions.update(set(sliced_reads.get_positions()))
				accessible_positions = list(accessible_positions)
				accessible_positions.sort()
				logger.info('Variants covered by at least one phase-informative read in at least individual after read selection: %d', len(accessible_positions))
				
				# Create genotype lists for each individual
				accessible_positions_set = set(accessible_positions)
				genotypesm, genotypesf, genotypesc = [], [], []
				for variant in variants:
					if variant.position not in accessible_positions_set:
						continue
					assert len(variant.genotype) == 3
					genotypesm.append(variant.genotype[0])
					genotypesf.append(variant.genotype[1])
					genotypesc.append(variant.genotype[2])
				assert len(genotypesm) == len(accessible_positions)
				
			with timers('phase'):
				#logger.info('Phasing the variants (using %d reads)...', len(sliced_reads))
				# Run the core algorithm: construct DP table ...
				#intersect=list(set(accessible_positionsc).intersection(set(accessible_positionsf).intersection(set(accessible_positionsm))))
				allreads = ReadSet()
				for sliced_reads in [sliced_readsm, sliced_readsf, sliced_readsc]:
					sliced_reads.sort()
					for read in sliced_reads:
						read.sort()
						allreads.add(read)
				allreads.sort()
				read_marks = [ read.source_id for read in allreads ]
				recombcost = recombination_cost_map(load_genetic_map(genmap), accessible_positions)
				dp_table = DPTable(allreads, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
				superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()
				assert len(superreadsm[0]) == len(superreadsf[0]) == len(superreadsc[0]) == len(transmission_vector)
				logger.info('PedMEC cost: %d', dp_table.get_optimal_cost())

			with timers('components'):
				components = find_components(accessible_positions, allreads)
				
				n_phased_blocks = len(set(components.values()))
				#stats.n_phased_blocks += n_phased_blocks
				logger.info('No. of phased blocks: %d', n_phased_blocks)
				n_recombination = find_recombination(transmission_vector, components, accessible_positions, recombcost, recombination_list_filename)
				logger.info('No. of detected recombination events: %d', n_recombination)

			with timers('write_vcf'):
				vcf_writerm.write(chromosomem, samplem, superreadsm, components)
				vcf_writerf.write(chromosomef, samplef, superreadsf, components)
				vcf_writerc.write(chromosomec, samplec, superreadsc, components)

			#if haplotype_bam_writer is not None:
				#logger.info('Writing used reads to haplotype-specific BAM files')
				#haplotype_bam_writer.write(sliced_reads, dp_table.get_optimal_partitioning(), chromosome)
			#logger.info('Chromosome %s finished', chromosome)
			#timers.start('parse_vcf')
		#timers.stop('parse_vcf')
	#logger.info('== SUMMARY ==')
	#logger.info('Best-case phasing would result in %d non-singleton phased blocks (%d in total)',
		#stats.n_best_case_nonsingleton_blocks, stats.n_best_case_blocks)
	#logger.info('... after read selection: %d non-singleton phased blocks (%d in total)',
		#stats.n_best_case_nonsingleton_blocks_cov, stats.n_best_case_blocks_cov)
	#if all_heterozygous:
		#assert stats.n_homozygous == 0
	#else:
		#logger.info('No. of heterozygous variants determined to be homozygous: %d', stats.n_homozygous)
	timers.stop('overall')
	logger.info('Time spent reading BAM:        %6.1f s', timers.elapsed('read_bam'))
	logger.info('Time spent parsing VCF:        %6.1f s', timers.elapsed('parse_vcf'))
	logger.info('Time spent filtering variants: %6.1f s', timers.elapsed('variant_filtering'))
	logger.info('Time spent slicing:            %6.1f s', timers.elapsed('slice'))
	logger.info('Time spent phasing:            %6.1f s', timers.elapsed('phase'))
	logger.info('Time spent finding components: %6.1f s', timers.elapsed('components'))
	logger.info('Time spent writing VCF:        %6.1f s', timers.elapsed('write_vcf'))
	logger.info('Time spent on rest:            %6.1f s', 2 * timers.elapsed('overall') - timers.total())
	logger.info('Total elapsed time:            %6.1f s', timers.elapsed('overall'))


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
	parser.add_argument('-om', '--outputm', default=sys.stdout,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('-of', '--outputf', default=sys.stdout,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('-oc', '--outputc', default=sys.stdout,
		help='Output VCF file. If omitted, use standard output.')
	parser.add_argument('--max-coveragem', '-M', metavar='MAXCOVM', default=5, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--max-coveragef', '-F', metavar='MAXCOVF', default=5, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--max-coveragec', '-C', metavar='MAXCOVC', default=5, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--mapping-quality', '--mapq', metavar='QUAL',
		default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
	parser.add_argument('--seed', default=123, type=int, help='Random seed (default: %(default)s)')
	parser.add_argument('--indels', dest='indels', default=False, action='store_true',
		help='Also phase indels (default: do not phase indels)')
	parser.add_argument('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM header and assume all reads come '
		'from the same sample.')
	parser.add_argument('--samplem', metavar='SAMPLEM', default=None,
		help='Name of a sample to phase. If not given, the first sample in the '
		'input VCF is phased.')
	parser.add_argument('--samplef', metavar='SAMPLEF', default=None,
		help='Name of a sample to phase. If not given, the first sample in the '
		'input VCF is phased.')
	parser.add_argument('--samplec', metavar='SAMPLEC', default=None,
		help='Name of a sample to phase. If not given, the first sample in the '
		'input VCF is phased.')
	parser.add_argument('--haplotype-bams', metavar='PREFIX', dest='haplotype_bams_prefix', default=None,
		help='Write reads that have been used for phasing to haplotype-specific BAM files. '
		'Creates PREFIX.1.bam and PREFIX.2.bam')
	parser.add_argument('--recombination-list', metavar='RECOMBLIST', dest='recombination_list_filename', default=None,
		help='Write putative recombination events to given filename.')
	parser.add_argument('chromosome', nargs=1, metavar='CHROM', help='Chromosome to work on')
	parser.add_argument('genmap', nargs=1, metavar='GENMAP', help='File with genetic map')
	parser.add_argument('vcfm', nargs=1, metavar='VCFM', help='VCF file')
	#parser.add_argument('bamm', nargs='+', metavar='BAMM', help='BAM file')
	parser.add_argument('bamm', nargs=1, metavar='BAMM', help='BAM file')
	parser.add_argument('vcff', nargs=1, metavar='VCFF', help='VCF file')
	parser.add_argument('bamf', nargs=1,metavar='BAMF', help='BAM file')
	parser.add_argument('vcfc', nargs=1,metavar='VCFC', help='VCF file')
	parser.add_argument('bamc',  nargs=1,metavar='BAMC', help='BAM file')

	args = parser.parse_args()
	logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
	del args.debug
	run_whatshap(**vars(args))

