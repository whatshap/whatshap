"""
Detect variants in reads.
"""
from collections import defaultdict, Counter
import logging
from .core import Read, ReadSet
from .bam import SampleBamReader, MultiBamReader

logger = logging.getLogger(__name__)


def covered_variants(variants, start, bam_read, source_id, numeric_sample_id):
	"""
	Find the variants that are covered by the given bam_read and return a
	core.Read instance that represents those variants. The instance may be
	empty.

	start -- index of the first variant (in the variants list) to check
	"""
	core_read = Read(bam_read.qname, bam_read.mapq, source_id, numeric_sample_id)
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
					qual = 30  # TODO average qualities of "not inserted" bases?
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
			raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))
	return core_read


class ReadSetError(Exception):
	pass


class ReadSetReader:
	"""
	Associate VCF variants with BAM reads.

	VCF file contain variants, and BAM file contain reads, but the
	information which read contains which variant is not available. This
	class re-discovers the variants in each read, using the
	knowledge in the VCF of where they should occur.
	"""
	def __init__(self, paths, numeric_sample_ids, mapq_threshold=20):
		"""
		paths -- list of BAM paths
		numeric_sample_ids -- ??
		mapq_threshold -- minimum mapping quality
		"""
		self._mapq_threshold = mapq_threshold
		self._numeric_sample_ids = numeric_sample_ids
		if len(paths) == 1:
			self._reader = SampleBamReader(paths[0])
		else:
			self._reader = MultiBamReader(paths)

	def read(self, chromosome, variants, sample):
		"""

		chromosome -- name of chromosome to work on
		variants -- list of vcf.VcfVariant objects
		sample -- name of sample to work on. If None, read group information is
			ignored and all reads in the file are used.

		Return a ReadSet object.
		"""
		# Since variants are identified by position, positions must be unique.
		if __debug__ and variants:
			varposc = Counter(variant.position for variant in variants)
			pos, count = varposc.most_common()[0]
			assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

		reads = self._fetch_all_reads(chromosome, variants, sample)

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

	def _fetch_all_reads(self, chromosome, variants, sample):
		"""
		Return a dict that maps read names to lists of Read objects.

		Each list has two entries paired-end reads, one entry for single-end reads.
		"""
		reads = defaultdict(list)

		# Retrieve all reads
		i = 0  # index into variants
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

			read = covered_variants(variants, i, alignment.bam_alignment, alignment.source_id,
				self._numeric_sample_ids[sample])
			# Only add if it covers at least one variant
			if read:
				reads[(alignment.source_id, alignment.bam_alignment.qname)].append(read)
		return reads

	@staticmethod
	def _merge_pair(read1, read2):
		"""
		Merge the two ends of a paired-end read into a single core.Read. Also
		takes care of self-overlapping read pairs.

		TODO this can be simplified as soon as a variant in a read can be
		modified.
		"""
		if read2:
			result = Read(read1.name, read1.mapqs[0], read1.source_id, read1.numeric_sample_id)
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
