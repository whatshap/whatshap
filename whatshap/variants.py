"""
Detect variants in reads.
"""
from collections import defaultdict, Counter
import logging
from .core import Read, ReadSet
from .bam import SampleBamReader, MultiBamReader
from .align import edit_distance, edit_distance_affine_gap
from ._variants import _iterate_cigar


logger = logging.getLogger(__name__)


class ReadSetError(Exception):
	pass


class ReadSetReader:
	"""
	Associate VCF variants with BAM reads.

	A VCF file contain variants, and a BAM file contain reads, but the
	information which read contains which variant is not available. This
	class re-discovers the variants in each read, using the
	knowledge in the VCF of where they should occur.
	"""

	def __init__(self, paths, reference, numeric_sample_ids, mapq_threshold=20, overhang=10, affine=False, gap_start=10, gap_extend=7, default_mismatch=15):
		"""
		paths -- list of BAM paths
		reference -- path to reference FASTA (can be None)
		numeric_sample_ids -- ??
		mapq_threshold -- minimum mapping quality
		overhang -- extend alignment by this many bases to left and right
		affine -- use affine gap costs
		gap_start, gap_extend, default_mismatch -- parameters for affine gap cost alignment
		"""
		self._mapq_threshold = mapq_threshold
		self._numeric_sample_ids = numeric_sample_ids
		self._use_affine = affine
		self._gap_start = gap_start
		self._gap_extend = gap_extend
		self._default_mismatch = default_mismatch
		self._overhang = overhang
		if len(paths) == 1:
			self._reader = SampleBamReader(paths[0], reference=reference)
		else:
			self._reader = MultiBamReader(paths, reference=reference)

	def read(self, chromosome, variants, sample, reference):
		"""
		Return a ReadSet object representing the given variants.

		If a reference is provided (reference is not None), alleles are
		detected by re-aligning	sections of the query to the REF and ALT
		sequence extended a few bases to the left and right.

		If reference is None, alleles are detected by inspecting the
		existing alignment (via the CIGAR).

		chromosome -- name of chromosome to work on
		variants -- list of vcf.VcfVariant objects
		sample -- name of sample to work on. If None, read group information is
			ignored and all reads in the file are used.
		reference -- reference sequence of the given chromosome (or None)
		"""
		# Since variants are identified by position, positions must be unique.
		if __debug__ and variants:
			varposc = Counter(variant.position for variant in variants)
			pos, count = varposc.most_common()[0]
			assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

		alignments = self._usable_alignments(chromosome, sample)
		reads = self._alignments_to_readdict(alignments, variants, sample, reference)
		return self._readdict_to_readset(reads)

	def _readdict_to_readset(self, reads):
		"""
		reads is a dict that maps read names to Read objects

		TODO this functionality should be within ReadSet
		"""
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

	def _usable_alignments(self, chromosome, sample):
		"""
		Retrieve usable (suficient mapping quality, not secondary etc.)
		alignments from the alignment file
		"""
		for alignment in self._reader.fetch(reference=chromosome, sample=sample):
			# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
			if alignment.bam_alignment.flag & 2048 != 0:
				continue
			if alignment.bam_alignment.mapping_quality < self._mapq_threshold:
				continue
			if alignment.bam_alignment.is_secondary:
				continue
			if alignment.bam_alignment.is_unmapped:
				continue
			if alignment.bam_alignment.is_duplicate:
				continue
			yield alignment

	def has_reference(self, chromosome):
		return self._reader.has_reference(chromosome)

	def _alignments_to_readdict(self, alignments, variants, sample, reference):
		"""
		Convert BAM alignments to Read objects.

		If reference is not None, alleles are detected through re-alignment.

		Return a dict that maps read names to lists of Read objects. Each list
		has two entries for paired-end reads, one entry for single-end reads.
		"""
		# FIXME hard-coded zero
		numeric_sample_id = 0 if sample is None else self._numeric_sample_ids[sample]
		reads = defaultdict(list)
		if reference is not None:
			# Copy the pyfaidx.FastaRecord into a str for faster access
			reference = reference[:]
			normalized_variants = variants
		else:
			normalized_variants = [ variant.normalized() for variant in variants ]

		i = 0  # index into variants
		for alignment in alignments:
			# Skip variants that are to the left of this read
			while i < len(normalized_variants) and normalized_variants[i].position < alignment.bam_alignment.reference_start:
				i += 1
			
			BX_tag = ''
			if alignment.bam_alignment.has_tag('BX'):
				BX_tag = alignment.bam_alignment.get_tag('BX')

			read = Read(alignment.bam_alignment.qname,
					alignment.bam_alignment.mapq, alignment.source_id,
					numeric_sample_id, alignment.bam_alignment.reference_start, BX_tag)

			if reference is None:
				detected = self.detect_alleles(normalized_variants, i, alignment.bam_alignment)

			else:
				detected = self.detect_alleles_by_alignment(variants, i, alignment.bam_alignment, reference, self._overhang, self._use_affine, self._gap_start, self._gap_extend, self._default_mismatch)
			for j, allele, quality in detected:
				read.add_variant(variants[j].position, allele, quality)
			if read:  # At least one variant covered and detected
				reads[(alignment.source_id, alignment.bam_alignment.qname, numeric_sample_id)].append(read)
		return reads

	@staticmethod
	def detect_alleles(variants, j, bam_read):
		"""
		Detect the correct alleles of the variants that are covered by the
		given bam_read.

		Yield tuples (index, allele, quality), where index is into the variants list.

		variants -- list of variants (VcfVariant objects)
		j -- index of the first variant (in the variants list) to check
		"""
		ref_pos = bam_read.reference_start  # position relative to reference
		query_pos = 0  # position relative to read

		seen_positions = set()
		for cigar_op, length in bam_read.cigartuples:
			# Skip variants that come before this region
			while j < len(variants) and variants[j].position < ref_pos:
				j += 1

			# The mapping of CIGAR operators to numbers is:
			# MIDNSHPX= => 012345678
			if cigar_op in (0, 7, 8):  # M, X, = operators (match)
				# Iterate over all variants that are in this region
				while j < len(variants) and variants[j].position < ref_pos + length:
					if len(variants[j].reference_allele) == len(variants[j].alternative_allele) == 1:
						# Variant is a SNV
						offset = variants[j].position - ref_pos
						base = bam_read.query_sequence[query_pos + offset]
						if base == variants[j].reference_allele:
							allele = 0
						elif base == variants[j].alternative_allele:
							allele = 1
						else:
							allele = None
						if allele is not None:
							# TODO
							# Fix this: we can actually have indel and SNV
							# calls at identical positions. For now, ignore the
							# second variant.
							if variants[j].position in seen_positions:
								logger.debug("Found two variants at identical positions. Ignoring the second one: %s",
									variants[j])
							else:
								# Do not use bam_read.qual here as it is extremely slow.
								# If we ever decide to be compatible with older pysam
								# versions, cache bam_read.qual somewhere - do not
								# access it within this loop (3x slower otherwise).
								if bam_read.query_qualities:
									qual = bam_read.query_qualities[query_pos + offset]
								else:
									qual = 30  # TODO
								yield (j, allele, qual)
								seen_positions.add(variants[j].position)
					elif len(variants[j].reference_allele) == 0:
						assert len(variants[j].alternative_allele) > 0
						# This variant is an insertion. Since we are in a region of
						# matches, the insertion was *not* observed (reference allele).
						qual = 30  # TODO average qualities of "not inserted" bases?
						yield (j, 0, qual)
						seen_positions.add(variants[j].position)
					elif len(variants[j].alternative_allele) == 0:
						assert len(variants[j].reference_allele) > 0
						# This variant is a deletion that was not observed.
						# Add it only if the next variant is not located 'within'
						# the deletion.
						deletion_end = variants[j].position + len(variants[j].reference_allele)
						if not (j + 1 < len(variants) and variants[j + 1].position < deletion_end):
							qual = 30  # TODO
							yield (j, 0, qual)
							seen_positions.add(variants[j].position)
						else:
							logger.info('Skipped a deletion overlapping another variant at pos. %d',
								variants[j].position)
							# Also skip all variants that this deletion overlaps
							while j + 1 < len(variants) and variants[j + 1].position < deletion_end:
								j += 1
							# One additional j += 1 is done below
					else:
						assert False, "Strange variant: {}".format(variants[j])
					j += 1
				query_pos += length
				ref_pos += length
			elif cigar_op == 1:  # I operator (insertion)
				if j < len(variants) and variants[j].position == ref_pos and \
						len(variants[j].reference_allele) == 0 and \
						variants[j].alternative_allele == bam_read.query_sequence[query_pos:query_pos + length]:
					qual = 30  # TODO
					assert variants[j].position not in seen_positions
					yield (j, 1, qual)
					seen_positions.add(variants[j].position)
					j += 1
				query_pos += length
			elif cigar_op == 2:  # D operator (deletion)
				# We only check the length of the deletion, not the sequence
				# that gets deleted since we don’t have the reference available.
				# (We could parse the MD tag if it exists.)
				if j < len(variants) and variants[j].position == ref_pos and \
						len(variants[j].alternative_allele) == 0 and \
						len(variants[j].reference_allele) == length:
					qual = 30  # TODO
					deletion_end = variants[j].position + len(variants[j].reference_allele)
					if not (j + 1 < len(variants) and variants[j + 1].position < deletion_end):
						qual = 30  # TODO
						assert variants[j].position not in seen_positions
						yield (j, 1, qual)
						seen_positions.add(variants[j].position)
					else:
						logger.info('Skipped a deletion overlapping another variant at pos. %d', variants[j].position)
						# Also skip all variants that this deletion overlaps
						while j + 1 < len(variants) and variants[j + 1].position < deletion_end:
							j += 1
						# One additional j += 1 is done below
					j += 1
				ref_pos += length
			elif cigar_op == 3:  # N operator (reference skip)
				ref_pos += length
			elif cigar_op == 4:  # S operator (soft clipping)
				query_pos += length
			elif cigar_op == 5 or cigar_op == 6:  # H or P (hard clipping or padding)
				pass
			else:
				logger.error("Unsupported CIGAR operation: %d", cigar_op)
				raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))

	@staticmethod
	def split_cigar(cigar, i, consumed):
		"""
		Split a CIGAR into two parts. i and consumed describe the split position.
		i is the element of the cigar list that should be split, and consumed says
		at how many operations to split within that element.

		The CIGAR is given as a list of (operation, length) pairs.

		i -- split at this index in cigar list
		consumed -- how many cigar ops at cigar[i] are to the *left* of the
			split position

		Return a tuple (left, right).

		Example:
		Assume the cigar is 3M 1D 6M 2I 4M.
		With i == 2 and consumed == 5, the cigar is split into
		3M 1D 5M and 1M 2I 4M.
		"""
		middle_op, middle_length = cigar[i]
		assert consumed <= middle_length
		if consumed > 0:
			left = cigar[:i] + [(middle_op, consumed)]
		else:
			left = cigar[:i]
		if consumed < middle_length:
			right = [(middle_op, middle_length-consumed)] + cigar[i+1:]
		else:
			right = cigar[i+1:]
		return left, right

	@staticmethod
	def cigar_prefix_length(cigar, reference_bases):
		"""
		Given a prefix of length reference_bases relative to the reference, how
		long is the prefix of the read? In other words: If reference_bases on
		the reference are consumed, how many bases on the query does that
		correspond to?

		If the position is within or at the end of an insertion (which do not
		consume bases on the reference), then the number of bases up to the
		beginning of the insertion is reported.

		Return a pair (reference_bases, query_bases) where the value for
		reference_bases may be smaller than the requested one if the CIGAR does
		not cover enough reference bases.

		Reference skips (N operators) are treated as the end of the read. That
		is, no positions beyond a reference skip are reported.
		"""
		ref_pos = 0
		query_pos = 0
		for op, length in cigar:
			if op in (0, 7, 8):  # M, X, =
				ref_pos += length
				query_pos += length
				if ref_pos >= reference_bases:
					return (reference_bases, query_pos + reference_bases - ref_pos)
			elif op == 2:  # D
				ref_pos += length
				if ref_pos >= reference_bases:
					return (reference_bases, query_pos)
			elif op == 1:  # I
				query_pos += length
			elif op == 4 or op == 5:  # soft or hard clipping
				pass
			elif op == 3:  # N
				# Always stop at reference skips
				return (reference_bases, query_pos)
			else:
				assert False, 'unknown CIGAR operator'
		assert ref_pos < reference_bases
		return (ref_pos, query_pos)

	@staticmethod
	def realign(variant, bam_read, cigartuples, i, consumed, query_pos, reference, overhang, use_affine, gap_start, gap_extend,default_mismatch):
		"""
		Realign a read to the two alleles of a single variant.
		i and consumed describe where to split the cigar into a part before the
		variant position and into a part starting at the variant position, see split_cigar().

		variant -- VcfVariant
		bam_read -- the AlignedSegment
		cigartuples -- the AlignedSegment.cigartuples property (accessing it is expensive, so re-use it)
		i, consumed -- see split_cigar method
		query_pos -- index of the query base that is at the variant position
		reference -- the reference as a str-like object (full chromosome)
		overhang -- extend alignment by this many bases to left and right
		use_affine -- if true, use affine gap costs for realignment
		gap_start, gap_extend -- if affine_gap=true, use these parameters for affine gap cost alignment
		default_mismatch -- if affine_gap=true, use this as mismatch cost in case no base qualities are in bam
		"""
		# Do not process symbolic alleles like <DEL>, <DUP>, etc.
		if variant.alternative_allele.startswith('<'):
			return None, None

		left_cigar, right_cigar = ReadSetReader.split_cigar(cigartuples, i, consumed)

		left_ref_bases, left_query_bases = ReadSetReader.cigar_prefix_length(left_cigar[::-1], overhang)
		right_ref_bases, right_query_bases = ReadSetReader.cigar_prefix_length(right_cigar,
			len(variant.reference_allele) + overhang)

		assert variant.position - left_ref_bases >= 0
		assert variant.position + right_ref_bases <= len(reference)

		query = bam_read.query_sequence[query_pos-left_query_bases:query_pos+right_query_bases]
		ref = reference[variant.position - left_ref_bases:variant.position + right_ref_bases]
		alt = reference[variant.position - left_ref_bases:variant.position] + variant.alternative_allele + \
				reference[variant.position+len(variant.reference_allele):variant.position + right_ref_bases]

		distance_ref = 0
		distance_alt = 0

		base_qual_score = 30

		if use_affine:
			assert gap_start != None
			assert gap_extend != None
			assert default_mismatch != None

			# get base qualities if present (to be used as mismatch costs)
			base_qualities = [default_mismatch] * len(query)
			#if bam_read.query_qualities != None:
			#	base_qualities = bam_read.query_qualities[query_pos-left_query_bases:query_pos+right_query_bases]

			# compute edit dist. with affine gap costs using base qual. as mismatch cost
			distance_ref = edit_distance_affine_gap(query,ref,base_qualities,gap_start,gap_extend)
			distance_alt = edit_distance_affine_gap(query,alt,base_qualities,gap_start,gap_extend)
			base_qual_score = abs(distance_ref-distance_alt)

		else:
			distance_ref = edit_distance(query, ref)
			distance_alt = edit_distance(query, alt)

		if distance_ref < distance_alt:
			return 0, base_qual_score  # detected REF
		elif distance_ref > distance_alt:
			return 1, base_qual_score  # detected ALT
		else:
			return None, None  # cannot decide

	@staticmethod
	def detect_alleles_by_alignment(variants, j, bam_read, reference, overhang=10, use_affine=False, gap_start=None, gap_extend=None, default_mismatch=None):
		"""
		Detect which alleles the given bam_read covers. Detect the correct
		alleles of the variants that are covered by the given bam_read.

		Yield tuples (position, allele, quality).

		variants -- list of variants (VcfVariant objects)
		j -- index of the first variant (in the variants list) to check
		"""
		# Accessing bam_read.cigartuples is expensive, do it only once
		cigartuples = bam_read.cigartuples

		# For the same reason, the following check is here instad of
		# in the _usable_alignments method
		if not cigartuples:
			return

		for index, i, consumed, query_pos in _iterate_cigar(variants, j, bam_read, cigartuples):
			allele, quality = ReadSetReader.realign(variants[index], bam_read, cigartuples, i,
			            consumed, query_pos, reference, overhang, use_affine, gap_start, gap_extend, default_mismatch)
			if allele in (0, 1):
				yield (index, allele, quality)  # TODO quality???

	@staticmethod
	def _merge_pair(read1, read2):
		"""
		Merge the two ends of a paired-end read into a single core.Read. Also
		takes care of self-overlapping read pairs.

		TODO this can be simplified as soon as a variant in a read can be
		modified.
		"""
		if read2:
			result = Read(read1.name, read1.mapqs[0], read1.source_id, read1.sample_id, read1.reference_start, read1.BX_tag)
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
