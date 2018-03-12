import os
from urllib.parse import urlparse

import pysam
import logging
import heapq
from collections import defaultdict, namedtuple
import subprocess

from .utils import detect_file_format


logger = logging.getLogger(__name__)

AlignmentWithSourceID = namedtuple('AlignmentWithSourceID', ['source_id', 'bam_alignment'])


class AlignmentFileNotIndexedError(Exception):
	pass


class SampleNotFoundError(Exception):
	pass


class ReferenceNotFoundError(Exception):
	pass


class EmptyAlignmentFileError(Exception):
	pass


def is_local(path):
	return urlparse(path).scheme == ''


class SampleBamReader:
	"""
	A wrapper for Samfile that provides only those reads from a BAM or CRAM file
	that belong to a specified sample. The BAM/CRAM file must have an index
	(bai/crai).
	"""
	def __init__(self, path, *, source_id=0, reference=None):
		"""
		path -- URL or path to BAM or CRAM file
		reference -- optional path to FASTA reference for CRAM
		"""
		self.source_id = source_id
		if reference:
			reference = os.path.abspath(reference)

		self._samfile = pysam.AlignmentFile(path, reference_filename=reference)
		try:
			# TODO
			# multiple_iterators should not be necessary -
			# perhaps a bug in pysam 0.13
			fetcher = self._samfile.fetch(multiple_iterators=True)
		except ValueError:
			raise AlignmentFileNotIndexedError(path)

		# For CRAM files, AlignmentFile.fetch() does not raise an error if the
		# reference could not be retrieved; it just returns an empty list.
		# To detect this situation, we treat all empty alignment files as an error.
		try:
			next(fetcher)
		except StopIteration:
			raise EmptyAlignmentFileError(path) from None
		self._references = frozenset(self._samfile.references)
		self._initialize_sample_to_group_ids()

	def has_reference(self, name):
		return name in self._references

	def _initialize_sample_to_group_ids(self):
		"""
		Create a dictionary that maps a sample name to a set of read group ids.
		"""
		read_groups = self._samfile.header.get('RG', [])  # a list of dicts
		logger.debug('Read groups in CRAM/BAM header: %s', read_groups)
		samples = defaultdict(list)
		for read_group in read_groups:
			if 'SM' in read_group:
				samples[read_group['SM']].append(read_group['ID'])
			else:
				logger.warning('Read group "%s" does not contain an SM field to assign it to a sample. Use --ignore-read-groups to use these alignments anyway.', read_group['ID'])
		self._sample_to_group_ids = {
			id: frozenset(values) for id, values in samples.items() }

	def has_sample(self, sample):
		"""Return whether this file contains reads for the given sample"""
		return sample in self._sample_to_group_ids

	def fetch(self, reference, sample):
		"""
		Yield instances of AlignmentWithSourceID, with source_id value given
		at construction time.
		Raise KeyError if sample not found among samples named in RG header.
		Raise ReferenceNotFoundError when reference is not in the BAM.
		"""
		if reference not in self._references:
			# fetch() (below) raises a ValueError when the reference is invalid,
			# but also on other errors.
			raise ReferenceNotFoundError(reference)
		if sample is None:
			for bam_read in self._samfile.fetch(reference):
				yield AlignmentWithSourceID(self.source_id, bam_read)
		else:
			try:
				read_groups = self._sample_to_group_ids[sample]
			except KeyError:
				raise SampleNotFoundError()
			# TODO
			# The multiple_iterators shouldn’t be necessary, but CRAM files
			# don’t work without it in pysam 0.13
			for bam_read in self._samfile.fetch(reference, multiple_iterators=True):
				if bam_read.opt('RG') in read_groups:
					yield AlignmentWithSourceID(self.source_id, bam_read)

	def close(self):
		self._samfile.close()


class ComparableAlignedSegment:
	"""
	Heapsort wants to be able to use the less than operator. Native
	AlignedSegment instances do not support this.
	"""
	def __init__(self, aligned_segment, source_id):
		self.segment = aligned_segment
		self.source_id = source_id

	def __lt__(self, other):
		self_pos = self.segment.reference_start
		other_pos = other.segment.reference_start
		return (self_pos < other_pos) or \
			(self_pos == other_pos and self.source_id < other.source_id)


class MultiBamReader:
	"""
	Read multiple sorted CRAM/BAM files and merge them on the fly.

	To avoid needing to handle renaming of duplicate read groups, this class
	just allows to specify a desired sample name. Doing that filtering here
	is much easier.
	"""

	def __init__(self, paths, *, reference=None):
		self._readers = []
		for source_id, path in enumerate(paths):
			self._readers.append(SampleBamReader(path, source_id=source_id, reference=reference))

	def fetch(self, reference=None, sample=None):
		"""
		Yield reads from the specified region in all the opened CRAM/BAM files,
		merging them on the fly. Each CRAM/BAM file must be indexed.

		Yields instances of AlignmentWithSourceID, where source_id corrsponds to 
		index of BAM file name given at construction time.

		If a sample name is given, only reads that belong to that sample are
		returned (the RG tags of each read and the RG header are used for that).
		"""
		assert reference is not None
		def make_comparable(reader):
			for alignment in reader.fetch(reference, sample):
				yield ComparableAlignedSegment(alignment.bam_alignment, alignment.source_id)
		iterators = []
		for reader in self._readers:
			if sample is None or reader.has_sample(sample):
				iterators.append(make_comparable(reader))
		if not iterators:
			raise SampleNotFoundError('Sample not found in any input CRAM/BAM file')
		for it in heapq.merge(*iterators):
			yield AlignmentWithSourceID(it.source_id, it.segment)

	def has_reference(self, name):
		return all(reader.has_reference(name) for reader in self._readers)

	def close(self):
		for f in self._readers:
			f.close()
