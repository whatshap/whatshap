import os
import pysam
import logging
import heapq
from collections import defaultdict

logger = logging.getLogger(__name__)

class SampleBamReader:
	"""
	Read only those reads from a BAM file that belong to a specified sample.
	"""
	def __init__(self, path):
		"""
		path -- path to BAM file
		"""
		bai1 = path + '.bai'
		bai2 = os.path.splitext(path)[0] + '.bai'
		if not os.path.exists(bai1) and not os.path.exists(bai2):
			logger.info('BAM index not found, creating it now.')
			pysam.index(path)
		self._samfile = pysam.Samfile(path)
		self._initialize_sample_to_group_ids()

	def _initialize_sample_to_group_ids(self):
		"""
		Create a dictionary that maps a sample name to a set of read group ids.
		"""
		read_groups = self._samfile.header.get('RG', [])  # a list of dicts
		logger.debug('Read groups in BAM header: %s', read_groups)
		samples = defaultdict(list)
		for read_group in read_groups:
			samples[read_group['SM']].append(read_group['ID'])
		self._sample_to_group_ids = {
			id: frozenset(values) for id, values in samples.items() }

	def fetch(self, reference, sample):
		"""
		Raise KeyError if sample not found among samples named in RG header.
		"""
		if sample is None:
			# PY32
			# Starting with Python 3.3, this loop could be replaced with
			# 'return self._samfile.fetch(reference)'
			for i in self._samfile.fetch(reference):
				yield i
			return
		read_groups = self._sample_to_group_ids[sample]
		for bam_read in self._samfile.fetch(reference):
			if bam_read.opt('RG') in read_groups:
				yield bam_read

	def close(self):
		self._samfile.close()


class ComparableAlignedSegment:
	"""
	Heapsort wants to be able to use the less than operator. Native
	AlignedSegment instances do not support this.
	"""
	def __init__(self, aligned_segment):
		self.segment = aligned_segment

	def __lt__(self, other):
		self_id = self.segment.reference_id
		self_pos = self.segment.reference_start
		other_id = other.segment.reference_id
		other_pos = other.segment.reference_start
		return self_id < other_id or (
			self_id == other_id and self_pos < other_pos)


class MultiBamReader:
	"""
	Read multiple sorted BAM files and merge them on the fly.

	To avoid needing to handle renaming of duplicate read groups, this class
	just allows to specify a desired sample name. Doing that filtering here
	is much easier.
	"""

	def __init__(self, paths):
		self._readers = []
		for path in paths:
			self._readers.append(SampleBamReader(path))

	def fetch(self, reference=None, sample=None):
		"""
		Yield reads from the specified region in all the opened BAM files,
		merging them on the fly. Each BAM file must have a BAI index.

		If a sample name is given, only reads that belong to that sample are
		returned (the RG tags of each read and the RG header are used for that).
		"""
		def make_comparable(reader):
			for segment in reader.fetch(reference, sample):
				yield ComparableAlignedSegment(segment)
		iterators = []
		for r in self._readers:
			iterators.append(make_comparable(r))
		for it in heapq.merge(*iterators):
			yield it.segment

	def close(self):
		for f in self._readers:
			f.close()


if __name__ == '__main__':
	import sys
	# merge given BAM files and write them into out.bam
	mb = MultiBam(sys.argv[1:])
	out = Samfile('out.bam', 'wb', template=mb._files[0])
	for r in mb.fetch():
		out.write(r)
	mb.close()
	out.close()
