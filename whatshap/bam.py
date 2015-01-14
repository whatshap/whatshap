import os
import pysam
import logging
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
		read_groups = self._samfile.header['RG']  # a list of dicts
		logger.debug('Read groups in BAM header: %s', read_groups)
		samples = defaultdict(list)
		for read_group in read_groups:
			samples[read_group['SM']].append(read_group['ID'])
		self._sample_to_group_ids = {
			id: frozenset(values) for id, values in samples.items() }

	def fetch(self, reference, sample):
		if sample is None:
			return self._file.fetch(reference)
		read_groups = self._sample_to_group_ids[sample]
		for bam_read in self._samfile.fetch(reference):
			if bam_read.opt('RG') in read_groups:
				yield bam_read

	def close(self):
		self._samfile.close()
