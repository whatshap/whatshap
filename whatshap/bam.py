import os
import pysam
import logging
import heapq
import gzip
from collections import defaultdict, namedtuple
import subprocess

logger = logging.getLogger(__name__)

AlignmentWithSourceID = namedtuple('AlignmentWithSourceID', ['source_id', 'bam_alignment'])


class BamIndexingError(Exception):
	pass


class SampleNotFoundError(Exception):
	pass


def index_bam(path):
	"""
	pysam.index fails silently on errors (such as when the input BAM file is not
	sorted). This function tries to always raise a BamIndexingError if something
	went wrong.
	"""
	# unfortunately, some versions of samtools index also fail silently on
	# corrupt BAM files. Before running it, we check the format manually.
	try:
		with gzip.GzipFile(path, 'rb') as bamfile:
			if bamfile.read(4) != b'BAM\1':
				raise BamIndexingError("{!r} is not a BAM file (header not found)".format(path))
	except OSError as e:
		raise BamIndexingError("{!r} is not a BAM file".format(path))
	po = subprocess.Popen(['samtools', 'index', path],
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	outs, errs = po.communicate()
	# samtools index also fails silently, at least in version 0.1.19 that comes
	# with Ubuntu 14.10, so to detect an error, we also inspect what it prints
	# on standard error.
	assert outs == ''
	if po.returncode != 0 or errs != '':
		raise BamIndexingError(errs)


class SampleBamReader:
	"""
	A wrapper for Samfile that provides only those reads from a BAM file that
	belong to a specified sample.
	"""
	def __init__(self, path, source_id=0):
		"""
		path -- path to BAM file
		"""
		# Raise an exception early if file does not exist or is not accessible.
		with open(path):
			pass
		bai1 = path + '.bai'
		bai2 = os.path.splitext(path)[0] + '.bai'
		if not os.path.exists(bai1) and not os.path.exists(bai2):
			logger.info('BAM index not found, creating it now.')
			index_bam(path)
		self.source_id = source_id
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

	def has_sample(self, sample):
		"""Return whether this file contains reads for the given sample"""
		return sample in self._sample_to_group_ids

	def fetch(self, reference, sample):
		"""
		Yield instances of AlignmentWithSourceID, with source_id value given
		at construction time.
		Raise KeyError if sample not found among samples named in RG header.
		"""
		if sample is None:
			for bam_read in self._samfile.fetch(reference):
				yield AlignmentWithSourceID(self.source_id, bam_read)
		else:
			try:
				read_groups = self._sample_to_group_ids[sample]
			except KeyError:
				raise SampleNotFoundError()
			for bam_read in self._samfile.fetch(reference):
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
		self_id = self.segment.reference_id
		self_pos = self.segment.reference_start
		other_id = other.segment.reference_id
		other_pos = other.segment.reference_start
		return self_id < other_id or (
			self_id == other_id and self_pos < other_pos) or (
			self_id == other_id and self_pos == other_pos and self.source_id < other.source_id)


class MultiBamReader:
	"""
	Read multiple sorted BAM files and merge them on the fly.

	To avoid needing to handle renaming of duplicate read groups, this class
	just allows to specify a desired sample name. Doing that filtering here
	is much easier.
	"""

	def __init__(self, paths):
		self._readers = []
		for source_id, path in enumerate(paths):
			self._readers.append(SampleBamReader(path, source_id))

	def fetch(self, reference=None, sample=None):
		"""
		Yield reads from the specified region in all the opened BAM files,
		merging them on the fly. Each BAM file must have a BAI index.

		Yields instances of AlignmentWithSourceID, where source_id corrsponds to 
		index of BAM file name given at construction time.

		If a sample name is given, only reads that belong to that sample are
		returned (the RG tags of each read and the RG header are used for that).
		"""
		def make_comparable(reader):
			for alignment in reader.fetch(reference, sample):
				yield ComparableAlignedSegment(alignment.bam_alignment, alignment.source_id)
		iterators = []
		for reader in self._readers:
			if reader.has_sample(sample):
				iterators.append(make_comparable(reader))
		if not iterators:
			raise SampleNotFoundError('Sample not found in any input BAM file')
		for it in heapq.merge(*iterators):
			yield AlignmentWithSourceID(it.source_id, it.segment)

	def close(self):
		for f in self._readers:
			f.close()


class HaplotypeBamWriter:
	'''Write haplotype specific BAM files of reads used for phasing.'''

	def __init__(self, input_bams, output_prefix, sample):
		# TODO: Take care of header: add all read groups and double check SQs, etc.
		assert len(input_bams) > 0
		template = pysam.AlignmentFile(input_bams[0], 'rb')
		self._outputfiles = [pysam.AlignmentFile('{}.{}.bam'.format(output_prefix,i), 'wb', template=template) for i in [1,2]]
		self.sample = sample
		if len(input_bams) == 1:
			self._reader = SampleBamReader(input_bams[0])
		else:
			logger.warning('Writing haplotype-specific BAM files for multiple input BAM files but proper resolution of name clashes not yet implemented.')
			self._reader = MultiBamReader(input_bams)

	def write(self, readset, bipartition, chromosome):
		'''Write all reads in the given readset to the two output files according to the 
		given bipartition. Expect reads to be from the given chromosome.'''
		assert len(readset) == len(bipartition)
		read_indices = { (r.source_id,r.name): i for i, r in enumerate(readset) }
		for alignment in self._reader.fetch(reference=chromosome, sample=self.sample):
			if alignment.bam_alignment.flag & 2048 != 0:
				continue
			if alignment.bam_alignment.is_secondary:
				continue
			read_id = (alignment.source_id, alignment.bam_alignment.qname)
			if read_id not in read_indices:
				continue
			target = bipartition[read_indices[read_id]]
			assert target in [0,1]
			# TODO: for multiple BAM files, translate read name and read group IDs
			self._outputfiles[target].write(alignment.bam_alignment)
