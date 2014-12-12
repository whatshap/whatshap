"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
import vcf
from . import __version__

logger = logging.getLogger(__name__)

#VcfVariant = namedtuple('VcfVariant', 'position reference_allele alternative_allele')

class VcfVariant:
	"""A variant in a VCF file"""
	def __init__(self, position, reference_allele, alternative_allele):
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele


class SampleNotFoundError(Exception):
	pass


def vcf_sample_reader(path, sample=None):
	"""
	Read only a single sample from a VCF.
	If sample is None, the first sample is used.

	Yield tuples (record, call).
	"""
	vcf_reader = vcf.Reader(filename=path)
	samples = vcf_reader.samples
	logger.info("Found %d samples in the VCF file.", len(samples))
	if sample is None:
		sample = samples[0]
		sample_index = 0
		if len(samples) > 1:
			logger.warn("More than one sample found in the VCF file, will work "
				"only on the first one (%s).", sample)
	else:
		try:
			sample_index = samples.index(sample)
		except ValueError:
			logger.error("Requested sample %r not found in VCF.", sample)
			raise SampleNotFoundError()
	for record in vcf_reader:
		call = record.samples[sample_index]
		yield sample, record, call


def parse_vcf(path, sample=None):
	"""
	Read a VCF and yield tuples (sample, chromosome, variants) for each
	chromosome for which there are variants in the VCF. chromosome is
	the name of the chromosome, and variants is a list of VcfVariant objects that
	represent the variants.

	path -- Path to VCF file
	sample -- The name of the sample whose calls should be extracted. If
		set to None, calls of the first sample are extracted.
	"""

	variants = []
	index = -1
	indices = None
	prev_chromosome = None
	for sample, record, call in vcf_sample_reader(path, sample):
		if record.CHROM != prev_chromosome:
			if prev_chromosome is not None:
				yield (sample, prev_chromosome, variants)
			prev_chromosome = record.CHROM
			variants = []
		if not record.is_snp:
			continue
		logger.debug("Call %s:%d %s→%s (Alleles: %s, %s; Het: %s; gt_bases; %s)",
			record.CHROM, record.start + 1,
			record.REF, record.ALT,
			record.alleles, call.gt_alleles, call.is_het, call.gt_bases)
		if len(record.ALT) != 1:
			logger.warn("Reading VCFs with multiple ALTs not implemented.")
			continue
		if not call.is_het:
			continue
		# found a heterozygous variant for the sample
		v = VcfVariant(
			position=record.start,
			reference_allele=record.REF,
			alternative_allele=record.ALT[0])
		variants.append(v)
	if prev_chromosome is not None:
		yield (sample, prev_chromosome, variants)


class PhasedVcfWriter:
	"""
	Read in a VCF file and write it back out with added phasing information.
	Phasing is written into HP and PQ tags, compatible with GATK’s
	ReadBackedPhasing.

	Avoid reading in full chromosomes as that uses too much memory for
	multi-sample VCFs.
	"""
	def __init__(self, in_path, command_line, out_file=sys.stdout):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry.
		out_file -- File-like object to which VCF is written.
		"""
		self._reader = vcf.Reader(filename=in_path)
		# FreeBayes adds phasing=none to its VCF output - remove that.
		self._reader.metadata['phasing'] = []
		if 'commandline' not in self._reader.metadata:
			self._reader.metadata['commandline'] = []
		self._reader.metadata['commandline'].append('"(whatshap ' + __version__ + ') ' + command_line + '"')
		self._reader.formats['HP'] = vcf.parser._Format(id='HP', num=None, type='String', desc='Phasing haplotype identifier')
		# TODO
		self._reader.formats['PQ'] = vcf.parser._Format(id='PQ', num=1, type='Float', desc='Phasing quality')

		self._writer = vcf.Writer(out_file, template=self._reader)
		logger.debug('Formats: %s', self._reader.formats)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)

	def _format_phasing_info(self, component, phase):
		assert phase in '01'
		phase = int(phase)
		return '{}-{},{}-{}'.format(component + 1, phase + 1, component + 1, 2 - phase)

	def write(self, chromosome, sample, superreads, components):
		"""
		Add phasing information to all variants on a single chromosome of a
		sample.
		"""
		assert self._unprocessed_record is None or (self._unprocessed_record.CHROM == chromosome)

		# TODO move sample parameter into the constructor
		sample_index = self._reader.samples.index(sample)

		# TODO don’t use dicts for *everything* ...
		phases = { v.position: v.allele for v in superreads[0].variants if v.allele in '01' }
		if self._unprocessed_record is not None:
			records_iter = itertools.chain([self._unprocessed_record], self._reader_iter)
		else:
			records_iter = self._reader_iter
		n = 0
		while True:
			try:
				record = next(records_iter)
			except StopIteration:
				break
			n += 1
			if record.CHROM != chromosome:
				# save it for later
				self._unprocessed_record = record
				assert n != 1
				break
			if record.start not in components:
				# Phasing info not available, just copy record
				self._writer.write_record(record)
				continue
			# Current PyVCF does not make it very easy to modify records/calls.
			record.add_format('HP')
			if record.FORMAT not in self._reader._format_cache:
				self._reader._format_cache[record.FORMAT] = self._reader._parse_sample_format(record.FORMAT)
			samp_fmt = self._reader._format_cache[record.FORMAT]

			# Set HP tag for all samples
			for i, call in enumerate(record.samples):
				if i == sample_index:
					phasing_info = self._format_phasing_info(components[record.start], phases[record.start])
				else:
					phasing_info = None
				call.data = samp_fmt(*(call.data + (phasing_info,)))
			self._writer.write_record(record)
