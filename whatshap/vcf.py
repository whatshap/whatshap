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

	def __str__(self):
		return "VcfVariant(pos={}, ref={}, alt={})".format(self.position+1,
			self.reference_allele, self.alternative_allele)


class SampleNotFoundError(Exception):
	pass


def vcf_sample_reader(path, sample=None):
	"""
	Read only a single sample from a VCF.
	If sample is None, the first sample is used.

	Yield tuples (sample, record, call).
	"""
	vcf_reader = vcf.Reader(filename=path)
	samples = vcf_reader.samples
	logger.info("Found %d sample(s) in the VCF file.", len(samples))
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
	assert sample is not None
	for record in vcf_reader:
		call = record.samples[sample_index]
		yield sample, record, call


def parse_vcf(path, indels=False, sample=None):
	"""
	Read a VCF and yield tuples (sample, chromosome, variants) for each
	chromosome for which there are variants in the VCF. chromosome is
	the name of the chromosome, and variants is a list of VcfVariant objects that
	represent the variants.

	path -- Path to VCF file
	indels -- Whether to include also insertions and deletions in the list of
		variants. Include only SNPs if set to False.
	sample -- The name of the sample whose calls should be extracted. If
		set to None, calls of the first sample are extracted.
	"""
	variants = []
	index = -1
	indices = None
	prev_chromosome = None
	n_indels = 0
	n_snps = 0
	for sample, record, call in vcf_sample_reader(path, sample):
		if record.CHROM != prev_chromosome:
			if prev_chromosome is not None:
				yield (sample, prev_chromosome, variants)
			prev_chromosome = record.CHROM
			variants = []
		alleles = [ str(record.alleles[int(s)]) for s in call.gt_alleles ]
		"""
		logger.debug("Call %s:%d %s→%s (Alleles: %s)",
			record.CHROM, record.start + 1,
			record.REF, record.ALT,
			alleles)
		"""
		if not call.is_het:
			continue
		assert len(alleles) == 2

		# Normalize variants in which the first two bases are identical,
		# such as CTG -> CTAAA (which is changed to TG -> TAAA).
		a0, a1 = alleles[0:2]
		pos = record.start
		while len(a0) >= 2 and len(a1) >= 2 and a0[0:2] == a1[0:2]:
			a0, a1 = a0[1:], a1[1:]
			pos += 1
		assert a0 != a1
		v = VcfVariant(position=pos, reference_allele=a0, alternative_allele=a1)
		if len(a0) == 1 and len(a1) == 1:
			n_snps += 1
			variants.append(v)
			continue
		if not indels:
			continue

		if a0[0] == a1[0] and ((len(a0) == 1) != (len(a1) == 1)):
			n_indels += 1
			variants.append(v)
			continue
		# Something like GCG -> TCT or CTCTC -> CA occurred.
		logger.warn('Complex variant found, skipping: %s', v)
	logger.debug("No. of SNPs on this chromosome: %s; no. of indels: %s", n_snps, n_indels)
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
		self._reader.formats['PQ'] = vcf.parser._Format(id='PQ', num=1, type='Float', desc='Phasing quality')

		self._writer = vcf.Writer(out_file, template=self._reader)
		logger.debug('Formats: %s', self._reader.formats)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)
		self._hp_found_warned = False

	def _format_phasing_info(self, component, phase):
		assert phase in [0,1]
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
		phases = { position: allele for position, base, allele, quality in superreads[0] if allele in [0,1] }
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
			if 'HP' not in record.FORMAT.split(':'):
				record.add_format('HP')
				if record.FORMAT not in self._reader._format_cache:
					self._reader._format_cache[record.FORMAT] = self._reader._parse_sample_format(record.FORMAT)
			samp_fmt = self._reader._format_cache[record.FORMAT]
			# Set HP tag for all samples
			for i, call in enumerate(record.samples):
				if i == sample_index:
					if (hasattr(call.data, 'HP') and call.data.HP is not None
							and not self._hp_found_warned):
						logger.warn('Ignoring existing phasing information '
							'found in input VCF (HP tag exists).')
						self._hp_found_warned = True
					# Set or overwrite HP tag
					phasing_info = self._format_phasing_info(components[record.start], phases[record.start])
					values = vars(call.data)
					values['HP'] = phasing_info
					call.data = samp_fmt(**values)
				elif not hasattr(call.data, 'HP'):
					# HP tag missing, set it to "."
					values = vars(call.data)
					values['HP'] = None
					call.data = samp_fmt(**values)
			self._writer.write_record(record)
