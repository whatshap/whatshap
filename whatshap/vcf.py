"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
from collections import defaultdict
import vcf

logger = logging.getLogger(__name__)


# TODO it is not quite accurate to call the attributes reference_allele and
# alternative_allele because there can be non-reference heterozygous variants
# For example, if the VCF file has an entry such as A -> C,T and the GT is 1/2,
# then the alleles are C and T, and both of them are 'alternative'.
class VcfCall:
	"""A called variant in a VCF file"""
	def __init__(self, position, reference_allele, alternative_allele, genotype):
		"""
		position -- 0-based start coordinate
		reference_allele -- string
		alternative_allele -- string
		genotype -- as in PyVCF:
			homozyguous reference is 0,
			heterozygous is 1,
			homozygous alternative is 2

		TODO seems unnecessary to store the genotype as an extra attribute
			because it is redundant with the alleles
		TODO reference_allele and alternative_allele are not the best names:
			it should be allele1 and allele2 (or a list 'alleles') because there
			can be non-reference heterozygous variants (such as C -> G,A with GT 1/2)
		"""
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele
		self.genotype = genotype

	def __repr__(self):
		return "VcfCall(pos={}, ref={}, alt={}, genotype={})".format(self.position+1,
			self.reference_allele, self.alternative_allele, self.genotype)

	def is_heterozygous(self):
		return self.reference_allele != self.alternative_allele


class SampleNotFoundError(Exception):
	pass


class VcfReader:
	"""
	Read a VCF file chromosome by chromosome.

	sample -- The name of the sample whose calls should be extracted. If
		set to None, calls of the first sample are extracted.

	TODO
	- Sites with multiple ALTs are currently skipped. This is easy to fix, but someone
	needs to understand what the expected output in the HP tag should then be.
	"""
	def __init__(self, path, indels=False, samples=None):
		"""
		path -- Path to VCF file

		indels -- Whether to include also insertions and deletions in the list of
		variants. Include only SNPs if set to False.
		TODO this should always be enabled since deletions can 'overlap' other variants
		"""
		self._indels = indels
		self._vcf_reader = vcf.Reader(filename=path)
		self.samples = self._vcf_reader.samples  # intentionally public

		if samples is None:
			self._samples_of_interest = frozenset(self.samples)
		else:
			for sample in samples:
				if sample not in self.samples:
					raise SampleNotFoundError("Requested sample %r not found in VCF.", sample)
			self._samples_of_interest = frozenset(samples)
		logger.debug("Found %d sample(s) in the VCF file.", len(self.samples))

	def _group_by_chromosome(self):
		"""
		Yield (chromosome, records) tuples, where records is a list of the
		VCF records on that chromosome.
		"""
		records = []
		prev_chromosome = None
		for record in self._vcf_reader:
			if record.CHROM != prev_chromosome:
				if prev_chromosome is not None:
					yield (prev_chromosome, records)
				prev_chromosome = record.CHROM
				records = []
			records.append(record)
		if records:
			yield (prev_chromosome, records)

	def _normalize(self, pos, ref, alt):
		"""
		Normalize variants that share a common prefix. The first shared base is
		kept. For example, GCTG -> GCTAAA is changed to TG -> TAAA.

		Return a (pos, ref, alt) tuple.
		"""
		while len(ref) >= 2 and len(alt) >= 2 and ref[0:2] == alt[0:2]:
			ref, alt = ref[1:], alt[1:]
			pos += 1
		return pos, ref, alt

	def __iter__(self):
		"""
		Yield tuples (chromosome, calls), where
		- chromosome is the chromosome name
		- calls is a dictionary that maps sample names to a list of VcfCall
		objects. Indels are normalized, (even the first shared base is removed),
		multi-ALT sites are skipped, and also complex variants.
		"""
		for chromosome, records in self._group_by_chromosome():
			yield (chromosome, self._process_single_chromosome(records))

	def _process_single_chromosome(self, records):
		n_snps = 0
		n_indels = 0
		n_complex = 0
		n_multi = 0
		samples = defaultdict(list)
		for record in records:
			if len(record.ALT) > 1:
				# Skip sites with multiple alternative alleles (see note in docstring)
				n_multi += 1
				continue

			for sample_index, call in enumerate(record.samples):
				sample = self.samples[sample_index]

				# Two PyVCF pecularities:
				#
				# gt_alleles is a list of the alleles in the GT field, but
				# as strings.
				# For example, when GT is 0/1, gt_alleles is ['0', '1'].
				# And when GT is 2|1, gt_alleles is ['2', '1'].
				#
				# record.alleles is a list where the first item is a string, but
				# the other objects are not (they are Substitution objects etc.) -
				# so do not remove the str() from the line below.
				#
				# TODO The sorted() should not be necessary
				alleles = [ str(record.alleles[int(s)]) for s in sorted(call.gt_alleles) ]
				"""
				logger.debug("Call %s:%d %s→%s with alleles %s, genotype %s",
					record.CHROM, record.start + 1,
					record.REF, record.ALT,
					alleles, call.gt_type)
				"""
				assert len(alleles) == 2
				pos, ref, alt = self._normalize(record.start, alleles[0], alleles[1])

				# Determine variant type: SNP, indel or complex
				if len(ref) == 1 and len(alt) == 1:
					n_snps += 1
					v = VcfCall(position=pos, reference_allele=ref, alternative_allele=alt, genotype=call.gt_type)
				elif ref[0] == alt[0] and ((len(ref) == 1) != (len(alt) == 1)):
					n_indels += 1
					if not self._indels:
						continue
					v = VcfCall(position=pos+1, reference_allele=ref[1:], alternative_allele=alt[1:], genotype=call.gt_type)
				else:
					# A complex variant such as GCG -> TCT or CTCTC -> CA occurred.
					n_complex += 1
					# TODO deal with complex variants
					# v = VcfCall(position=pos, reference_allele=a0, alternative_allele=a1, genotype=call.gt_type)
					continue
				samples[sample].append(v)

		logger.debug("No. of SNPs on this chromosome: %s; no. of indels: %s. "
			"Skipped %s complex variants. Skipped %s multi-ALTs.", n_snps, n_indels, n_complex, n_multi)

		nonoverlapping = dict()
		for sample, calls in samples.items():
			if sample in self._samples_of_interest:
				nonoverlapping[sample] = remove_overlapping_calls(calls)
		return nonoverlapping


def remove_overlapping_calls(calls):
	"""
	Filter a list of variants such that no variants overlap each other.
	This applies mainly to deletions: If they occur too close to another
	variant, the deletion and the other variant are removed.

	This function also guarantees that the positions of the returned variants
	are unique. For that, it may also remove other variants (not necessarily
	involved in a deletion).

	calls -- a list of VcfCall objects

	Return a list of VcfCall objects.
	"""
	# TODO obviously, this is not implemented ...
	return calls


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
		command_line = command_line.replace('"', '')
		self._reader.metadata['commandline'].append('"' + command_line + '"')
		self._reader.formats['HP'] = vcf.parser._Format(id='HP', num=None, type='String', desc='Phasing haplotype identifier')
		self._reader.formats['PQ'] = vcf.parser._Format(id='PQ', num=1, type='Float', desc='Phasing quality')

		self._writer = vcf.Writer(out_file, template=self._reader)
		logger.debug('Formats: %s', self._reader.formats)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)
		self._hp_found_warned = False

	def _format_phasing_info(self, component, phase):
		"""
		component -- name of the component
		phase -- 0 or 1
		"""
		assert phase in [0,1]
		return '{}-{},{}-{}'.format(component + 1, phase + 1, component + 1, 2 - phase)

	def write(self, chromosome, sample_superreads, sample_components):
		"""
		Add phasing information to all variants on a single chromosome.

		chromosome -- name of chromosome
		sample_superreads -- dictionary that maps each sample to a superread
		sample_components -- a dictionary that maps each sample to its connected components

			Each component in turn is a dict that maps each variant position to a
			component, where a component is identified by the position of its
			left-most variant

		# TODO introduce a PhasingResult class that combines superread and component
		"""
		assert self._unprocessed_record is None or (self._unprocessed_record.CHROM == chromosome)

		if self._unprocessed_record is not None:
			records_iter = itertools.chain([self._unprocessed_record], self._reader_iter)
		else:
			records_iter = self._reader_iter
		allowed_alleles = frozenset({(0, 1), (1, 0), (0, 0), (1, 1)})
		sample_phases = dict()
		for sample, superreads in sample_superreads.items():
			sample_phases[sample] = {
				v1.position: v1.allele for v1, v2 in zip(*superreads)
					if (v1.allele, v2.allele) in allowed_alleles
			}
		n = 0
		for record in records_iter:
			n += 1
			if record.CHROM != chromosome:
				# save it for later
				self._unprocessed_record = record
				assert n != 1
				break

			# Determine whether the variant is phased in any sample
			is_phased = True
			for sample in self._reader.samples:
				if sample in sample_superreads:
					components = sample_components[sample]
					phases = sample_phases[sample]
					if record.start in components and record.start in phases:
						break
			else:
				is_phased = False

			if is_phased:
				# Current PyVCF does not make it very easy to modify records/calls.
				# Add HP tag to FORMAT
				if 'HP' not in record.FORMAT.split(':'):
					record.add_format('HP')
					if record.FORMAT not in self._reader._format_cache:
						self._reader._format_cache[record.FORMAT] = self._reader._parse_sample_format(record.FORMAT)
				samp_fmt = self._reader._format_cache[record.FORMAT]

				# Set HP tag for all samples
				for i, call in enumerate(record.samples):
					sample = self._reader.samples[i]
					if sample not in sample_superreads:
						# This sample has not been phased, leave unchanged
						continue
					components = sample_components[sample]
					phases = sample_phases[sample]

					if (hasattr(call.data, 'HP') and call.data.HP is not None
							and not self._hp_found_warned):
						logger.warning('Ignoring existing phasing information '
							'found in input VCF (HP tag exists).')
						self._hp_found_warned = True

					values = call.data._asdict()
					if record.start in components and record.start in phases and call.is_het:
						# Set or overwrite HP tag
						values['HP'] = self._format_phasing_info(
							components[record.start], phases[record.start])
					else:
						# Unphased - set HP to '.'
						values['HP'] = None
					call.data = samp_fmt(**values)
			self._writer.write_record(record)
