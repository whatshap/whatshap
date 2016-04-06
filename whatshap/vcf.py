"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
from array import array
from collections import defaultdict
import vcf

logger = logging.getLogger(__name__)


class VcfVariant:
	"""A variant in a VCF file"""
	def __init__(self, position, reference_allele, alternative_allele):
		"""
		position -- 0-based start coordinate
		reference_allele -- string
		alternative_allele -- string
		genotype -- equal to sum of genotypes:
			0 is 0/0 (homozygous reference)
			1 is 0/1 or 1/0 (heterozygous)
			2 is 1/1 (homozygous alternative)

		Multi-ALT sites are not modelled.
		"""
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele

	def __repr__(self):
		return "VcfVariant(pos={}, ref={}, alt={})".format(self.position+1,
			self.reference_allele, self.alternative_allele)


class VariantTable:
	"""
	For a single chromosome, store variants and their genotypes.

	Each column contains the genotypes of a single sample.

	chromosome -- chromosome name
	samples -- list of sample names

	TODO We are re-implementing a pandas.DataFrame here.
	"""
	def __init__(self, chromosome, samples):
		self.chromosome = chromosome
		self.samples = samples
		self.genotypes = [ array('b', []) for _ in samples ]
		self.variants = []
		self._sample_to_index = { sample: index for index, sample in enumerate(samples) }

	#def add_sample(self, name, genotypes):
		#"Add a column to the table"
		#if len(genotypes) != len(self.variants):
			#raise ValueError('Expecting as many genotypes as there are variants')
		#self._name_to_index[name] = len(self.samples)
		#self.samples.append(name)
		#self.genotypes.append(genotypes)

	def add_variant(self, variant, genotypes):
		if len(genotypes) != len(self.genotypes):
			raise ValueError('Expecting as many genotypes as there are samples')
		self.variants.append(variant)
		for i, genotype in enumerate(genotypes):
			self.genotypes[i].append(genotype)

	def genotypes_of(self, sample):
		"""Retrieve genotypes by sample name"""
		return self.genotypes[self._sample_to_index[sample]]

	def id_of(self, sample):
		"""Return a unique int id of a sample given by name"""
		return self._sample_to_index[sample]

	def remove_rows_by_index(self, indices):
		"""Remove variants given by their index in the variant list"""
		for i in sorted(indices, reverse=True):
			del self.variants[i]
			for gt in self.genotypes:
				del gt[i]

		for gt in self.genotypes:
			assert len(self.variants) == len(gt)
		assert len(self.samples) == len(self.genotypes)

	def subset_rows_by_position(self, positions):
		"""Keep only rows given in positions, discard the rest"""
		positions = frozenset(positions)
		to_discard = [ i for i, v in enumerate(self.variants) if v.position not in positions ]
		self.remove_rows_by_index(to_discard)


class SampleNotFoundError(Exception):
	pass


class VcfReader:
	"""
	Read a VCF file chromosome by chromosome.

	TODO
	- Sites with multiple ALTs are skipped.
	- Always include deletions since they can 'overlap' other variants
	"""
	def __init__(self, path, indels=False, samples=None):
		"""
		path -- Path to VCF file
		samples -- Names of samples to extract. If None, extract all.
		indels -- Whether to include also insertions and deletions in the list of
			variants.
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
		Yield VariantTable objects for each chromosome.

		Indels are normalized; multi-ALT sites and complex variants are skipped.
		"""
		for chromosome, records in self._group_by_chromosome():
			yield self._process_single_chromosome(chromosome, records)

	def _process_single_chromosome(self, chromosome, records):
		n_snps = 0
		n_indels = 0
		n_complex = 0
		n_multi = 0
		table = VariantTable(chromosome, self.samples)
		for record in records:
			if len(record.ALT) > 1:
				# Skip sites with multiple alternative alleles (see note in docstring)
				n_multi += 1
				continue

			ref, alt = str(record.REF), str(record.ALT[0])
			pos, ref, alt = self._normalize(record.start, ref, alt)

			# PyVCF pecularity: gt_alleles is a list of the alleles in the
			# GT field, but as strings.
			# For example, when GT is 0/1, gt_alleles is ['0', '1'].
			# And when GT is 2|1, gt_alleles is ['2', '1'].
			"""
			logger.debug("Call %s:%d %s→%s",
				record.CHROM, record.start + 1, record.REF, record.ALT)
			"""

			# Determine variant type: SNP, indel or complex
			if len(ref) == len(alt) == 1:
				n_snps += 1
				variant = VcfVariant(position=pos, reference_allele=ref, alternative_allele=alt)
			elif ref[0] == alt[0] and ((len(ref) == 1) != (len(alt) == 1)):
				n_indels += 1
				if not self._indels:
					continue
				variant = VcfVariant(position=pos+1, reference_allele=ref[1:], alternative_allele=alt[1:])
			else:
				# A complex variant such as GCG -> TCT or CTCTC -> CA occurred.
				n_complex += 1
				continue

			genotypes = [ call.gt_type for call in record.samples ]
			table.add_variant(variant, genotypes)

		logger.debug("No. of SNPs on this chromosome: %s; no. of indels: %s. "
			"Skipped %s complex variants. Skipped %s multi-ALTs.", n_snps,
			n_indels, n_complex, n_multi)

		# TODO remove overlapping variants
		return table


def remove_overlapping_calls(calls):
	"""
	Filter a list of variants such that no variants overlap each other.
	This applies mainly to deletions: If they occur too close to another
	variant, the deletion and the other variant are removed.

	This function also guarantees that the positions of the returned variants
	are unique. For that, it may also remove other variants (not necessarily
	involved in a deletion).

	calls -- a list of VcfVariant objects

	Return a list of VcfVariant objects.
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
