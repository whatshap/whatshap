"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
import math
from array import array
from collections import namedtuple, defaultdict
import vcf
from .core import Read, PhredGenotypeLikelihoods

logger = logging.getLogger(__name__)


class VcfNotSortedError(Exception):
	pass


class VcfVariant:
	"""A variant in a VCF file (not to be confused with core.Variant)"""

	__slots__ = ('position', 'reference_allele', 'alternative_allele')

	def __init__(self, position, reference_allele, alternative_allele):
		"""
		position -- 0-based start coordinate
		reference_allele -- string
		alternative_allele -- string

		Multi-ALT sites are not modelled.
		"""
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele

	def __repr__(self):
		return "VcfVariant({}, {!r}, {!r})".format(self.position,
			self.reference_allele, self.alternative_allele)

	def __hash__(self):
		return hash((self.position, self.reference_allele, self.alternative_allele))

	def __eq__(self, other):
		return (self.position == other.position) and \
		       (self.reference_allele == other.reference_allele) and \
		       (self.alternative_allele == other.alternative_allele)

	def __lt__(self, other):
		return (self.position, self.reference_allele, self.alternative_allele) < (other.position, other.reference_allele, other.alternative_allele)

	def is_snv(self):
		return (self.reference_allele != self.alternative_allele) and (
			len(self.reference_allele) == len(self.alternative_allele) == 1)

	def normalized(self):
		"""
		Return a normalized version of this variant.

		Common prefixes and/or suffixes between the reference and alternative allele are removed,
		and the position is adjusted as necessary.

		>>> VcfVariant(100, 'GCTGTT', 'GCTAAATT').normalized()
		VcfVariant(103, 'G', 'AAA')
		"""
		pos, ref, alt = self.position, self.reference_allele, self.alternative_allele
		while len(ref) >= 1 and len(alt) >= 1 and ref[-1] == alt[-1]:
			ref, alt = ref[:-1], alt[:-1]

		while len(ref) >= 1 and len(alt) >= 1 and ref[0] == alt[0]:
			ref, alt = ref[1:], alt[1:]
			pos += 1

		return VcfVariant(pos, ref, alt)


class GenotypeLikelihoods:
	__slots__ = ('log_prob_g0', 'log_prob_g1', 'log_prob_g2')

	def __init__(self, log_prob_g0, log_prob_g1, log_prob_g2):
		"""Likelihoods of the three genotypes 0, 1, 2 to be given
		as log10 of the original probability."""
		self.log_prob_g0 = log_prob_g0
		self.log_prob_g1 = log_prob_g1
		self.log_prob_g2 = log_prob_g2

	def __repr__(self):
		return "GenotypeLikelihoods({}, {}, {})".format(self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)

	def log10_probs(self):
		return ( self.log_prob_g0, self.log_prob_g1, self.log_prob_g2 )

	def log10_prob_of(self, genotype):
		return self.log10_probs()[genotype]

	def as_phred(self, regularizer=None):
		if regularizer is None:
			# shift log likelihoods such that the largest one is zero
			m = max(self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)
			return PhredGenotypeLikelihoods(
				round((self.log_prob_g0-m) * -10),
				round((self.log_prob_g1-m) * -10),
				round((self.log_prob_g2-m) * -10)
			)
		else:
			p = [ 10**x for x in (self.log_prob_g0, self.log_prob_g1, self.log_prob_g2) ]
			s = sum(p)
			p = [ x/s + regularizer for x in p ]
			m = max(p)
			return PhredGenotypeLikelihoods( *(round(-10*math.log10(x/m)) for x in p) )


class VariantTable:
	"""
	For a single chromosome, store variants and their genotypes.
	Each row of this table contains a variant, each column
	contains the genotypes of a single sample.

	chromosome -- chromosome name
	samples -- list of sample names
	"""
	def __init__(self, chromosome, samples):
		self.chromosome = chromosome
		self.samples = samples
		self.genotypes = [array('b', []) for _ in samples]
		self.phases = [[] for _ in samples]
		self.genotype_likelihoods = [[] for _ in samples]
		self.variants = []
		self._sample_to_index = {sample: index for index, sample in enumerate(samples)}

	def __len__(self):
		return len(self.variants)

	#def add_sample(self, name, genotypes):
		#"Add a column to the table"
		#if len(genotypes) != len(self.variants):
			#raise ValueError('Expecting as many genotypes as there are variants')
		#self._name_to_index[name] = len(self.samples)
		#self.samples.append(name)
		#self.genotypes.append(genotypes)

	def add_variant(self, variant, genotypes, phases, genotype_likelihoods):
		"""
		Add a row to the table

		variant -- a VcfVariant
		genotypes -- iterable of ints that encode the genotypes of the samples:
			-1 represents an unknown genotype
			0 represents 0/0 (homozygous reference)
			1 represents 0/1 or 1/0 (heterozygous)
			2 represents 1/1 (homozygous alternative)
		phases -- iterable of VariantCallPhase objects
		genotype_likelihoods -- iterable of GenotypeLikelihoods objects
		"""
		if len(genotypes) != len(self.genotypes):
			raise ValueError('Expecting as many genotypes as there are samples')
		if len(phases) != len(self.phases):
			raise ValueError('Expecting as many phases as there are samples')
		self.variants.append(variant)
		for i, genotype in enumerate(genotypes):
			self.genotypes[i].append(genotype)
		for i, phase in enumerate(phases):
			self.phases[i].append(phase)
		for i, gl in enumerate(genotype_likelihoods):
			self.genotype_likelihoods[i].append(gl)

	def genotypes_of(self, sample):
		"""Retrieve genotypes by sample name"""
		return self.genotypes[self._sample_to_index[sample]]

	def set_genotypes_of(self, sample, genotypes):
		"""Set genotypes by sample name"""
		assert len(genotypes) == len(self.variants)
		self.genotypes[self._sample_to_index[sample]] = genotypes


	def genotype_likelihoods_of(self, sample):
		"""Retrieve genotype likelihoods by sample name"""
		return self.genotype_likelihoods[self._sample_to_index[sample]]

	def set_genotype_likelihoods_of(self, sample, genotype_likelihoods):
		"""Set genotype likelihoods by sample name"""
		assert len(genotype_likelihoods) == len(self.variants)
		self.genotype_likelihoods[self._sample_to_index[sample]] = genotype_likelihoods

	def phases_of(self, sample):
		"""Retrieve phases by sample name"""
		return self.phases[self._sample_to_index[sample]]

	def num_of_blocks_of(self, sample):
		""" Retrieve the number of blocks of the sample"""
		return len(set([i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None]))

	def id_of(self, sample):
		"""Return a unique int id of a sample given by name"""
		return self._sample_to_index[sample]

	def remove_rows_by_index(self, indices):
		"""Remove variants given by their index in the variant list"""
		for i in sorted(indices, reverse=True):
			del self.variants[i]
			for gt in self.genotypes:
				del gt[i]
			for ph in self.phases:
				del ph[i]
			for gl in self.genotype_likelihoods:
				del gl[i]

		for gt in self.genotypes:
			assert len(self.variants) == len(gt)
		for ph in self.phases:
			assert len(self.variants) == len(ph)
		for gl in self.genotype_likelihoods:
			assert len(self.variants) == len(gl)
		assert len(self.samples) == len(self.genotypes) == len(self.phases) == len(self.genotype_likelihoods)

	def subset_rows_by_position(self, positions):
		"""Keep only rows given in positions, discard the rest"""
		positions = frozenset(positions)
		to_discard = [ i for i, v in enumerate(self.variants) if v.position not in positions ]
		self.remove_rows_by_index(to_discard)

	def phased_blocks_as_reads(self, sample, input_variants, source_id, numeric_sample_id, default_quality=20, mapq=100):
		"""
		Yields one sorted core.Read object per phased block, encoding the phase information as
		if this block was a single sequencing read. Reads are yielded in arbitrary order.

		sample -- name of sample to retrieve
		input_variants -- variants of interest, i.e. only these variants will be retrieved
		source_id -- source_id to be assigned to each read
		numeric_sample_id -- sample id to be stored in generated reads
		default_quality -- quality assigned to heterozygous with missing phasing quality
		mapq -- mapping quality for generated reads
		"""
		try:
			sample_index = self._sample_to_index[sample]
		except KeyError:
			return
		input_variant_set = set(input_variants)
		read_map = {} # maps block_id core.Read objects
		assert len(self.variants) == len(self.genotypes[sample_index]) == len(self.phases[sample_index])
		for variant, genotype, phase in zip(self.variants, self.genotypes[sample_index], self.phases[sample_index]):
			if not variant in input_variant_set:
				continue
			if genotype != 1:
				continue
			if phase is None:
				continue
			if phase.quality is None:
				quality = default_quality
			else:
				quality = phase.quality
			if phase.block_id in read_map:
				read_map[phase.block_id].add_variant(variant.position, phase.phase, quality)
			else:
				r = Read('{}_block_{}'.format(sample, phase.block_id), mapq, source_id, numeric_sample_id)
				r.add_variant(variant.position, phase.phase, quality)
				read_map[phase.block_id] = r
		for key, read in read_map.items():
			read.sort()
			if len(read) > 1:
				yield read


class MixedPhasingError(Exception):
	pass


VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'phase', 'quality'])
VariantCallPhase.__doc__ = \
	"""
	block_id is a numeric id of the phased block and
	phase is either 0 or 1 (indicating whether the REF allele is on haplotype 0 or 1).
	"""


class VcfReader:
	"""
	Read a VCF file chromosome by chromosome.
	"""
	def __init__(self, path, indels=False, phases=False, genotype_likelihoods=False, ignore_genotypes=False):
		"""
		path -- Path to VCF file
		indels -- Whether to include also insertions and deletions in the list of
			variants.
		ignore_genotypes: in case of genotyping algorithm, no genotypes may be given in
								vcf, so ignore all genotypes
		"""
		# TODO Always include deletions since they can 'overlap' other variants
		self._indels = indels
		self._vcf_reader = vcf.Reader(filename=path)
		self._phases = phases
		self._genotype_likelihoods = genotype_likelihoods
		self.samples = self._vcf_reader.samples  # intentionally public
		self.ignore_genotypes = ignore_genotypes
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

	def __iter__(self):
		"""
		Yield VariantTable objects for each chromosome.

		Multi-ALT sites are skipped.
		"""
		for chromosome, records in self._group_by_chromosome():
			yield self._process_single_chromosome(chromosome, records)

	@staticmethod
	def _extract_HP_phase(call):
		HP = getattr(call.data, 'HP', None)
		if HP is None:
			return None
		assert len(HP) == 2
		fields = [[int(x) for x in s.split('-')] for s in HP]
		assert fields[0][0] == fields[1][0]
		block_id = fields[0][0]
		phase1, phase2 = fields[0][1]-1, fields[1][1]-1
		assert ((phase1, phase2) == (0, 1)) or ((phase1, phase2) == (1, 0))
		return VariantCallPhase(block_id=block_id, phase=phase1, quality=getattr(call.data, 'PQ', None))

	@staticmethod
	def _extract_GT_PS_phase(call):
		if not call.is_het:
			return None
		if not call.phased:
			return None
		block_id = getattr(call.data, 'PS', 0)
		if block_id is None:
			block_id = 0
		assert call.data.GT in ['0|1','1|0']
		phase = int(call.data.GT[0])
		return VariantCallPhase(block_id=block_id, phase=phase, quality=getattr(call.data, 'PQ', None))

	def _process_single_chromosome(self, chromosome, records):
		phase_detected = None
		n_snvs = 0
		n_other = 0
		n_multi = 0
		table = VariantTable(chromosome, self.samples)
		prev_position = None
		for record in records:
			if len(record.ALT) > 1:
				# Multi-ALT sites are not supported, yet
				n_multi += 1
				continue

			pos, ref, alt = record.start, str(record.REF), str(record.ALT[0])
			if len(ref) == len(alt) == 1:
				n_snvs += 1
			else:
				n_other += 1
				if not self._indels:
					continue

			if (prev_position is not None) and (prev_position > pos):
				raise VcfNotSortedError('VCF not ordered: {}:{} appears before {}:{}'.format(chromosome, prev_position+1, chromosome, pos+1))

			if prev_position == pos:
				logger.warning('Skipping duplicated position %s on chromosome %r', pos+1, chromosome)
				continue
			prev_position = pos

			# Read phasing information (allow GT/PS or HP phase information, but not both),
			# if requested
			if self._phases:
				phases = []
				for call in record.samples:
					phase = None
					for extract_phase, phase_name in [(self._extract_HP_phase, 'HP'), (self._extract_GT_PS_phase, 'GT_PS')]:
						p = extract_phase(call)
						if p is not None:
							if phase_detected is None:
								phase_detected = phase_name
							elif phase_detected != phase_name:
								raise MixedPhasingError('Mixed phasing information in input VCF (e.g. mixing PS and HP fields)')
							phase = p
					phases.append(phase)
			else:
				phases = [ None ] * len(record.samples)

			# Read genotype likelihoods, if requested
			if self._genotype_likelihoods:
				genotype_likelihoods = []
				for call in record.samples:
					GL = getattr(call.data, 'GL', None)
					PL = getattr(call.data, 'PL', None)
					# Prefer GLs (floats) over PLs (ints) if both should be present
					if GL is not None:
						assert len(GL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods(*GL))
					elif PL is not None:
						assert len(PL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods( *(pl/-10 for pl in PL) ))
					else:
						genotype_likelihoods.append(None)
			else:
				genotype_likelihoods = [ None ] * len(record.samples)

			# PyVCF pecularity: gt_alleles is a list of the alleles in the
			# GT field, but as strings.
			# For example, when GT is 0/1, gt_alleles is ['0', '1'].
			# And when GT is 2|1, gt_alleles is ['2', '1'].
			GT_TO_INT = { 0: 0, 1: 1, 2: 2, None: -1 }
			if not self.ignore_genotypes:
				genotypes = array('b', (GT_TO_INT[call.gt_type] for call in record.samples))
			else:
				genotypes = array('b', ([-1] * len(self.samples)))
				phases = [None] * len(self.samples)
			variant = VcfVariant(position=pos, reference_allele=ref, alternative_allele=alt)
			table.add_variant(variant, genotypes, phases, genotype_likelihoods)

		logger.debug("Parsed %s SNVs and %s non-SNVs. Also skipped %s multi-ALTs.", n_snvs,
			n_other, n_multi)

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

GenotypeChange = namedtuple('GenotypeChange', ['sample', 'chromosome', 'variant', 'old_gt', 'new_gt'])

class PhasedVcfWriter:
	"""
	Read in a VCF file and write it back out with added phasing information.
	Phasing is written into HP and PQ tags, compatible with GATK’s
	ReadBackedPhasing.

	Avoid reading in full chromosomes as that uses too much memory for
	multi-sample VCFs.
	"""
	def __init__(self, in_path, command_line, out_file=sys.stdout, tag='PS'):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry
			(use None to not add this to the VCF header)
		out_file -- Open file-like object to which VCF is written.
		tag -- which type of tag to write, either 'PS' or 'HP'
		"""
		self._reader = vcf.Reader(filename=in_path)
		# FreeBayes adds phasing=none to its VCF output - remove that.
		self._reader.metadata['phasing'] = []
		if command_line is not None:
			if 'commandline' not in self._reader.metadata:
				self._reader.metadata['commandline'] = []
			command_line = command_line.replace('"', '')
			self._reader.metadata['commandline'].append('"' + command_line + '"')
		if tag not in ('HP', 'PS'):
			raise ValueError('Tag must be either "HP" or "PS"')
		if tag == 'HP':
			desc = 'Phasing haplotype identifier'
			fmt = vcf.parser._Format(id=tag, num=None, type='String', desc=desc)
		else:
			desc = 'Phase set identifier'
			fmt = vcf.parser._Format(id=tag, num=1, type='Integer', desc=desc)
		self._reader.formats[tag] = fmt

		# self._reader.formats['PQ'] = vcf.parser._Format(id='PQ', num=1, type='Float', desc='Phasing quality')

		self._writer = vcf.Writer(out_file, template=self._reader)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)
		self._phase_tag_found_warned = False
		self.tag = tag
		self._set_phasing_tags = self._set_HP if tag == 'HP' else self._set_PS

	@property
	def samples(self):
		return self._reader.samples

	def _set_HP(self, values, component, phase):
		"""
		values -- tag dict to update
		component -- name of the component
		phase -- 0 or 1
		"""
		assert phase in [0, 1]
		values['HP'] = '{}-{},{}-{}'.format(component + 1, phase + 1, component + 1, 2 - phase)

	def _set_PS(self, values, component, phase):
		"""
		values -- tag dict to update
		component -- name of the component
		phase -- 0 or 1
		"""
		assert phase in [0, 1]
		values['PS'] = str(component + 1)
		values['GT'] = '0|1' if phase == 0 else '1|0'

	def write(self, chromosome, sample_superreads, sample_components):
		"""
		Add phasing information to all variants on a single chromosome.

		chromosome -- name of chromosome
		sample_superreads -- dictionary that maps each sample to a superread
		sample_components -- a dictionary that maps each sample to its connected components

			Each component in turn is a dict that maps each variant position to a
			component, where a component is identified by the position of its
			left-most variant

		Since coordinates within the superreads are used to identify variants,
		variants at duplicate positions (allowed by the VCF spec) are currently
		not supported.

		Returns a list of changed genotyes (i.e. a list of GenotypeChange objects)
		"""
		assert self._unprocessed_record is None or (self._unprocessed_record.CHROM == chromosome)

		genotype_changes = []

		if self._unprocessed_record is not None:
			records_iter = itertools.chain([self._unprocessed_record], self._reader_iter)
		else:
			records_iter = self._reader_iter
		allowed_alleles = frozenset({(0, 1), (1, 0), (0, 0), (1, 1)})
		sample_phases = dict()
		sample_genotypes = dict()
		for sample, superreads in sample_superreads.items():
			sample_phases[sample] = {
				v1.position: v1.allele for v1, v2 in zip(*superreads)
					if (v1.allele, v2.allele) in allowed_alleles
			}
			sample_genotypes[sample] = {
				v1.position: v1.allele + v2.allele for v1, v2 in zip(*superreads)
					if (v1.allele, v2.allele) in allowed_alleles
			}
		n = 0
		prev_pos = None
		INT_TO_UNPHASED_GT = { 0: '0/0', 1: '0/1', 2: '1/1', -1: '.' }
		for record in records_iter:
			n += 1
			pos, ref, alt = record.start, str(record.REF), str(record.ALT[0])

			if record.CHROM != chromosome:
				# save it for later
				self._unprocessed_record = record
				assert n != 1
				break

			if len(record.ALT) > 1:
				# we do not phase multiallelic sites currently
				is_phased = False
			elif pos == prev_pos:
				# duplicate position, skip it
				is_phased = False
			else:
				# Determine whether the variant is phased in any sample
				is_phased = True
				for sample in self._reader.samples:
					if sample in sample_superreads:
						components = sample_components[sample]
						phases = sample_phases[sample]
						if pos in components and pos in phases:
							break
				else:
					is_phased = False

			if self.tag == 'PS':
				# Remove any existing phasing from the GT field, no matter
				# whether we have phasing info for it or not
				samp_fmt = self._reader._format_cache[record.FORMAT]
				for sample, call in zip(self._reader.samples, record.samples):
					if sample not in sample_superreads:
						continue

					if '|' in call.data.GT:
						values = call.data._asdict()
						gt_fields = call.data.GT.split('|')
						values['GT'] = '/'.join(sorted(gt_fields))
						call.data = samp_fmt(**values)

			if is_phased:
				# Current PyVCF does not make it very easy to modify records/calls.
				# Add HP tag to FORMAT
				if self.tag not in record.FORMAT.split(':'):
					record.add_format(self.tag)
					if record.FORMAT not in self._reader._format_cache:
						self._reader._format_cache[record.FORMAT] = self._reader._parse_sample_format(record.FORMAT)
				samp_fmt = self._reader._format_cache[record.FORMAT]

				# Set phase tag for all samples
				for i, call in enumerate(record.samples):
					sample = self._reader.samples[i]
					if sample not in sample_superreads:
						# This sample has not been phased, leave unchanged
						continue
					components = sample_components[sample]
					phases = sample_phases[sample]
					genotypes = sample_genotypes[sample]

					if (hasattr(call.data, self.tag) and getattr(call.data, self.tag) is not None
							and not self._phase_tag_found_warned):
						logger.warning('Ignoring existing phasing information '
							'found in input VCF ({} tag exists).'.format(self.tag))
						self._phase_tag_found_warned = True

					values = call.data._asdict()
					is_het = call.is_het

					# is genotype to be changed?
					if (pos in genotypes) and (genotypes[pos] != call.gt_type):
						values['GT'] = INT_TO_UNPHASED_GT[genotypes[pos]]
						variant = VcfVariant(record.POS, record.REF, record.ALT[0])
						genotype_changes.append(GenotypeChange(sample, chromosome, variant, call.gt_type, genotypes[pos]))
						is_het = genotypes[pos] == 1

					if pos in components and pos in phases and is_het:
						self._set_phasing_tags(values, components[pos], phases[pos])
					else:
						# Unphased - set phase tag to '.'
						values[self.tag] = None
					call.data = samp_fmt(**values)
			self._writer.write_record(record)
			prev_pos = pos
		return genotype_changes


############ class to print computed genotypes,likelihoods (still needs to be improved...) ###########
# in input vcf, currently GT is still required..


class GenotypeVcfWriter:
	"""
	Read in a VCF file and write it back out with added genotyping information.

	Avoid reading in full chromosomes as that uses too much memory for
	multi-sample VCFs.
	"""
	def __init__(self, in_path, command_line, out_file=sys.stdout):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry.
		out_file -- Open file-like object to which VCF is written.
		"""
		self._reader = vcf.Reader(filename=in_path)
		
		if command_line is not None:
			if 'commandline' not in self._reader.metadata:
				self._reader.metadata['commandline'] = []
			command_line = command_line.replace('"', '')
			self._reader.metadata['commandline'] = [command_line]

		# add tag for genotype
		fmt = vcf.parser._Format(id='GT', num=1, type='String', desc='Genotype computed by whatshap genotyping algorithm.')
		self._reader.formats['GT'] = fmt

		# add tag for genotype quality
		fmt = vcf.parser._Format(id='GQ', num=1, type='Integer', desc='Phred scaled genotype quality computed by whatshap genotyping algorithm.')
		self._reader.formats['GQ'] = fmt

		# add tag for genotype likelihoods
		fmt = vcf.parser._Format(id='GL', num='G', type='Float', desc='log10-scaled likelihoods for genotypes: 0/0,0/1,1/1, computed by whatshap genotyping algorithm.')
		self._reader.formats['GL'] = fmt

		self._writer = vcf.Writer(out_file, template=self._reader)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)

	@property
	def samples(self):
		return self._reader.samples

	def write_genotypes(self, chromosome, variant_table, indels, leave_unchanged=False):
		"""
		Add genotyping information to all variants on a single chromosome.

		chromosome -- name of chromosome
		variant_table -- contains genotyping information for all accessible variant positions
		leave_unchanged -- if True, leaves records of current chromosome unchanged
		"""

		assert self._unprocessed_record is None or (self._unprocessed_record.CHROM == chromosome)

		if self._unprocessed_record is not None:
			records_iter = itertools.chain([self._unprocessed_record], self._reader_iter)
		else:
			records_iter = self._reader_iter

		# map positions to index
		genotyped_variants = dict()
		for i in range(len(variant_table)):
			genotyped_variants[variant_table.variants[i].position] = i

		n = 0
		prev_pos = None
		INT_TO_UNPHASED_GT = { 0: '0/0', 1: '0/1', 2: '1/1', -1: '.' }
		for record in records_iter:

			n += 1
			pos, ref, alt = record.start, str(record.REF), str(record.ALT[0])

			if record.CHROM != chromosome:
				# save it for later
				self._unprocessed_record = record
				assert n != 1
				break

			# if current chromosome was genotyped, write this new information to vcf
			if not leave_unchanged:
				no_format_field = False

				# add GT,GQ,GL fields in case they are not present yet
				if record.FORMAT == None:
					record.FORMAT = 'GT'
					no_format_field = True

				else:
					if 'GT' not in record.FORMAT.split(':'):
						record.add_format('GT')

				for field in ['GQ', 'GL']:
					if field not in record.FORMAT.split(':'):
						record.add_format(field)

				if record.FORMAT not in self._reader._format_cache:
					self._reader._format_cache[record.FORMAT] = self._reader._parse_sample_format(record.FORMAT)

				samp_fmt = self._reader._format_cache[record.FORMAT]

				# create call objects in case the format field was empty before
				if no_format_field:
					for sample in self._reader.samples:
						call = vcf.model._Call(record, sample, [])
						values = {'GT': '.', 'GQ': '.', 'GL': '.'}
						call.data = samp_fmt(**values)
						record.samples.append(call)


				for i, call in enumerate(record.samples):
					sample = self._reader.samples[i]
					values = call.data._asdict()

					geno = -1
					geno_l = [1/3.0] * 3
					geno_q = '.'

					# for genotyped variants, get computed likelihoods/genotypes (for all others, give uniform likelihoods)
					if pos in genotyped_variants:
						likelihoods = variant_table.genotype_likelihoods_of(sample)[genotyped_variants[pos]]
						# likelihoods can be 'None' if position was not accessible
						if not likelihoods==None:
							geno_l = likelihoods
							geno = variant_table.genotypes_of(sample)[genotyped_variants[pos]]

					# compute GQ
					if geno == 0:
						geno_q = geno_l[1] + geno_l[2]
					elif geno == 1:
						geno_q = geno_l[0] + geno_l[2]
					elif geno == 2:
						geno_q = geno_l[0] + geno_l[1]

					# store genotype
					values['GT'] = INT_TO_UNPHASED_GT[geno]
					# store quality as phred score
					if not geno == -1:
						# TODO default value ok?
						if geno_q > 0:
							values['GQ'] = min(round(-10.0 * math.log10(geno_q)), 10000)
						else:
							values['GQ'] = 10000
					else:
						values['GQ'] = '.'

					# TODO default value ok?
					# store likelihoods log10-scaled
					values['GL'] = [max(math.log10(j),-1000) if j>0 else -1000 for j in geno_l]

					record.QUAL = '.'

					# delete all other genotype information that might have been present before
					for tag in record.FORMAT.split(':'):
						if tag not in ['GT', 'GL', 'GQ']:
							values[tag] = '.'

					call.data = samp_fmt(**values)
			self._writer.write_record(record)
			prev_pos = pos
