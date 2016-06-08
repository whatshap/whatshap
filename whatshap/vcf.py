"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
from array import array
from collections import defaultdict, namedtuple
import vcf
from .core import Read

logger = logging.getLogger(__name__)


class VcfVariant:
	__slots__ = ('position', 'reference_allele', 'alternative_allele')

	"""A variant in a VCF file (not to be confused with core.Variant)"""
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
		return "VcfVariant(pos={}, ref={!r}, alt={!r})".format(self.position,
			self.reference_allele, self.alternative_allele)

	def __hash__(self):
		return hash((self.position, self.reference_allele, self.alternative_allele))

	def __eq__(self, other):
		return (self.position == other.position) and (self.reference_allele == other.reference_allele) and (self.alternative_allele == other.alternative_allele)

	def __lt__(self, other):
		return (self.position, self.reference_allele, self.alternative_allele) < (other.position, other.reference_allele, other.alternative_allele)


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
		self.phases = [ [] for _ in samples ]
		self.variants = []
		self._sample_to_index = { sample: index for index, sample in enumerate(samples) }

	def __len__(self):
		return len(self.variants)

	#def add_sample(self, name, genotypes):
		#"Add a column to the table"
		#if len(genotypes) != len(self.variants):
			#raise ValueError('Expecting as many genotypes as there are variants')
		#self._name_to_index[name] = len(self.samples)
		#self.samples.append(name)
		#self.genotypes.append(genotypes)

	def add_variant(self, variant, genotypes, phases):
		"""
		Add a row to the table

		variant -- a VcfVariant
		genotypes -- iterable of ints that encode the genotypes of the samples:
			-1 represents an unknown genotype
			0 represents 0/0 (homozygous reference)
			1 represents 0/1 or 1/0 (heterozygous)
			2 represents 1/1 (homozygous alternative)
		phases -- iterable yielding a pair (block_id, phase)
			for each sample, where and block_id is a numeric id of the phased block and
			phase is either 0 or 1 (indicating whether the REF allele is on haplotype 0 or 1).
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

	def genotypes_of(self, sample):
		"""Retrieve genotypes by sample name"""
		return self.genotypes[self._sample_to_index[sample]]

	def phases_of(self, sample):
		"""Retrieve phases by sample name"""
		return self.phases[self._sample_to_index[sample]]

	def num_of_blocks_of(self, sample):
		""" Retrieve the number of blocks of the sample"""
		return len(set([i.block_id for i in self.phases[self._sample_to_index[sample]] if not i is None]))

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

		for gt in self.genotypes:
			assert len(self.variants) == len(gt)
		for ph in self.phases:
			assert len(self.variants) == len(ph)
		assert len(self.samples) == len(self.genotypes) == len(self.phases)

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
				r = Read('{}_block_{}'.format(sample,phase.block_id), mapq, source_id, numeric_sample_id)
				r.add_variant(variant.position, phase.phase, quality)
				read_map[phase.block_id] = r
		for key, read in read_map.items():
			read.sort()
			if len(read) > 1:
				yield read



class MixedPhasingError(Exception):
	pass


VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'phase', 'quality'])


class VcfReader:
	"""
	Read a VCF file chromosome by chromosome.
	"""
	def __init__(self, path, indels=False, normalize=False):
		"""
		path -- Path to VCF file
		indels -- Whether to include also insertions and deletions in the list of
			variants.
		normalize -- Whether to normalize variants
		"""
		# TODO Always include deletions since they can 'overlap' other variants
		self._indels = indels
		self._vcf_reader = vcf.Reader(filename=path)
		self._normalize = normalize
		self.samples = self._vcf_reader.samples  # intentionally public
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

	@staticmethod
	def normalize(pos, ref, alt):
		"""
		Normalize variants that share a common prefix and/or suffix.
		For example, GCTGTT -> GCTAAATT is changed to G -> AAA.

		Return a (pos, ref, alt) tuple.
		"""
		while len(ref) >= 1 and len(alt) >= 1 and ref[-1] == alt[-1]:
			ref, alt = ref[:-1], alt[:-1]

		while len(ref) >= 1 and len(alt) >= 1 and ref[0] == alt[0]:
			ref, alt = ref[1:], alt[1:]
			pos += 1

		return pos, ref, alt

	def __iter__(self):
		"""
		Yield VariantTable objects for each chromosome.

		Multi-ALT sites are skipped.
		"""
		for chromosome, records in self._group_by_chromosome():
			yield self._process_single_chromosome(chromosome, records)


	@staticmethod
	def _read_PQ(call):
		if (not hasattr(call.data,'PQ')) or (call.data.PQ is None):
			return None
		return call.data.PQ

	@staticmethod
	def _extract_HP_phase(call):
		if (not hasattr(call.data,'HP')) or (call.data.HP is None):
			return None
		HP = call.data.HP
		assert len(HP) == 2
		fields = [[int(x) for x in s.split('-')] for s in HP]
		assert fields[0][0] == fields[1][0]
		block_id = fields[0][0]
		phase1, phase2 = fields[0][1]-1, fields[1][1]-1
		assert ((phase1, phase2) == (0, 1)) or ((phase1, phase2) == (1, 0))
		return VariantCallPhase(block_id=block_id, phase=phase1, quality=VcfReader._read_PQ(call))

	@staticmethod
	def _extract_GT_PS_phase(call):
		if not call.is_het:
			return None
		if not call.phased:
			return None
		block_id = 0
		if (hasattr(call.data,'PS')) and (not call.data.PS is None):
			block_id = call.data.PS
		assert call.data.GT in ['0|1','1|0']
		phase = int(call.data.GT[0])
		return VariantCallPhase(block_id=block_id, phase=phase, quality=VcfReader._read_PQ(call))

	def _process_single_chromosome(self, chromosome, records):
		phase_detected = None
		n_snps = 0
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
				n_snps += 1
			else:
				n_other += 1
				if not self._indels:
					continue

			if prev_position == pos:
				logger.warning('Skipping duplicated position %s on chromosome %r', pos+1, chromosome)
				continue
			prev_position = pos
			if self._normalize:
				pos, ref, alt = self.normalize(pos, ref, alt)

			# PyVCF pecularity: gt_alleles is a list of the alleles in the
			# GT field, but as strings.
			# For example, when GT is 0/1, gt_alleles is ['0', '1'].
			# And when GT is 2|1, gt_alleles is ['2', '1'].
			"""
			logger.debug("Call %s:%d %s→%s",
				record.CHROM, record.start + 1, record.REF, record.ALT)
			"""

			# Read phasing information (allow GT/PS or HP phase information, but not both)
			phases = []
			for call in record.samples:
				phase = None
				for extract_phase, phase_name in [(self._extract_HP_phase, 'HP'), (self._extract_GT_PS_phase, 'GT_PS')]:
					p = extract_phase(call)
					if not p is None:
						if phase_detected is None:
							phase_detected = phase_name
						elif phase_detected != phase_name:
							raise MixedPhasingError('Mixed phasing information in input VCF (e.g. mixing PS and HP fields)')
						phase = p
				phases.append(phase)

			GT_TO_INT = { 0: 0, 1: 1, 2: 2, None: -1 }
			genotypes = array('b', (GT_TO_INT[call.gt_type] for call in record.samples))
			variant = VcfVariant(position=pos, reference_allele=ref, alternative_allele=alt)
			table.add_variant(variant, genotypes, phases)

		logger.debug("Parsed %s SNPs and %s non-SNPs. Also skipped %s multi-ALTs.", n_snps,
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


class PhasedVcfWriter:
	"""
	Read in a VCF file and write it back out with added phasing information.
	Phasing is written into HP and PQ tags, compatible with GATK’s
	ReadBackedPhasing.

	Avoid reading in full chromosomes as that uses too much memory for
	multi-sample VCFs.
	"""
	def __init__(self, in_path, command_line, normalized, out_file=sys.stdout):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry.
		out_file -- File-like object to which VCF is written.
		normalized -- whether the phased variants have been normalized
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
		self._normalized = normalized

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

		Coordinates within the superreads are used to identify variants. Since
		variant normalization changes coordinates, make sure that the PhasedVcfWriter
		has been initialized with the correct 'normalized' setting!

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
		prev_pos = None
		for record in records_iter:
			n += 1
			pos, ref, alt = record.start, str(record.REF), str(record.ALT[0])
			if self._normalized:
				pos, ref, alt = VcfReader.normalize(pos, ref, alt)

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
					if pos in components and pos in phases and call.is_het:
						# Set or overwrite HP tag
						values['HP'] = self._format_phasing_info(
							components[pos], phases[pos])
					else:
						# Unphased - set HP to '.'
						values['HP'] = None
					call.data = samp_fmt(**values)
			self._writer.write_record(record)
			prev_pos = pos
