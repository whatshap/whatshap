"""
Functions for reading VCFs.
"""
import sys
import logging
import itertools
import math
from abc import ABC, abstractmethod
from array import array
from collections import namedtuple
from pysam import VariantFile
from typing import List

from .core import Read, PhredGenotypeLikelihoods

logger = logging.getLogger(__name__)


class VcfError(Exception):
	pass


class VcfNotSortedError(VcfError):
	pass


class VcfVariant:
	"""A variant in a VCF file (not to be confused with core.Variant)"""

	__slots__ = ('position', 'reference_allele', 'alternative_allele')

	def __init__(self, position: int, reference_allele: str, alternative_allele: str):
		"""
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
		return (self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)

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
			p = [10**x for x in (self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)]
			s = sum(p)
			p = [x / s + regularizer for x in p]
			m = max(p)
			return PhredGenotypeLikelihoods(*(round(-10*math.log10(x/m)) for x in p))


class VariantTable:
	"""
	For a single chromosome, store variants and their genotypes.
	Each row of this table contains a variant, each column
	contains the genotypes of a single sample.

	chromosome -- chromosome name
	samples -- list of sample names
	"""
	def __init__(self, chromosome: str, samples: List[str]):
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

	def add_variant(self, variant: VcfVariant, genotypes, phases, genotype_likelihoods):
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

	def genotypes_of(self, sample: str):
		"""Retrieve genotypes by sample name"""
		return self.genotypes[self._sample_to_index[sample]]

	def set_genotypes_of(self, sample: str, genotypes):
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

	def phases_of(self, sample: str):
		"""Retrieve phases by sample name"""
		return self.phases[self._sample_to_index[sample]]

	def num_of_blocks_of(self, sample: str):
		""" Retrieve the number of blocks of the sample"""
		return len(set([i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None]))

	def id_of(self, sample: str):
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
		to_discard = [i for i, v in enumerate(self.variants) if v.position not in positions]
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
			if variant not in input_variant_set:
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
		self._vcf_reader = VariantFile(path)
		self._phases = phases
		self._genotype_likelihoods = genotype_likelihoods
		self.samples = list(self._vcf_reader.header.samples)  # intentionally public
		self.ignore_genotypes = ignore_genotypes
		logger.debug("Found %d sample(s) in the VCF file.", len(self.samples))

	def _fetch(self, chromosome: str):
		"""
		Return VariantTable object for a given chromosome.
		"""
		records = list(self._vcf_reader.fetch(chromosome))
		return self._process_single_chromosome(chromosome, records)

	def __iter__(self):
		"""
		Yield VariantTable objects for each chromosome.

		Multi-ALT sites are skipped.
		"""
		for chromosome, records in itertools.groupby(self._vcf_reader, lambda record: record.chrom):
			yield self._process_single_chromosome(chromosome, records)

	@staticmethod
	def _extract_HP_phase(call):
		hp = call.get('HP')
		if hp is None or hp == ('.', ):
			return None
		assert len(hp) == 2, hp
		fields = [[int(x) for x in s.split('-')] for s in hp]
		assert fields[0][0] == fields[1][0]
		block_id = fields[0][0]
		phase1, phase2 = fields[0][1]-1, fields[1][1]-1
		assert (phase1, phase2) == (0, 1) or (phase1, phase2) == (1, 0)
		return VariantCallPhase(block_id=block_id, phase=phase1, quality=call.get('PQ', None))

	@staticmethod
	def _extract_GT_PS_phase(call):
		is_het = call['GT'] == (0, 1) or call['GT'] == (1, 0)
		if not is_het:
			return None
		if not call.phased:
			return None
		block_id = call.get('PS', 0)
		assert call['GT'] in ((0, 1), (1, 0))
		phase = call['GT'][0]
		return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get('PQ', None))

	def _process_single_chromosome(self, chromosome, records):
		phase_detected = None
		n_snvs = 0
		n_other = 0
		n_multi = 0
		table = VariantTable(chromosome, self.samples)
		prev_position = None
		for record in records:
			if len(record.alts) > 1:
				# Multi-ALT sites are not supported, yet
				n_multi += 1
				continue

			pos, ref, alt = record.start, str(record.ref), str(record.alts[0])
			if len(ref) == len(alt) == 1:
				n_snvs += 1
			else:
				n_other += 1
				if not self._indels:
					continue

			if (prev_position is not None) and (prev_position > pos):
				raise VcfNotSortedError('VCF not ordered: {}:{} appears before {}:{}'.format(
					chromosome, prev_position + 1, chromosome, pos + 1))

			if prev_position == pos:
				logger.warning('Skipping duplicated position %s on chromosome %r', pos + 1, chromosome)
				continue
			prev_position = pos

			# Read phasing information (allow GT/PS or HP phase information, but not both),
			# if requested
			if self._phases:
				phases = []
				for sample_name, call in record.samples.items():
					phase = None
					for extract_phase, phase_name in [
						(self._extract_HP_phase, 'HP'),
						(self._extract_GT_PS_phase, 'GT_PS')
					]:
						p = extract_phase(call)
						if p is not None:
							if phase_detected is None:
								phase_detected = phase_name
							elif phase_detected != phase_name:
								raise MixedPhasingError(
									"Mixed phasing information in input VCF (e.g. mixing PS "
									"and HP fields)")
							phase = p
					phases.append(phase)
			else:
				phases = [None] * len(record.samples)

			# Read genotype likelihoods, if requested
			if self._genotype_likelihoods:
				genotype_likelihoods = []
				for call in record.samples.values():
					GL = call.get('GL', None)
					PL = call.get('PL', None)
					# Prefer GLs (floats) over PLs (ints) if both should be present
					if GL is not None:
						assert len(GL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods(*GL))
					elif PL is not None:
						assert len(PL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods(*(pl/-10 for pl in PL)))
					else:
						genotype_likelihoods.append(None)
			else:
				genotype_likelihoods = [None] * len(record.samples)

			if not self.ignore_genotypes:
				genotypes = array('b', (genotype_code(call['GT'], unknown=-1) for call in record.samples.values()))
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


class VcfHeader:
	def __init__(self, format_or_info, id, number, typ, description):
		self.format_or_info = format_or_info
		self.id = id
		self.number =  number
		self.typ = typ
		self.description = description

	def line(self):
		return '##{format_or_info}=<ID={id},Number={number},Type={typ},'\
			'Description="{description}">'.format(
				format_or_info=self.format_or_info,
				id=self.id,
				number=self.number,
				typ=self.typ,
				description=self.description
			)


PREDEFINED_FORMATS = {
	'GL': VcfHeader('FORMAT', 'GL', 'G', 'Float', 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy'),
	'GQ': VcfHeader('FORMAT', 'GQ', 1, 'Integer', 'Phred-scaled genotype quality'),
	'GT': VcfHeader('FORMAT', 'GT', 1, 'String', 'Genotype'),
	'HP': VcfHeader('FORMAT', 'HP', '.', 'String', 'Phasing haplotype identifier'),
	'PQ': VcfHeader('FORMAT', 'PQ', 1, 'Float', 'Phasing quality'),
	'PS': VcfHeader('FORMAT', 'PS', 1, 'Integer', 'Phase set identifier'),
}

PREDEFINED_INFOS = {
	'AC': VcfHeader('INFO', 'AC', 'A', 'Integer', 'Allele count in genotypes, for each ALT allele, in the same order as listed'),
	'AN': VcfHeader('INFO', 'AN', 'A', 'Integer', 'Total number of alleles in called genotypes'),
	'END': VcfHeader('INFO', 'END', 1, 'Integer', 'Stop position of the interval'),
	'SVLEN': VcfHeader('INFO', 'SVLEN', '.', 'Integer', 'Difference in length between REF and ALT alleles'),
	'SVTYPE': VcfHeader('INFO', 'SVTYPE', 1, 'String', 'Type of structural variant'),
}


def augment_header(header, contigs, formats, infos):
	"""
	Add contigs, formats and infos to a VariantHeader.

	formats and infos are given as a list of strings, where each item is the ID of the header
	line to add. The full header info (Number, Type, Description) is taken from the PREDEFINED_*
	constants above. Any other FORMATs or INFOs that are not predefined will raise a VcfError.

	The header is modified in place.
	"""
	for contig in contigs:
		header.contigs.add(contig)

	for fmt in formats:
		if fmt in header.formats:
			header.formats[fmt].remove_header()
		try:
			h = PREDEFINED_FORMATS[fmt]
		except KeyError:
			raise VcfError("FORMAT {!r} not defined in VCF header".format(fmt)) from None
		header.add_line(h.line())

	for info in infos:
		try:
			h = PREDEFINED_INFOS[info]
		except KeyError:
			raise VcfError("INFO {!r} not defined in VCF header".format(info)) from None
		header.add_line(h.line())


def missing_headers(path):
	"""
	Find contigs, FORMATs and INFOs that are used within the body of a VCF file, but are
	not listed in the header or that have an incorrect type.

	Return a tuple (contigs, formats, infos) where each of the items are lists of
	strings.

	The reason this function exists is that pysam.VariantFile crashes when we
	try to write a VCF record to it that uses contigs, INFOs or FORMATs that
	are missing from the header. See also
	<https://github.com/pysam-developers/pysam/issues/771>
	"""
	with VariantFile(path) as variant_file:
		header = variant_file.header.copy()
		# Check for FORMATs that do not have the expected type
		incorrect_formats = []
		for fmt, v in variant_file.header.formats.items():
			if fmt not in PREDEFINED_FORMATS:
				continue
			h = PREDEFINED_FORMATS[fmt]
			if v.number != h.number or v.type != h.typ:
				incorrect_formats.append(fmt)

		# Iterate through entire file and check which contigs, formats and
		# info fields are used
		contigs = []  # contigs encountered, in the proper order
		seen_contigs = set()
		formats = []  # FORMATs encountered, in the proper order
		seen_formats = set()
		seen_infos = set()  # INFOs encountered

		for record in variant_file:
			seen_infos.update(record.info)
			if record.alts is not None:
				for alt in record.alts:
					# If there are "vague" ALT alleles such as <INS>, <DEL> etc, then
					# the header needs to contain a LEN info entry even if LEN
					# is never used
					if alt.startswith('<'):
						seen_infos.add('END')

			# For the contigs, we maintain a list *and* a set because we want to
			# keep track of the order of the contigs.
			if record.contig not in seen_contigs:
				contigs.append(record.contig)
			seen_contigs.add(record.contig)

			for fmt in record.format:
				if fmt not in seen_formats:
					formats.append(fmt)
				seen_formats.add(fmt)

	# Determine which contigs are missing from the header
	header_contigs = set(header.contigs)
	missing_contigs = []
	for contig in contigs:
		if contig not in header_contigs:
			missing_contigs.append(contig)

	# Determine which FORMATs are missing from the header
	header_formats = set(header.formats)
	missing_formats = []
	for fmt in formats:
		if fmt in header_formats:
			continue
		missing_formats.append(fmt)

	# Determine which INFOs are missing from the header
	missing_infos = list(set(seen_infos) - set(header.info))

	return (missing_contigs, incorrect_formats + missing_formats, missing_infos)


GenotypeChange = namedtuple('GenotypeChange', ['sample', 'chromosome', 'variant', 'old_gt', 'new_gt'])


class VcfAugmenter(ABC):

	def __init__(self, in_path, command_line, out_file=sys.stdout):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry
			(use None to not add this to the VCF header)
		out_file -- Open file-like object to which VCF is written.
		tag -- which type of tag to write, either 'PS' or 'HP'. 'PS' is standardized;
		    'HP' is compatible with GATK’s ReadBackedPhasing.
		"""
		# TODO This is slow because it reads in the entire VCF one extra time
		contigs, formats, infos = missing_headers(in_path)
		# We repair the header (adding missing contigs, formats, infos) of the *input* VCF because
		# we will modify the records that we read, and these are associated with the input file.
		self._reader = VariantFile(in_path)
		augment_header(self._reader.header, contigs, formats, infos)
		if command_line is not None:
			command_line = '"' + command_line.replace('"', '') + '"'
			self._reader.header.add_meta('commandline', command_line)
		self.setup_header(self._reader.header)
		self._writer = VariantFile(out_file, mode='w', header=self._reader.header)
		self._unprocessed_record = None
		self._reader_iter = iter(self._reader)

	@abstractmethod
	def setup_header(self, header):
		pass

	def close(self):
		self._writer.close()

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()

	@property
	def samples(self):
		return list(self._reader.header.samples)

	def _iterrecords(self, chromosome):
		"""Yield all records for the target chromosome"""
		n = 0
		if self._unprocessed_record is not None:
			assert self._unprocessed_record.chrom == chromosome
			yield self._unprocessed_record
			n += 1
		for record in self._reader_iter:
			n += 1
			if record.chrom != chromosome:
				# save it for later
				self._unprocessed_record = record
				assert n != 1
				return
			yield record


class PhasedVcfWriter(VcfAugmenter):
	"""
	Read in a VCF file and write it back out with added phasing information.

	Avoid reading in full chromosomes as that uses too much memory for
	multi-sample VCFs.
	"""
	def __init__(self, in_path, command_line, out_file=sys.stdout, tag='PS'):
		"""
		in_path -- Path to input VCF, used as template.
		command_line -- A string that will be added as a VCF header entry
			(use None to not add this to the VCF header)
		out_file -- Open file-like object to which VCF is written.
		tag -- which type of tag to write, either 'PS' or 'HP'. 'PS' is standardized;
			'HP' is compatible with GATK’s ReadBackedPhasing.
		"""
		if tag not in ('HP', 'PS'):
			raise ValueError('Tag must be either "HP" or "PS"')
		self.tag = tag
		super().__init__(in_path, command_line, out_file)
		self._phase_tag_found_warned = False
		self._set_phasing_tags = self._set_HP if tag == 'HP' else self._set_PS

	def setup_header(self, header):
		"""Called by baseclass constructor"""

		# FreeBayes adds phasing=none to its VCF output - remove that.
		for hr in header.records:
			if hr.key == 'phasing':
				hr.remove()
				break

		header.add_line(PREDEFINED_FORMATS[self.tag].line())

	def _set_HP(self, call, component, phase):
		"""
		values -- tag dict to update
		component -- name of the component
		phase -- 0 or 1
		"""
		assert phase in [0, 1]
		call['HP'] = '{}-{},{}-{}'.format(component + 1, phase + 1, component + 1, 2 - phase)

	def _set_PS(self, call, component, phase):
		"""
		values -- tag dict to update
		component -- name of the component
		phase -- 0 or 1
		"""
		assert phase in [0, 1]
		call['PS'] = component + 1
		call['GT'] = (0, 1) if phase == 0 else (1, 0)
		call.phased = True

	def write(self, chromosome: str, sample_superreads, sample_components):
		"""
		Add phasing information to all variants on a single chromosome.

		chromosome -- name of chromosome
		sample_superreads -- dictionary that maps a sample name to a superread
		sample_components -- a dictionary that maps a sample to its connected components

			Each component in turn is a dict that maps each variant position to a
			component, where a component is identified by the position of its
			left-most variant

		Since coordinates within the superreads are used to identify variants,
		variants at duplicate positions (allowed by the VCF spec) are currently
		not supported.

		Returns a list of changed genotyes (i.e. a list of GenotypeChange objects)
		"""
		genotype_changes = []
		allowed_alleles = frozenset({(0, 1), (1, 0), (0, 0), (1, 1)})
		sample_phases = dict()
		sample_genotypes = dict()
		for sample, superreads in sample_superreads.items():
			alleles = {
				(v1.position, v1.allele, v2.allele) for v1, v2 in zip(*superreads)
				if (v1.allele, v2.allele) in allowed_alleles
			}
			sample_phases[sample] = {position: allele1 for position, allele1, allele2 in alleles}
			# Map position to encoded genotype
			sample_genotypes[sample] = {
				position: allele1 + allele2 for position, allele1, allele2 in alleles}
		prev_pos = None
		INT_TO_UNPHASED_GT = {0: (0, 0), 1: (0, 1), 2: (1, 1), -1: None}
		for record in self._iterrecords(chromosome):
			pos, ref, alt = record.start, record.ref, record.alts[0]

			if len(record.alts) > 1:
				# we do not phase multiallelic sites currently
				is_phased = False
			elif pos == prev_pos:
				# duplicate position, skip it
				is_phased = False
			else:
				# Determine whether the variant is phased in any sample
				is_phased = True
				for sample in self.samples:
					if sample in sample_superreads:
						components = sample_components[sample]
						phases = sample_phases[sample]
						if pos in components and pos in phases:
							break
				else:
					is_phased = False

			if self.tag == 'PS':
				# Remove any existing phasing from the GT field
				for sample, call in record.samples.items():
					if sample in sample_superreads:
						call.phased = False
						if call['GT'] is not None and call['GT'][0] is not None and call['GT'][1] is not None:
							call['GT'] = sorted(call['GT'])

			if is_phased:
				# Set phase tag for all target samples
				for sample, call in record.samples.items():
					if sample not in sample_superreads:
						# This sample has not been phased, leave unchanged
						continue
					components = sample_components[sample]
					phases = sample_phases[sample]
					genotypes = sample_genotypes[sample]

					if (self.tag in call and call[self.tag] is not None
							and not self._phase_tag_found_warned):
						logger.warning('Ignoring existing phasing information '
							'found in input VCF ({} tag exists).'.format(self.tag))
						self._phase_tag_found_warned = True

					gt_type = genotype_code(call['GT'])
					is_het = gt_type == 1

					# is genotype to be changed?
					if pos in genotypes and genotypes[pos] != gt_type:
						call['GT'] = INT_TO_UNPHASED_GT[genotypes[pos]]
						variant = VcfVariant(record.start, record.ref, record.alts[0])
						genotype_changes.append(
							GenotypeChange(sample, chromosome, variant, gt_type, genotypes[pos]))
						is_het = genotypes[pos] == 1

					if pos in components and pos in phases and is_het:
						self._set_phasing_tags(call, components[pos], phases[pos])
					else:
						# Unphased
						call[self.tag] = None
			self._writer.write(record)
			prev_pos = pos
		return genotype_changes


def genotype_code(gt, unknown=None):
	"""Return genotype encoded as PyVCF-compatible number"""
	if gt is None:
		return unknown
	if gt[0] is None or gt[1] is None:
		return unknown
	if gt == (0, 0):
		return 0  # homozygous ref
	if gt == (1, 1):
		return 2  # homozygous alt
	assert gt == (0, 1) or gt == (1, 0)
	return 1  # heterozygous


############ class to print computed genotypes,likelihoods (still needs to be improved...) ###########
# in input vcf, currently GT is still required..


class GenotypeVcfWriter(VcfAugmenter):
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
		super().__init__(in_path, command_line, out_file)

	def setup_header(self, header):
		"""Called by baseclass constructor"""
		header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype computed by WhatsHap genotyping algorithm">')
		header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled genotype quality computed by WhatsHap genotyping algorithm">')
		header.add_line('##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled likelihoods for genotypes: 0/0, 0/1, 1/1, computed by WhatsHap genotyping algorithm">')

	def write_genotypes(self, chromosome, variant_table, indels, leave_unchanged=False):
		"""
		Add genotyping information to all variants on a single chromosome.

		chromosome -- name of chromosome
		variant_table -- contains genotyping information for all accessible variant positions
		leave_unchanged -- if True, leaves records of current chromosome unchanged
		"""

		# map positions to index
		genotyped_variants = dict()
		for i in range(len(variant_table)):
			genotyped_variants[variant_table.variants[i].position] = i

		INT_TO_UNPHASED_GT = {0: (0, 0), 1: (0, 1), 2: (1, 1), -1: None}
		GT_GL_GQ = frozenset(['GT', 'GL', 'GQ'])
		for record in self._iterrecords(chromosome):
			pos = record.start

			# if current chromosome was genotyped, write this new information to VCF
			if not leave_unchanged:
				for sample, call in record.samples.items():
					geno = -1
					n_alleles = len(record.alts) + 1
					n_genotypes = (n_alleles*n_alleles + n_alleles)/2
					geno_l = [1/n_genotypes] * int(n_genotypes)
					geno_q = None

					# for genotyped variants, get computed likelihoods/genotypes (for all others, give uniform likelihoods)
					if pos in genotyped_variants:
						likelihoods = variant_table.genotype_likelihoods_of(sample)[genotyped_variants[pos]]
						# likelihoods can be 'None' if position was not accessible
						if likelihoods is not None:
							assert n_alleles == 2
							geno_l = likelihoods
							geno = variant_table.genotypes_of(sample)[genotyped_variants[pos]]

					# Compute GQ
					if geno == 0:
						geno_q = geno_l[1] + geno_l[2]
					elif geno == 1:
						geno_q = geno_l[0] + geno_l[2]
					elif geno == 2:
						geno_q = geno_l[0] + geno_l[1]

					# TODO default value ok?
					# store likelihoods log10-scaled

					# Temporarily overwrite the GT field with a (fake) genotype that indicates a
					# diploid sample. Otherwise, if the GT field happens to be empty, pysam
					# complains that we are setting an incorrect number of GL values.
					call['GT'] = (0, 1)

					call['GL'] = [max(math.log10(j), -1000) if j > 0 else -1000 for j in geno_l]
					call['GT'] = INT_TO_UNPHASED_GT[geno]

					# store quality as phred score
					if not geno == -1:
						# TODO default value ok?
						assert geno_q is not None
						if geno_q > 0:
							call['GQ'] = min(round(-10.0 * math.log10(geno_q)), 10000)
						else:
							call['GQ'] = 10000
					else:
						call['GQ'] = None

					record.qual = None

					# delete all other genotype information that might have been present before
					for tag in set(call.keys()) - GT_GL_GQ:
						del call[tag]

			self._writer.write(record)
