"""
Print phasing statistics of a single VCF file
"""
import logging
from collections import defaultdict, namedtuple
from contextlib import ExitStack
from .math import median
from .vcf import VcfReader

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('--gtf', default=None, help='Write phased blocks to GTF file.')
	add('--sample', metavar='SAMPLE', default=None, help='Name of the sample '
			'to process. If not given, use first sample found in VCF.')
	add('--chr-lengths', default=None, help='File with chromosome lengths '
			'(one line per chromosome, tab separated "<chr> <length>") '
			'needed to compute N50 values.')
	add('--tsv', metavar='TSV', default=None, help='Filename to write '
		'statistics to (tab-separated).')
	add('--only-snvs', default=False, action="store_true", help='Only process SNVs '
		'and ignore all other variants.')
	add('--block-list', action="store", default=None, help='Filename to write '
		'list of all blocks to (one block per line).')
	add('--chromosome', dest='chromosome', metavar='CHROMOSOME', default=None,
		help='Name of chromosome to process. If not given, all chromosomes in the '
		'input VCF are considered.')
	add('vcf', metavar='VCF', help='Phased VCF file')


class PhasedBlock:
	def __init__(self, chromosome=None):
		self.phases = {}
		self.leftmost_variant = None
		self.rightmost_variant = None
		self.chromosome = chromosome

	def add(self, variant, phase):
		if len(self.phases) == 0:
			self.leftmost_variant = variant
			self.rightmost_variant = variant
		else:
			if variant < self.leftmost_variant:
				self.leftmost_variant = variant
			if self.rightmost_variant < variant:
				self.rightmost_variant = variant
		self.phases[variant] = phase

	def span(self):
		"""Returns the length of the covered genomic region in bp."""
		return self.rightmost_variant.position - self.leftmost_variant.position

	def variants(self):
		return list(sorted(self.phases.keys()))

	def count_snvs(self):
		return sum(int(variant.is_snv()) for variant in self.phases)

	def __repr__(self):
		return "PhasedBlock({})".format(str(self.phases))

	def __len__(self):
		return len(self.phases)

	def __lt__(self, other):
		return (self.leftmost_variant, self.rightmost_variant) < (other.leftmost_variant, other.rightmost_variant)


class GtfWriter:
	def __init__(self, file):
		self._file = file

	def write(self, chromosome, start, stop, name):
		"""
		Write a feature to the GTF. start is 0-based.
		"""
		assert start < stop
		print(chromosome, 'Phasing', 'exon', start + 1, stop, '.', '+', '.',
			'gene_id "{}"; transcript_id "{}.1";'.format(name, name),
			sep='\t', file=self._file)


detailed_stats_fields = [
	'variants',
	'phased',
	'unphased',
	'singletons',
	'blocks',
	'variant_per_block_median',
	'variant_per_block_avg',
	'variant_per_block_min',
	'variant_per_block_max',
	'variant_per_block_sum',
	'bp_per_block_median',
	'bp_per_block_avg',
	'bp_per_block_min',
	'bp_per_block_max',
	'bp_per_block_sum',
	'heterozygous_variants',
	'heterozygous_snvs',
	'phased_snvs',
	'block_n50'
]
DetailedStats = namedtuple('DetailedStats', detailed_stats_fields)


def compute_n50(blocks, chr_lengths):
	chromosomes = set(b.chromosome for b in blocks)
	target_length = 0
	for chromosome in sorted(chromosomes):
		try:
			target_length += chr_lengths[chromosome]
		except KeyError:
			logger.warning('Not able to compute N50 because chromosome length of %s not available', chromosome)
			return float('nan')

	# Cut interleaved blocks to avoid inflating N50 in this case
	pos_sorted = sorted(blocks, key=lambda b: (b.chromosome, b.leftmost_variant.position) )
	block_lengths = []
	for i, block in enumerate(pos_sorted):
		if len(block) < 2:
			continue
		start, end = block.leftmost_variant.position, block.rightmost_variant.position
		if i+1 < len(pos_sorted):
			next_block = pos_sorted[i+1]
			if end > next_block.leftmost_variant.position:
				#logger.warning('Blocks are interleaved, cutting first block: end=%s --> %s',  end, next_block.leftmost_variant.position)
				end = next_block.leftmost_variant.position
		block_lengths.append(end-start)
	block_lengths.sort(reverse=True)
	s = 0
	for l in block_lengths:
		s += l
		if s >= 0.5*target_length:
			return l

	return 0


class PhasingStats:
	def __init__(self):
		self.blocks = []
		self.unphased = 0
		self.variants = 0
		self.heterozygous_variants = 0
		self.heterozygous_snvs = 0
		self.phased_snvs = 0

	def __iadd__(self, other):
		self.blocks.extend(other.blocks)
		self.unphased += other.unphased
		self.variants += other.variants
		self.heterozygous_variants += other.heterozygous_variants
		self.heterozygous_snvs += other.heterozygous_snvs
		self.phased_snvs += other.phased_snvs
		return self

	def add_blocks(self, blocks):
		self.blocks.extend(blocks)

	def add_unphased(self, unphased: int=1):
		self.unphased += unphased

	def add_variants(self, variants: int):
		self.variants += variants

	def add_heterozygous_variants(self, variants: int):
		self.heterozygous_variants += variants

	def add_heterozygous_snvs(self, snvs: int):
		self.heterozygous_snvs += snvs

	def get(self, chr_lengths=None):
		"""Return DetailedStats
		"""
		block_sizes = sorted(len(block) for block in self.blocks)
		n_singletons = sum(1 for size in block_sizes if size == 1)
		block_sizes = [size for size in block_sizes if size > 1]
		block_lengths = sorted(block.span() for block in self.blocks if len(block) > 1)
		phased_snvs = sum(block.count_snvs() for block in self.blocks if len(block) > 1)
		if block_sizes:
			return DetailedStats(
				variants = self.variants,
				phased = sum(block_sizes),
				unphased = self.unphased,
				singletons = n_singletons,
				blocks = len(block_sizes),
				variant_per_block_median = median(block_sizes),
				variant_per_block_avg = sum(block_sizes) / len(block_sizes),
				variant_per_block_min = block_sizes[0],
				variant_per_block_max = block_sizes[-1],
				variant_per_block_sum = sum(block_sizes),
				bp_per_block_median = median(block_lengths),
				bp_per_block_avg = sum(block_lengths) / len(block_lengths),
				bp_per_block_min = block_lengths[0],
				bp_per_block_max = block_lengths[-1],
				bp_per_block_sum = sum(block_lengths),
				heterozygous_variants = self.heterozygous_variants,
				heterozygous_snvs = self.heterozygous_snvs,
				phased_snvs = phased_snvs,
				block_n50 = compute_n50(self.blocks, chr_lengths) if chr_lengths is not None else float('nan')
			)
		else:
			return DetailedStats(
				variants = self.variants,
				phased = 0,
				unphased = self.unphased,
				singletons = n_singletons,
				blocks = 0,
				variant_per_block_median = float('nan'),
				variant_per_block_avg = float('nan'),
				variant_per_block_min = 0,
				variant_per_block_max = 0,
				variant_per_block_sum = 0,
				bp_per_block_median = float('nan'),
				bp_per_block_avg = float('nan'),
				bp_per_block_min = 0,
				bp_per_block_max = 0,
				bp_per_block_sum = 0,
				heterozygous_variants = self.heterozygous_variants,
				heterozygous_snvs = self.heterozygous_snvs,
				phased_snvs = 0,
				block_n50 = float('nan')
			)

	def print(self, chr_lengths=None):
		stats = self.get(chr_lengths)
		WIDTH = 21
		print('Variants in VCF:'.rjust(WIDTH), '{:8d}'.format(stats.variants))
		print('Heterozygous:'.rjust(WIDTH), '{:8d} ({:8d} SNVs)'.format(stats.heterozygous_variants, stats.heterozygous_snvs))
		print('Phased:'.rjust(WIDTH), '{:8d} ({:8d} SNVs)'.format(stats.phased, stats.phased_snvs))
		print('Unphased:'.rjust(WIDTH), '{:8d}'.format(stats.unphased), '(not considered below)')
		print('Singletons:'.rjust(WIDTH), '{:8d}'.format(stats.singletons), '(not considered below)')
		print('Blocks:'.rjust(WIDTH), '{:8d}'.format(stats.blocks))
		print()
		print('Block sizes (no. of variants)')
		print('Median block size:'.rjust(WIDTH), '{:11.2f} variants'.format(stats.variant_per_block_median))
		print('Average block size:'.rjust(WIDTH), '{:11.2f} variants'.format(stats.variant_per_block_avg))
		print('Largest block:'.rjust(WIDTH), '{:8d}    variants'.format(stats.variant_per_block_max))
		print('Smallest block:'.rjust(WIDTH), '{:8d}    variants'.format(stats.variant_per_block_min))
		print()
		print('Block lengths (basepairs)')
		print('Sum of lengths:'.rjust(WIDTH), '{:8d}    bp'.format(stats.bp_per_block_sum))
		print('Median block length:'.rjust(WIDTH), '{:11.2f} bp'.format(stats.bp_per_block_median))
		print('Average block length:'.rjust(WIDTH), '{:11.2f} bp'.format(stats.bp_per_block_avg))
		print('Longest block:'.rjust(WIDTH), '{:8d}    bp'.format(stats.bp_per_block_max))
		print('Shortest block:'.rjust(WIDTH), '{:8d}    bp'.format(stats.bp_per_block_min))
		print('Block N50:'.rjust(WIDTH), '{:8.0f}    bp'.format(stats.block_n50))
		assert stats.phased + stats.unphased + stats.singletons == stats.heterozygous_variants


def parse_chr_lengths(filename):
	chr_lengths = {}
	with open(filename) as f:
		for line in f:
			fields = line.split('\t')
			assert len(fields) == 2
			chr_lengths[fields[0]] = int(fields[1])
	return chr_lengths


def main(args):
	gtfwriter = tsv_file = block_list_file = None
	with ExitStack() as stack:
		if args.gtf:
			gtf_file = stack.enter_context(open(args.gtf, 'wt'))
			gtfwriter = GtfWriter(gtf_file)
		if args.tsv:
			tsv_file = stack.enter_context(open(args.tsv, 'w'))
		if args.block_list:
			block_list_file = stack.enter_context(open(args.block_list, 'w'))

		if args.chr_lengths:
			chr_lengths = parse_chr_lengths(args.chr_lengths)
			logger.info('Read length of %d chromosomes from %s', len(chr_lengths), args.chr_lengths)
		else:
			chr_lengths = None

		vcf_reader = VcfReader(args.vcf, phases=True, indels=not args.only_snvs)
		if len(vcf_reader.samples) == 0:
			logger.error('Input VCF does not contain any sample')
			return 1
		else:
			logger.info('Found {} sample(s) in input VCF'.format(len(vcf_reader.samples)))
		if args.sample:
			if args.sample in vcf_reader.samples:
				sample = args.sample
			else:
				logger.error('Requested sample ({}) not found'.format(args.sample))
				return 1
		else:
			sample = vcf_reader.samples[0]
			logger.info('Reporting results for sample {}'.format(sample))

		if tsv_file:
			print('#sample', 'chromosome', 'file_name', *detailed_stats_fields, sep='\t', file=tsv_file)

		if block_list_file:
			print('#sample', 'chromosome', 'phase_set', 'from', 'to', 'variants', sep='\t', file=block_list_file)

		print('Phasing statistics for sample {} from file {}'.format(sample, args.vcf))
		total_stats = PhasingStats()
		chromosome_count = 0
		for variant_table in vcf_reader:
			if args.chromosome:
				if variant_table.chromosome != args.chromosome:
					continue
			chromosome_count += 1
			chromosome = variant_table.chromosome
			stats = PhasingStats()
			print('---------------- Chromosome {} ----------------'.format(chromosome))
			genotypes = variant_table.genotypes_of(sample)
			phases = variant_table.phases_of(sample)
			assert len(genotypes) == len(phases) == len(variant_table.variants)
			blocks = defaultdict(PhasedBlock)
			prev_block_id = None
			prev_block_fragment_start = None
			prev_block_fragment_end = None
			for variant, genotype, phase in zip(variant_table.variants, genotypes, phases):
				stats.add_variants(1)
				if genotype != 1:
					continue
				stats.add_heterozygous_variants(1)
				if variant.is_snv():
					stats.add_heterozygous_snvs(1)
				if phase is None:
					stats.add_unphased()
				else:
					# a phased variant
					blocks[phase.block_id].add(variant, phase)
					if gtfwriter:
						if prev_block_id is None:
							prev_block_fragment_start = variant.position
							prev_block_fragment_end = variant.position + 1
							prev_block_id = phase.block_id
						else:
							if prev_block_id != phase.block_id:
								gtfwriter.write(chromosome, prev_block_fragment_start, prev_block_fragment_end, prev_block_id)
								prev_block_fragment_start = variant.position
								prev_block_id = phase.block_id
							prev_block_fragment_end = variant.position + 1

			# Add chromosome information to each block. This is needed to
			# sort blocks later when we compute N50s
			for block_id, block in blocks.items():
				block.chromosome = chromosome

			if gtfwriter and prev_block_id is not None:
				gtfwriter.write(chromosome, prev_block_fragment_start, prev_block_fragment_end, prev_block_id)

			if block_list_file:
				block_ids = sorted(blocks.keys())
				for block_id in block_ids:
					print(sample, chromosome, block_id, blocks[block_id].leftmost_variant.position + 1, blocks[block_id].rightmost_variant.position + 1, len(blocks[block_id]), sep='\t', file=block_list_file)

			stats.add_blocks(blocks.values())
			stats.print(chr_lengths)
			if tsv_file:
				print(sample, chromosome, args.vcf, sep='\t', end='\t', file=tsv_file)
				print(*stats.get(chr_lengths), sep='\t', file=tsv_file)

			total_stats += stats

		if chromosome_count > 1:
			print('---------------- ALL chromosomes (aggregated) ----------------')
			total_stats.print(chr_lengths)
			if tsv_file:
				print(sample, 'ALL', args.vcf, sep='\t', end='\t', file=tsv_file)
				print(*total_stats.get(chr_lengths), sep='\t', file=tsv_file)
