"""
Print phasing statistics
"""
import logging
from collections import defaultdict, namedtuple
from .math import median
from .vcf import VcfReader

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('--gtf', default=None, help='Write phased blocks to GTF file.')
	add('--sample', metavar='SAMPLE', default=None, help='Name of the sample '
			'to process. If not given, use first sample found in VCF.')
	add('--tsv', metavar='TSV', default=None, help='Filename to write '
		'statistics to (tab-separated).')
	add('vcf', metavar='VCF', help='Phased VCF file')


class PhasedBlock:
	def __init__(self):
		self.phases = {}
		self.leftmost_variant = None
		self.rightmost_variant = None

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

	def variants():
		return list(sorted(self.phases.keys()))

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

detailed_stats_fields = ['variants','phased','unphased','singletons','blocks',
	'variant_per_block_median','variant_per_block_avg','variant_per_block_min','variant_per_block_max','variant_per_block_sum',
	'bp_per_block_median','bp_per_block_avg','bp_per_block_min','bp_per_block_max','bp_per_block_sum'
]
DetailedStats = namedtuple('DetailedStats', detailed_stats_fields)


class PhasingStats:
	def __init__(self):
		self.blocks = []
		self.unphased = 0
		self.variants = 0

	def __iadd__(self, other):
		self.blocks.extend(other.blocks)
		self.unphased += other.unphased
		self.variants += other.variants
		return self

	def add_blocks(self, blocks):
		self.blocks.extend(blocks)

	def add_unphased(self, unphased=1):
		self.unphased += unphased

	def add_variants(self, variants):
		self.variants += variants

	def get(self):
		block_sizes = sorted(len(block) for block in self.blocks)
		n_singletons = sum(1 for size in block_sizes if size == 1)
		block_sizes = [ size for size in block_sizes if size > 1 ]
		block_lengths = sorted(block.span() for block in self.blocks if len(block) > 1)
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
				bp_per_block_sum = sum(block_lengths)
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
				variant_per_block_min = float('nan'),
				variant_per_block_max = float('nan'),
				variant_per_block_sum = 0,
				bp_per_block_median = float('nan'),
				bp_per_block_avg = float('nan'),
				bp_per_block_min = float('nan'),
				bp_per_block_max = float('nan'),
				bp_per_block_sum = 0
			)

	def print(self):
		stats = self.get()

		WIDTH = 21
		print('Variants in VCF:'.rjust(WIDTH), '{:8d}'.format(stats.variants))
		print('Phased:'.rjust(WIDTH), '{:8d}'.format(stats.phased))
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


def main(args):
	gtfwriter = None
	if args.gtf:
		gtf_file = open(args.gtf, 'wt')
		gtfwriter = GtfWriter(gtf_file)
	if args.tsv:
		tsv_file = open(args.tsv, 'w')
	else:
		tsv_file = None

	vcf_reader = VcfReader(args.vcf, indels=False)  # TODO: also indels
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
		print('#sample', 'chromosome', 'file_name', sep='\t', end='\t', file=tsv_file)
		print(*detailed_stats_fields, sep='\t', file=tsv_file)

	print('Phasing statistics for sample {} from file {}'.format(sample, args.vcf))
	total_stats = PhasingStats()
	chromosome_count = 0
	for variant_table in vcf_reader:
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
			if phase is None:
				stats.add_unphased()
			else:
				blocks[phase.block_id].add(variant, phase)
				if gtfwriter:
					if prev_block_id is None:
						prev_block_fragment_start = variant.position
						prev_block_fragment_end = variant.position + 1
						prev_block_id = phase.block_id
					else:
						if (prev_block_id != phase.block_id):
							gtfwriter.write(chromosome, prev_block_fragment_start, prev_block_fragment_end, prev_block_id)
							prev_block_fragment_start = variant.position
							prev_block_id = phase.block_id
						prev_block_fragment_end = variant.position + 1
		if gtfwriter and (not prev_block_id is None):
			gtfwriter.write(chromosome, prev_block_fragment_start, prev_block_fragment_end, prev_block_id)

		stats.add_blocks(blocks.values())
		stats.print()
		if tsv_file:
			print(sample, chromosome, args.vcf, sep='\t', end='\t', file=tsv_file)
			print(*stats.get(), sep='\t', file=tsv_file)

		total_stats += stats

	if chromosome_count > 1:
		print('---------------- ALL chromosomes (aggregated) ----------------'.format(chromosome))
		total_stats.print()
		if tsv_file:
			print(sample, 'ALL', args.vcf, sep='\t', end='\t', file=tsv_file)
			print(*total_stats.get(), sep='\t', file=tsv_file)

	if gtfwriter:
		gtf_file.close()

	if tsv_file:
		tsv_file.close()
