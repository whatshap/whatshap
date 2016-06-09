"""
Compare two or more phasings
"""
import sys
import logging
from collections import defaultdict, namedtuple
from .vcf import VcfReader

logger = logging.getLogger(__name__)

count_width = 9

def add_arguments(parser):
	add = parser.add_argument
	add('--sample', metavar='SAMPLE', default=None, help='Name of the sample '
			'to process. If not given, use first sample found in VCF.')
	add('--names', metavar='NAMES', default=None, help='Comma-separated list '
			'of data set names to be used in the report (in same order as VCFs).')
	# TODO: what's the best way to request "two or more" VCFs?
	add('vcf', nargs='+', metavar='VCF', help='At least two phased VCF files to be compared.')


def validate(args, parser):
	if len(args.vcf) < 2:
		parser.error('At least two VCFs need to be given.')


class SwitchFlips:
	def __init__(self, switches=0, flips=0):
		self.switches = switches
		self.flips = flips
	def __iadd__(self, other):
		self.switches += other.switches
		self.flips += other.flips
		return self
	def __repr__(self):
		return 'SwitchFlips(switches={}, flips={})'.format(self.switches, self.flips)


class PhasingErrors:
	def __init__(self, switches=0, hamming=0, switch_flips=None):
		self.switches = switches
		self.hamming = hamming
		self.switch_flips = SwitchFlips() if switch_flips is None else switch_flips
	def __iadd__(self, other):
		self.switches += other.switches
		self.hamming += other.hamming
		self.switch_flips += other.switch_flips
		return self
	def __repr__(self):
		return 'PhasingErrors(switches={}, hamming={}, switch_flips={})'.format(self.switches, self.hamming, self.switch_flips)


def complement(s):
	t = { '0': '1', '1':'0' }
	return ''.join(t[c] for c in s)


def hamming(s0, s1):
	assert len(s0) == len(s1)
	return sum( c0!=c1 for c0, c1 in zip(s0, s1) )


def switch_encoding(phasing):
	return ''.join( ('0' if phasing[i-1]==phasing[i] else '1') for i in range(1,len(phasing)) )


def compute_switch_flips(phasing0, phasing1):
	assert len(phasing0) == len(phasing1)
	s0 = switch_encoding(phasing0)
	s1 = switch_encoding(phasing1)
	result = SwitchFlips()
	switches_in_a_row = 0
	for i, (p0, p1) in enumerate(zip(s0, s1)):
		if p0 != p1:
			switches_in_a_row += 1
		if (i + 1 == len(s0)) or (p0 == p1):
			result.flips += switches_in_a_row // 2
			result.switches += switches_in_a_row % 2
			switches_in_a_row = 0
	if False:
		print('switch_flips():')
		print('   phasing0={}'.format(phasing0))
		print('   phasing1={}'.format(phasing1))
		print('         s0={}'.format(s0))
		print('         s1={}'.format(s1))
		print('   switches={}, flips={}'.format(result.switches, result.flips))
	return result


def compare_block(phasing0, phasing1):
	"""Input are two strings over {0,1}. Output is a PhasingErrors object."""
	return PhasingErrors(
		switches = hamming(switch_encoding(phasing0), switch_encoding(phasing1)),
		hamming = min(hamming(phasing0, phasing1), hamming(phasing0, complement(phasing1))),
		switch_flips = compute_switch_flips(phasing0, phasing1)
	)


def fraction2str(nominator, denominator):
	if denominator == 0:
		return '--'
	else:
		return '{:.2f}%'.format(nominator*100.0/denominator)


def print_errors(errors, phased_pairs, print_hamming=False):
	print('    phased pairs of variants assessed:', str(phased_pairs).rjust(count_width))
	print('                        switch errors:', str(errors.switches).rjust(count_width))
	print('                    switch error rate:', fraction2str(errors.switches, phased_pairs).rjust(count_width))
	print('            switch/flip decomposition:', '{}/{}'.format(errors.switch_flips.switches,errors.switch_flips.flips).rjust(count_width) )
	print('                     switch/flip rate:', fraction2str(errors.switch_flips.switches+errors.switch_flips.flips, phased_pairs).rjust(count_width))


def compare(variant_tables, sample, dataset_names):
	assert len(variant_tables) > 1
	variants = None
	for variant_table in variant_tables:
		het_variants = [ v for v, gt in zip(variant_table.variants, variant_table.genotypes_of(sample)) if gt == 1 ]
		if variants is None:
			variants = set(het_variants)
		else:
			variants.intersection_update(het_variants)
	print('         common heterozygous variants:', str(len(variants)).rjust(count_width))
	print('         (restricting to these below)')
	phases = []
	for variant_table in variant_tables:
		p = [ phase for variant, phase in zip(variant_table.variants, variant_table.phases_of(sample)) if variant in variants ]
		assert len(p) == len(variants)
		phases.append(p)

	blocks = [ defaultdict(list) for _ in variant_tables ]
	block_intersection = defaultdict(list)
	for variant_index in range(len(variants)):
		any_none = False
		for i in range(len(phases)):
			phase = phases[i][variant_index]
			if phase is None:
				any_none = True
			else:
				blocks[i][phase.block_id].append(variant_index)
		if not any_none:
			joint_block_id = tuple( phases[i][variant_index].block_id for i in range(len(phases)) )
			block_intersection[joint_block_id].append(variant_index)
	for i in range(len(phases)):
		print('non-singleton blocks in {}:'.format(dataset_names[i]).rjust(38), str(len([b for b in blocks[i].values() if len(b) > 1])).rjust(count_width))
		print('                 --> covered variants:', str(sum(len(b) for b in blocks[i].values() if len(b) > 1)).rjust(count_width))
	print('    non-singleton intersection blocks:', str(len([b for b in block_intersection.values() if len(b) > 1])).rjust(count_width))
	print('                 --> covered variants:', str(sum(len(b) for b in block_intersection.values() if len(b) > 1)).rjust(count_width))
	longest_block = None
	longest_block_errors = None
	phased_pairs = 0
	if len(variant_tables) == 2:
		total_errors = PhasingErrors()
		for block in block_intersection.values():
			if len(block) < 2:
				continue
			phasing0 = ''.join( str(phases[0][i].phase) for i in block )
			phasing1 = ''.join( str(phases[1][i].phase) for i in block )
			errors = compare_block(phasing0, phasing1)
			total_errors += errors
			phased_pairs += len(block) - 1
			if (longest_block is None) or (len(block) > longest_block):
				longest_block = len(block)
				longest_block_errors = errors
		print('              ALL INTERSECTION BLOCKS:', '-'*count_width)
		print_errors(total_errors, phased_pairs)
		print('           LARGEST INTERSECTION BLOCK:', '-'*count_width)
		print_errors(longest_block_errors, longest_block-1)
		print('                     Hamming distance:', str(longest_block_errors.hamming).rjust(count_width))
		print('                 Hamming distance [%]:', fraction2str(longest_block_errors.hamming, longest_block).rjust(count_width))
	else:
		histogram = defaultdict(int)
		total_compared = 0
		for block in block_intersection.values():
			if len(block) < 2:
				continue
			total_compared += len(block) - 1
			phasings = [ ''.join( str(phases[j][i].phase) for i in block ) for j in range(len(phases)) ]
			switch_encodings = [ switch_encoding(p) for p in phasings ]
			for i in range(len(block)-1):
				s = ''.join(switch_encodings[j][i] for j in range(len(switch_encodings)) )
				s = min(s, complement(s))
				histogram[s] += 1
		print('           Compared pairs of variants:', str(total_compared).rjust(count_width))
		bipartitions = list(histogram.keys())
		bipartitions.sort()
		for i, s in enumerate(bipartitions):
			count = histogram[s]
			if i == 0:
				assert set(c for c in s) == set('0')
				print('                            ALL AGREE:')
			elif i == 1:
				print('                         DISAGREEMENT:')
			left, right = [], []
			for name, leftright in zip(dataset_names,s):
				if leftright == '0':
					left.append(name)
				else:
					right.append(name)
			print(
				('{%s} vs. {%s}:' % (','.join(left),','.join(right))).rjust(38),
				str(count).rjust(count_width),
				fraction2str(count, total_compared).rjust(8)
			)


def main(args):
	vcf_readers = [VcfReader(f, indels=False) for f in args.vcf]  # TODO: also indels
	if args.names:
		dataset_names = args.names.split(',')
		if len(dataset_names) != len(args.vcf):
			logger.error('Number of names given with --names does not equal number of VCFs.')
			sys.exit(1)
	else:
		dataset_names = ['file{}'.format(i) for i in range(len(args.vcf))]
	longest_name = max(len(n) for n in dataset_names)

	all_samples = set()
	sample_intersection = None
	for vcf_reader in vcf_readers:
		if sample_intersection is None:
			sample_intersection = set(vcf_reader.samples)
		else:
			sample_intersection.intersection_update(vcf_reader.samples)
		all_samples.update(vcf_reader.samples)

	if args.sample:
		sample_intersection.intersection_update([args.sample])
		if len(sample_intersection) == 0:
			logger.error('Sample %r requested on command-line not found in all VCFs', args.sample)
			sys.exit(1)
		sample = args.sample
	else:
		if len(sample_intersection) == 0:
			logger.error('None of the samples is present in all VCFs')
			sys.exit(1)
		elif len(sample_intersection) == 1:
			sample = list(sample_intersection)[0]
		else:
			logger.error('More than one sample is present in all VCFs, please use --sample to specify which sample to work on.')
			sys.exit(1)
	print('Comparing phasings for sample', sample)

	chromosomes = None
	vcfs = []
	for reader, filename in zip(vcf_readers, args.vcf):
		# create dict mapping chromsome names to VariantTables
		m = dict()
		logger.info('Reading phasing from %r', filename)
		for variant_table in reader:
			m[variant_table.chromosome] = variant_table
		vcfs.append(m)
		if chromosomes is None:
			chromosomes = set(m.keys())
		else:
			chromosomes.intersection_update(m.keys())
	if len(chromosomes) == 0:
		logger.error('No chromosome is contained in all VCFs. Aborting.')
		sys.exit(1)

	logger.info('Chromosomes present in all VCFs: %s', ', '.join(sorted(chromosomes)))

	print('FILENAMES')
	for name, filename in zip(dataset_names, args.vcf):
		print(name.rjust(longest_name+2), '=', filename)

	width = max(longest_name, 15) + 5
	for chromosome in sorted(chromosomes):
		print('---------------- Chromosome {} ----------------'.format(chromosome))
		variant_tables = [ vcf[chromosome] for vcf in vcfs ]
		all_variants_union = set()
		all_variants_intersection = None
		het_variants_union = set()
		het_variants_intersection = None
		het_variant_sets = []
		print('VARIANT COUNTS (heterozygous / all): ')
		for variant_table, name in zip(variant_tables, dataset_names):
			all_variants_union.update(variant_table.variants)
			het_variants = [ v for v, gt in zip(variant_table.variants, variant_table.genotypes_of(sample)) if gt == 1 ]
			het_variants_union.update(het_variants)
			if all_variants_intersection is None:
				all_variants_intersection = set(variant_table.variants)
				het_variants_intersection = set(het_variants)
			else:
				all_variants_intersection.intersection_update(variant_table.variants)
				het_variants_intersection.intersection_update(het_variants)
			het_variant_sets.append(set(het_variants))
			print('{}:'.format(name).rjust(width), str(len(het_variants)).rjust(count_width), '/', str(len(variant_table.variants)).rjust(count_width))
		print('UNION:'.rjust(width), str(len(het_variants_union)).rjust(count_width), '/', str(len(all_variants_union)).rjust(count_width))
		print('INTERSECTION:'.rjust(width), str(len(het_variants_intersection)).rjust(count_width), '/', str(len(all_variants_intersection)).rjust(count_width))

		for i in range(len(vcfs)):
			for j in range(i+1, len(vcfs)):
				print('PAIRWISE COMPARISON: {} <--> {}:'.format(dataset_names[i],dataset_names[j]))
				compare([variant_tables[i], variant_tables[j]], sample, dataset_names)
		
		if len(vcfs) > 2:
			print('MULTIWAY COMPARISON OF ALL PHASINGS:')
			compare(variant_tables, sample, dataset_names)
