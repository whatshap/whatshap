#!/usr/bin/env python3

from argparse import ArgumentParser
import sys
import vcf
from collections import defaultdict


def read_truth_vcf(filename):
	"""Reads VCF with true phasing and returns a dict mapping (chr,pos) --> genotype. That 
	dict will contain one entry per heterozygous SNP."""
	vcf_reader = vcf.Reader(filename)
	samples = vcf_reader.samples
	assert len(samples) == 1, 'Expected exactly one sample in truth VCF'
	result = {}
	for record in vcf_reader:
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			continue
		call = record.samples[0]
		if not call.is_het:
			continue
		genotype = call.data.GT
		assert genotype in ['0|1', '1|0'], 'Expecting all SNPs to be phased (in "|" notation) in truth file'
		result[(record.CHROM, record.POS)] = genotype
	return result


def parse_hp(HP):
	"""Parse HP field from VCF record and returns (id,genotype), where
	genotype is in "|"-notation."""
	assert len(HP) == 2
	fields = [[int(x) for x in s.split('-')] for s in HP]
	assert fields[0][0] == fields[1][0]
	block_id = fields[0][0]
	genotype = '{}|{}'.format(fields[0][1]-1,fields[1][1]-1)
	assert genotype in ['0|1', '1|0']
	return block_id, genotype


def read_phased_vcf(filename, target_snps = None):
	"""Reads VCF as output by whatshap and returns a tuple: (phased, unphased, blocks), where
	"phased" and "unphased" give the number of phased/unphased heterozygous SNPs and "blocks"
	is a dict with one entry for each phased block, mapping (chromosome,block-ID) --> [(pos,gt),..].
	If true_snps is given, then only entries that are also in true snps are retained and FPs and FNs
	are computed."""
	print('Reading VCF:', filename, file=sys.stderr)
	vcf_reader = vcf.Reader(open(filename))
	samples = vcf_reader.samples
	chromosomes = set()
	snps = set()
	assert len(samples) == 1, 'Expected exactly one sample in phased VCF, found {}'.format(str(samples))
	blocks = defaultdict(list)
	global_phasing = not ('HP' in vcf_reader.formats)
	if global_phasing:
		print('  ... no HP tags, assuming global phasing given through genotypes {0|1, 1|0}.', file = sys.stderr)
	else:
		print('  ... block-wise phasing in HP tags.', file = sys.stderr)
	for record in vcf_reader:
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			continue
		call = record.samples[0]
		if not call.is_het:
			continue
		chromosomes.add(record.CHROM)
		genotype = call.data.GT
		snp = (record.CHROM, record.POS)
		if target_snps != None:
			if not snp in target_snps:
				continue
		snps.add( snp )
		if global_phasing:
			assert genotype in ['0|1', '1|0'], 'Expecting all SNPs to be phased (in "|" notation) in truth file'
			blocks[(record.CHROM,0)].append((record.POS,genotype))
		else:
			if hasattr(call.data,'HP'):
				block_id, genotype = parse_hp(call.data.HP)
				blocks[(record.CHROM, block_id)].append((record.POS,genotype))
	non_singleton_blocks = 0
	phased_variants = 0
	for (chromosome, block_id), block in blocks.items():
		if len(block) >= 2:
			non_singleton_blocks += 1
			phased_variants += len(block)
	print('  ... read {} heterozygous variants'.format(len(snps)), file = sys.stderr)
	print('  ... found {} non-singleton blocks (containing >=2 variants)'.format(non_singleton_blocks), file = sys.stderr)
	print('  ... {} variants are in a non-singleton block (and are thus phased)'.format(phased_variants), file = sys.stderr)
	return samples[0], chromosomes, snps, blocks


def create_phase_dict(blocks):
	result = {}
	for (chromosome, block_id), block in blocks.items():
		for i in range(1, len(block)):
			pos0, genotype0 = block[i-1]
			pos1, genotype1 = block[i]
			result[(pos0,pos1)] = (genotype0 == genotype1)
	return result


def main():
	parser = ArgumentParser(prog='evaluate-phasing', description=__doc__)
	#parser.add_argument('-v', dest='verbose', action='store_true', default=False,
		#help='Be (very) verbose and output statistics on every single phased block.')
	#parser.add_argument('truthvcf', metavar='truthvcf',
		#help='VCF with true phasing')
	parser.add_argument('vcf', metavar='vcf', nargs='+',
		help='VCF files (at least two) with phasings to be compared.')

	args = parser.parse_args()
	
	if len(args.vcf) < 2:
		print('Expected at least two VCF files as input.', file=sys.stderr)
		return 1

	datasets = [read_phased_vcf(filename) for filename in args.vcf]
	
	all_samples = set()
	all_chromosomes = set()
	all_snps = set()
	snp_intersection = None
	for i, (sample, chromosomes, snps, blocks) in enumerate(datasets):
		if i == 0:
			snp_intersection = set(snps)
		else:
			snp_intersection.intersection_update(snps)
		all_snps.update(snps)
		all_samples.add(sample)
		all_chromosomes.update(chromosomes)
	if len(all_samples) > 1:
		print('All VCF files need to contain data on the same (one) sample, but found samples', all_samples, file=sys.stderr)
		return 1
	if len(all_chromosomes) > 1:
		print('All VCF files need to contain data on the same (one) chromosome, but found chromosomes', all_chromosomes, file=sys.stderr)
		return 1

	print('Found {} different SNPs across all datasets'.format(len(all_snps)), file=sys.stderr)
	print(' ... out of which {} are present in all datasets'.format(len(snp_intersection)), file=sys.stderr)
	if all_snps != snp_intersection:
		print('Re-reading all data sets, restricting them to the common SNPs', file=sys.stderr)
		datasets = [read_phased_vcf(filename, snp_intersection) for filename in args.vcf]

	phases = [ create_phase_dict(blocks) for (sample, chromosomes, snps, blocks) in datasets ]
	
	positions = [pos for (chromosome, pos) in snp_intersection]
	positions.sort()
	
	histogram = defaultdict(int)
	for i in range(1, len(positions)):
		pos0 = positions[i-1]
		pos1 = positions[i]
		s = ''
		for phase_dict in phases:
			if (pos0,pos1) in phase_dict:
				s += str(int(phase_dict[(pos0,pos1)]))
			else:
				s += '-'
		histogram[s] += 1
		#print(pos0, pos1, s)

	for s, count in histogram.items():
		print(s, count)

	return 0


if __name__ == '__main__':
	sys.exit(main())
