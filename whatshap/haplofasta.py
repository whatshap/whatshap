"""
Write haplotypes in FASTA format from a phased VCF

This tool is provisional! Output file format and command-line
interface may change.

The phased variants from the input VCF are read and two new
references (haplotypes) are created by incorporating the
alleles into the reference genome. For each sample in the
input VCF, two FASTA files are written to the target directory.

With FOLDER being the target folder, the following files are written
for each SAMPLE:

- FOLDER/SAMPLE.1.fasta
- FOLDER/SAMPLE.2.fasta
- FOLDER/SAMPLE.1.liftover (TODO not a proper chain file)
- FOLDER/SAMPLE.2.liftover (TODO not a proper chain file)
"""
import sys
from collections import defaultdict
import logging

import pyfaidx

__author__ = "Alex Schoenhuth, Tobias Marschall, Marcel Martin"


logger = logging.getLogger(__name__)

allowed_dna_chars = set('ACGTNacgtn')


def valid_dna_string(s):
	chars = set(c for c in s)
	return chars.issubset(allowed_dna_chars)


def add(variants_dict, individual, chromosome, variant, genotype):
	if genotype == '.':
		pass
	if genotype == '0|0' or genotype == '0/0':
		pass
	elif genotype in ['1|0', '1|0']:
		variants_dict[(individual, chromosome, 1)].append(variant)
	elif genotype in ['0|1', '0|1']:
		variants_dict[(individual, chromosome, 2)].append(variant)
	elif genotype == '1|1' or genotype == '1/1':
		variants_dict[(individual, chromosome, 1)].append(variant)
		variants_dict[(individual, chromosome, 2)].append(variant)
	else:
		logger.error('Invalid genotype: %s', genotype)
		sys.exit(1)


def make_chromosome(chr_out, liftover_out, chromosome, reference, variants):
	logger.info('  %d variants', len(variants))
	modrefchrom = [x for x in reference]
	logger.info("Refchromlen: %s", len(modrefchrom))  # at the beginning modrefchrom is refchrom
	delcount = inscount = mixcount = snvcount = mnpcount = invcount = 0
	inversions = []
	run = 0
	for vartype, left, right, seq in variants:
		if vartype == 'SNV':
			snvcount += 1
			if modrefchrom[left] != '':
				modrefchrom[left] = modrefchrom[left][:-1] + seq
		elif vartype == 'MNP':
			mnpcount += 1
			modrefchrom[right] = seq + modrefchrom[right]
			run += len(seq)
			for i in range(left, right):
				if modrefchrom[i] != '':
					modrefchrom[i] = modrefchrom[i][:-1]
					run -= 1
		elif vartype == 'INS':
			inscount += 1
			modrefchrom[left] = seq + modrefchrom[left]
			run += len(seq)
		elif vartype == 'DEL':
			delcount += 1
			for i in range(left, right):
				if modrefchrom[i] != '':
					modrefchrom[i] = modrefchrom[i][:-1] 
					run -= 1
		elif vartype == 'MIX':
			mixcount += 1
			modrefchrom[right] = seq + modrefchrom[right]
			run += len(seq)
			for i in range(left, right):
				if modrefchrom[i] != '':
					modrefchrom[i] = modrefchrom[i][:-1]
					run -= 1
		elif vartype == 'INV':
			invcount += 1
			inversions.append((left,right))
		else:
			assert False
	for left, right in inversions:
		invseq = modrefchrom[left:right]
		invseq.reverse()
		modrefchrom[left:right] = invseq
	logger.info("SNV: %s MNP: %s MIX: %s DEL: %s INS: %s INV: %s",
		snvcount, mnpcount, mixcount, delcount, inscount, invcount)
	logger.info(run)
	# construct simchromstring
	simchrom = ''.join(modrefchrom)
	diff = len(simchrom) - len(modrefchrom) - run
	logger.info("Len simchrom: %s run: %s diff: %s", len(simchrom), run, diff)
	if diff != 0:
		logger.warning("diff not equal 0!")
	print(">%s" % (chromosome), file=chr_out)
	i = 0
	while i < len(simchrom):
		print(simchrom[i:i+50], file=chr_out)
		i += 50
	# write liftover file
	modind = 0
	numsame = 0
	diff = 0
	for ind, mer in enumerate(modrefchrom):
		for x in mer:
			if diff == modind - ind:
				numsame += 1
			else: # diff != modind - ind
				print(chromosome, numsame, diff, file=liftover_out)
				numsame = 1
				diff = modind - ind
			modind += 1
	if numsame > 0:
		print(chromosome, numsame, diff, file=liftover_out)


def run_haplofasta(vcf_path, reference_path, destination_folder, chromosome=None):
	"""

	:param vcf_path:
	:param reference_path:
	:param destination_folder:
	:param chromosome: Chromosome to process. If None, process all.
	:return:
	"""
	with pyfaidx.Fasta(reference_path, as_raw=True, sequence_always_upper=True) as fasta:
		if chromosome is not None:
			chromosomes = [chromosome]
		else:
			chromosomes = list(fasta.records)
		reference = {}
		for chromosome in chromosomes:
			reference[chromosome] = str(fasta[chromosome])
			logger.info('Loaded chromosome %r', chromosome)

	# read variants
	# mapping (individual, chromosome, allelenr) to lists of tuples (vartype, coord1, coord2, seq)
	variants = defaultdict(list)
	linenr = 0
	header = None
	individuals = None
	with open(vcf_path) as vcf:
		for line in vcf:
			line = line.strip()
			linenr += 1
			if line.startswith('##'):
				continue
			if line.startswith('#'):
				header = line[1:].split('\t')
				assert len(header) >= 10
				individuals = header[9:]
				continue
			assert header is not None
			fields = line.split('\t')
			assert len(fields) >= 10
			chrom = fields[0]
			variant_start = int(fields[1]) - 1
			variant_id = fields[2]
			variant_ref = fields[3]
			variant_alt = fields[4]
			# if chrom not in reference:
			# 	logger.warning('Skipping variant for unknown reference "%s" in line %d', chrom, linenr)
			# 	continue
			ref = reference[chrom]
			if variant_alt == '<INV>':
				# INVERSION
				if not valid_dna_string(variant_ref):
					logger.warning('Skipping invalid variant in line %s', linenr)
					continue
				inversion_start = variant_start
				inversion_end = variant_start+len(variant_ref)
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('INV', inversion_start, inversion_end, ''), genotype)
			elif (len(variant_ref) == 1) and (len(variant_alt) == 1):
				# SNV
				if not valid_dna_string(variant_ref) or not valid_dna_string(variant_alt):
					logger.warning('Skipping invalid variant in line %s', linenr)
					continue
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('SNV', variant_start, None, variant_alt), genotype)
			elif (len(variant_ref) > 1) and (len(variant_ref) == len(variant_alt)):
				# MNP
				if not valid_dna_string(variant_ref) or not valid_dna_string(variant_alt):
					logger.warning('Skipping invalid variant in line %s', linenr)
					continue
				while (len(variant_ref) > 0) and (len(variant_alt) > 0) and (variant_ref[0] == variant_alt[0]):
					variant_ref = variant_ref[1:]
					variant_alt = variant_alt[1:]
					variant_start += 1
				if len(variant_ref) == 0:
					continue
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('MNP', variant_start, variant_start+len(variant_ref), variant_alt), genotype)
			elif (len(variant_ref) > 1) and (len(variant_alt) == 1):
				# DELETION
				if not valid_dna_string(variant_ref) or not valid_dna_string(variant_alt):
					logger.warning('Skipping invalid variant in line', linenr)
					continue
				variant_end = variant_start + len(variant_ref)
				if variant_alt != variant_ref[0]:
					logger.error('ALT not equal to first character of REF in line %s', linenr)
					exit(1)
				del_start = variant_start + 1
				del_end = variant_end
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('DEL', del_start, del_end, ''), genotype)
			elif (len(variant_ref) == 1) and (len(variant_alt) > 1):
				# INSERTION
				if not valid_dna_string(variant_ref) or not valid_dna_string(variant_alt):
					logger.warning('Skipping invalid variant in line %s', linenr)
					continue
				if variant_alt[0] != variant_ref:
					logger.error('REF not equal to first character of ALT in line %s', linenr)
					exit(1)
				# position directly BEFORE breakpoint
				insertion_pos = variant_start + 1
				insertion_seq = variant_alt[1:]
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('INS', insertion_pos, None, insertion_seq), genotype)
			else:
				# MIX
				if not valid_dna_string(variant_ref) or not valid_dna_string(variant_alt):
					logger.warning('Skipping invalid variant in line %s', linenr)
					continue
				while (len(variant_ref) > 0) and (len(variant_alt) > 0) and (variant_ref[0] == variant_alt[0]):
					variant_ref = variant_ref[1:]
					variant_alt = variant_alt[1:]
					variant_start += 1
				for i,individual in enumerate(individuals):
					genotype = fields[9+i]
					add(variants, individual, chrom, ('MIX', variant_start, variant_start+len(variant_ref), variant_alt), genotype)
	logger.info('Read %s', vcf_path)

	# produce new alleles
	for individual in individuals:
		for allelenr in [1, 2]:
			logger.info('Processing individual %r, allele %s.', individual, allelenr)
			basepath = '{}/{}.{}'.format(destination_folder, individual, allelenr)
			with open(basepath + '.fasta', 'w') as chr_out, \
					open(basepath + '.liftover', 'w') as liftover_out:
				for chromosome in chromosomes:
					logger.info('... chromosome %r', chromosome)
					ref = reference[chromosome]
					make_chromosome(chr_out, liftover_out, chromosome, ref, variants[(individual, chromosome, allelenr)])


def add_arguments(parser):
	arg = parser.add_argument
	arg('-c', '--chromosome', action='store', default=None, help='Process only CHROMOSOME')
	arg('vcf', metavar='VCF', help='VCF file with phased variants (can be gzip-compressed)')
	arg('reference', metavar='FASTA', help='Full reference in FASTA format')
	arg('destination', metavar='FOLDER', help='Write output files to FOLDER/...')


def main(args):
	run_haplofasta(
		vcf_path=args.vcf,
		reference_path=args.reference,
		destination_folder=args.destination,
		chromosome=args.chromosome)


if __name__ == '__main__':
	# This script has currently no dependency on WhatsHap.
	# The following code makes it possible to run it stand-alone.
	logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
	from argparse import ArgumentParser
	parser = ArgumentParser(description=__doc__)
	add_arguments(parser)
	main(parser.parse_args())
