#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import sys
import vcf
from collections import defaultdict

# TODO: diffstats and evaluatephasings can be simplified a lot. We don't need to process "|" and "-" symbols any more

def diffstats(predhap, truehap, truehap_offset):
	"""true hap is supposed to be the true pair of haplotypes,
	gapless, that is everything is correctly phased"""

	#flips = 0
	switches = 0
	#switches = []
	phasings = [] # 'C' is correct, 'S' is switch error, 'U' is unphasable, 'HOM' is homozygosity error, 'HET' is heterozygosity error 
	homoerrors = 0 # truth: hetero, prediction: homo
	heteroerrors = 0 # truth: homo, prediction: hetero
	gaperrors = 0 # truth: phased, prediction: remains unphased
	breaks = 0 # number of "|" symbols, i.e. "unbridgeable" positions
	correct = 0 # number of correctly phased (pairs of neighboring) SNPs
	unphasable = 0 # number of (pairs of neighboring) SNPs that cannot be phased due to gaps or breaks
	ambiguous = 0
	
	oldx = ('-','-')
	truehap_pos = truehap_offset
	predhap0_pos = 0
	predhap1_pos = 0
	while True:
		# if break symbol found in both haplotypes, then increment the counter
		# if it is in only one, then the other haplotype has a gap which will be conted as such
		if (predhap[0][predhap0_pos] == '|') and (predhap[1][predhap1_pos] == '|'):
			breaks += 1
		found_break = False
		if predhap[0][predhap0_pos] == '|':
			found_break = True
			predhap0_pos += 1
		if predhap[1][predhap1_pos] == '|':
			found_break = True
			predhap1_pos += 1
		if found_break:
			oldx = ('-','-')
			continue
	
		if (oldx[0] == '-') or (oldx[1] == '-') or (predhap[0][predhap0_pos] == '-') or (predhap[1][predhap1_pos] == '-'):
			unphasable += 1
			phasings.append('U')
			if (predhap[0][predhap0_pos] == '-') or (predhap[1][predhap1_pos] == '-'): 
				# gap! In our data, this means SNP desert, so forget about switch errors across a gap.
				gaperrors += 1
				oldx = ('-','-')
			else:
				oldx = (predhap[0][predhap0_pos], truehap[0][truehap_pos])

		elif (oldx[0] == 'X') or (oldx[1] == 'X') or (predhap[0][predhap0_pos] == 'X') or (predhap[1][predhap1_pos] == 'X'):
			ambiguous += 1
			phasings.append('A')
			if (predhap[0][predhap0_pos] == 'X') or (predhap[1][predhap1_pos] == 'X'): 
				# gap! In our data, this means SNP desert, so forget about switch errors across a gap.
				# gaperrors += 1
				oldx = ('X','X')
			else:
				oldx = (predhap[0][predhap0_pos], truehap[0][truehap_pos])

		elif (predhap[0][predhap0_pos] == predhap[1][predhap1_pos]) and (truehap[0][truehap_pos] != truehap[1][truehap_pos]): 
			# truth: hetero, prediction: homo, leave oldx!
			phasings.append('O')
			homoerrors += 1
			#print("Homo-erro at SNP: ", truehap_pos)

		elif (predhap[0][predhap0_pos] != predhap[1][predhap1_pos]) and (truehap[0][truehap_pos] == truehap[1][truehap_pos]): 
			# truth: homo, prediction: hetero, leave oldx!
			phasings.append('E')
			heteroerrors += 1
			#print("Hetero-error at SNP: ", truehap_pos)
		# NOTE: the below conditions lead to that switch errors are only counted within stretches of predictions
		#       i.e. the first position after a gap is never considered to be switch error
		elif ((predhap[0][predhap0_pos] == oldx[0] and truehap[0][truehap_pos] != oldx[1]) or 
		      (predhap[0][predhap0_pos] != oldx[0] and truehap[0][truehap_pos] == oldx[1])):
			switches += 1
			phasings.append('S')
			#print("Switch at SNP: ", truehap_pos)
			oldx = (predhap[0][predhap0_pos], truehap[0][truehap_pos])
		else:
			# No error detected
			oldx = (predhap[0][predhap0_pos], truehap[0][truehap_pos])
			correct += 1
			phasings.append('C')

		truehap_pos += 1
		predhap0_pos += 1
		predhap1_pos += 1

		if predhap0_pos == len(predhap[0]) : break
		if predhap1_pos == len(predhap[1]) : break
		if truehap_pos == len(truehap[0]) : break

	return switches, homoerrors, heteroerrors, gaperrors, breaks, unphasable, ambiguous, correct, phasings

def evaluatephasings(phasings):

	"""expects a list with phasing information ['C' is for
	correct, 'S' is for switch error, 'U' is for unphasable, 'O'
	is for homozygosity error, 'E' is for heterozygosity error]
	and outputs statistics on flip, switch, homo- and
	heterozygosity errors, as well as unphasable positions"""

	flips = 0
	switches = 0
	deserts = 0
	ambiguous = 0
	homo = 0
	hetero = 0
	correct = 0
	oldswitches = 0

#	nozygosity = []
#	for phasing in phasings:
#		if phasing == 'O':
#			homo += 1
#		elif phasing == 'E':
#			hetero += 1
#		else:
#			nozygosity.append(phasing)
	
	
	consecutive = 0
	unphasable = True
	for i, phasing in enumerate(phasings):
		if phasing == 'O' or phasing == 'E':
			if phasing == 'O':
				homo += 1
			if phasing == 'E':
				hetero += 1
			continue
		if phasing == 'A':
			ambiguous += 1
			continue
		if phasing == 'C':
			#print(int(consecutive/2))
			if unphasable:
				flips += (consecutive+1) // 2
				oldswitches += consecutive
			elif not unphasable:
				switches += consecutive % 2
				flips += consecutive // 2
				oldswitches += consecutive
			correct += 1
			consecutive = 0
			unphasable = False

		elif phasing == 'U':
			flips += (consecutive+1) // 2
			oldswitches += consecutive
			deserts += 1
			consecutive = 0
			unphasable = True
		
		elif phasing == 'S':
			consecutive += 1
			#oldswitches += 1
			
	return flips, switches, homo, hetero, deserts, ambiguous, correct, oldswitches

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

def read_phased_vcf(filename, true_snps = None):
	"""Reads VCF as output by whatshap and returns a tuple: (phased, unphased, blocks), where
	"phased" and "unphased" give the number of phased/unphased heterozygous SNPs and "blocks"
	is a dict with one entry for each phased block, mapping (chromosome,block-ID) --> [(pos,gt),..].
	If true_snps is given, then only entries that are also in true snps are retained and FPs and FNs
	are computed."""
	vcf_reader = vcf.Reader(filename)
	samples = vcf_reader.samples
	assert len(samples) == 1, 'Expected exactly one sample in phased VCF'
	blocks = defaultdict(list)
	phased = 0
	unphased = 0
	tp_snps = set()
	FP = 0
	FN = 0
	for record in vcf_reader:
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			continue
		call = record.samples[0]
		if not call.is_het:
			continue
		if true_snps != None:
			if (record.CHROM, record.POS) not in true_snps:
				FP += 1
				continue
		tp_snps.add( (record.CHROM, record.POS) )
		if hasattr(call.data,'HP'):
			phased += 1
			block_id, genotype = parse_hp(call.data.HP)
			blocks[(record.CHROM, block_id)].append((record.POS,genotype))
		else:
			unphased += 1
	if true_snps != None:
		FN = len(true_snps) - len(tp_snps)
	return phased, unphased, FP, FN, blocks

def main():
	parser = ArgumentParser(prog='evaluate-phasing', description=__doc__)
	parser.add_argument('-v', dest='verbose', action='store_true', default=False,
		help='Be (very) verbose and output statistics on every single phased block.')
	parser.add_argument('truthvcf', metavar='truthvcf',
		help='VCF with true phasing')
	parser.add_argument('predictedvcf', metavar='predictedvcf',
		help='VCF with predicted phasing as output by Whatshap')

	args = parser.parse_args()

	truth = read_truth_vcf(open(args.truthvcf))
	
	print('Read {} heterozygous SNPs from file with true phasing'.format(len(truth)), file=sys.stderr)
	
	phased, unphased, FP, FN, phased_blocks = read_phased_vcf(open(args.predictedvcf), truth)
	
	print('Found {} false positive (FP) heterozygous SNPs, i.e. heterozygous in predictions but absent or homozygous in truth'.format(FP), file=sys.stderr)
	print('There are {} false negative (FN) heterozygous SNPs, i.e. absent or homozygous in predictions but heterozygous in truth'.format(FN), file=sys.stderr)
	print('Retained {} heterozygous SNPs from file with predicted phasing, out of which:'.format(phased+unphased), file=sys.stderr)
	
	print('  phased SNPs: {}'.format(phased), file=sys.stderr)
	print('  unphased SNPs: {}'.format(unphased), file=sys.stderr)
	print('  number of blocks: {}'.format(len(phased_blocks)), file=sys.stderr)
	print('    --> unphasable due to being first in a block: {}'.format(len(phased_blocks)), file=sys.stderr)
	
	block_list = list(phased_blocks.keys())
	block_list.sort()
	correct_total = 0
	switches_total = 0
	flips_total = 0
	homoerrors_total = 0
	heteroerrors_total = 0
	blocks_larger_1 = 0
	for (chromosome,block_id) in block_list:
		block = phased_blocks[(chromosome,block_id)]
		if len(block) >= 2:
			blocks_larger_1 += 1
		truehap = ['','']
		predhap = ['','']
		for pos, pred_gt in block:
			predhap[0] += pred_gt[0]
			predhap[1] += pred_gt[2]
			assert (chromosome, pos) in truth, 'Positions in truth VCF and predictions VCF do not match'
			true_gt = truth[(chromosome, pos)]
			truehap[0] += true_gt[0]
			truehap[1] += true_gt[2]
		switches0, homoerrors0, heteroerrors0, gaperrors0, breaks0, unphasable0, ambiguous0, correct0, phasings = diffstats(predhap, truehap, 0)
		flips, switches1, homoerrors1, heteroerrors1, unphasable1, ambiguous1, correct1, oldswitches = evaluatephasings(phasings)
		assert unphasable1 == 1
		assert ambiguous1 == 0
		if args.verbose:
			print('-'*100, file=sys.stderr)
			print('Block ID:', block_id, file=sys.stderr)
			print('True block:     ', truehap[0], file=sys.stderr)
			print('Predicted block:', predhap[0], file=sys.stderr)
			print("Phasings:", ''.join(phasings), file=sys.stderr)
			print("Correctly phased: ", correct1, file=sys.stderr)
			print("Unphasable (due to gaps or breaks): ", unphasable1, file=sys.stderr)
			print("Switch errors: ", switches1, file=sys.stderr)
			print("Flip errors: ", flips, file=sys.stderr)
			#print("Truth hetero, but predicted homo: ", homoerrors1)
			#print("Truth homo, but predicted hetero: ", heteroerrors1)
		correct_total += correct1
		switches_total += switches1
		flips_total += flips
		homoerrors_total += homoerrors1
		heteroerrors_total += heteroerrors1
	print('='*100, file=sys.stderr)
	print('Evaluation of blocks:', file=sys.stderr)
	print('Number of blocks with at least 2 SNPs:', blocks_larger_1, file=sys.stderr)
	print("Correctly phased: ", correct_total, file=sys.stderr)
	print("Switch errors: ", switches_total, file=sys.stderr)
	print("Flip errors: ", flips_total, file=sys.stderr)
	assert homoerrors_total == heteroerrors_total == 0
	print('#filename', 'phased', 'unphased', ' len(phased_blocks)', 'blocks_larger_1', 'correct_total', 'switches_total', 'flips_total')
	print(args.predictedvcf, phased, unphased,  len(phased_blocks), blocks_larger_1, correct_total, switches_total, flips_total)
	#print("Truth hetero, but predicted homo: ", homoerrors_total)
	#print("Truth homo, but predicted hetero: ", heteroerrors_total)

if __name__ == '__main__':
	sys.exit(main())
