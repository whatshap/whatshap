#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys

usage = """%prog [options] <predicted.haplo> <true.haplo>

Compares a predicted to a true haplotype. The predicted haplotype may contain
gap ("-") or break ("|") symbols while the true haplotype is expected not to."""

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
		
def main():
	parser = OptionParser(usage=usage)

	(options, args) = parser.parse_args()
	if len(args) != 2:
		parser.print_help()
		sys.exit(1)

	hap1 = open(args[0], 'r').readlines()
	hap2 = open(args[1], 'r').readlines()

	hap1 = [x.strip() for x in hap1]
	hap2 = [x.strip() for x in hap2]

	assert len(hap1[0].replace('|','')) == len(hap1[1].replace('|',''))
	assert len(hap2[0]) == len(hap2[1])

	print(len(hap1[0].replace('|','')), len(hap2[0]))

	# again, we need this offset in case predhap starts further down than first pos
	#start_pos = int(open(args[2],'r').readline().split()[0])
	#positions = [int(x) for x in open(args[3],'r')]
	offset = 0
	#if len(hap1[0].replace('|','')) < len(positions) : # haplo is shorter
		#for p in positions :
			#if p == start_pos : break
			#offset += 1

	switches0, homoerrors0, heteroerrors0, gaperrors0, breaks0, unphasable0, ambiguous0, correct0, phasings = diffstats(hap1, hap2, offset)
	#phasings = diffstats(hap1, hap2, offset)[-1]
	
	flips, switches1, homoerrors1, heteroerrors1, unphasable1, ambiguous1, correct1, oldswitches = evaluatephasings(phasings)
	# I think this is a more meaningful definition of fragments, no?
	fragments1 = filter(lambda x:len(x)>0,hap1[0].replace('|','-').split('-'))
	fragments2 = filter(lambda x:len(x)>0,hap1[1].replace('|','-').split('-'))
	# was :
	#fragments1 = hap1[0].replace('|','-').split('-')
	#fragments2 = hap1[1].replace('|','-').split('-')
	# which is not meaningful, i.e., '00--0-'.split('-') = ['00', '', '0', '']
	print("Number of fragments in predicted haplotype1:", len(fragments1))
	print("Number of fragments in predicted haplotype2:", len(fragments2))
	print("Average fragment length in predicted haplotype1:", float(sum(len(f) for f in fragments1)) / float(len(fragments1)))
	print("Average fragment length in predicted haplotype2:", float(sum(len(f) for f in fragments2)) / float(len(fragments2)))
	print("Breaks: ", breaks0)
	print("Gaps in predicted haplotype: ", gaperrors0)
	print("======= Breakdown of SNPs into categories =======")
	print("Total SNP positions: ", len(hap2[0]))
	print("Correctly phased: ", correct1)
	print("Unphasable (due to gaps or breaks): ", unphasable1)
	print("Ambiguous (due to equal amounts of 0's and 1's in parts of partition): ", ambiguous1)
	print("Switch errors: ", switches1)
	print("Flip errors: ", flips)
	print("Truth hetero, but predicted homo: ", homoerrors1)
	print("Truth homo, but predicted hetero: ", heteroerrors1) 
	
	

if __name__ == '__main__':
	sys.exit(main())
