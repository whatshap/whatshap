# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import logging
import itertools
from itertools import chain
import vcf
from cyvcf2 import VCF, Writer
from collections import defaultdict
import collections
from collections import deque
from operator import itemgetter
from .customcontainer import DefaultOrderedDict
import gc
import math
from sys import getsizeof, stderr
import numpy as np
try:
    from reprlib import repr
except ImportError:
    pass


def create_blocks(target_file, out_blockends_file, variant_set):
	vcf_reader=vcf.Reader(open(target_file, 'r'))
	i = 0
	blocks, ps = [], []
	intE, E, pairs = [], [], []
	count = 0
	phasedcount = 0
	haplo = {}
		 	
	for record in vcf_reader:
		for sample in vcf_reader.samples:
			#if record.genotype(sample).phased: #and record.POS in variant_list ? (von compute_ref vorher übergeben)
			if (record.POS in variant_set):		
				#phasing information for each variant position is stored in dict
				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
				phasedcount += 1
				#variants without a PS tag get an "unknown" value (this should be equal to the unphased ones)
			#	if not record.genotype(sample)['PS']:
				if not (record.genotype(sample).phased ):		
					count += 1
					ps.append('U'+str(count))
				else:
					phase_set = record.genotype(sample)['PS']
					ps.append(phase_set)
####  homozygous positions are currently left out #######
#			 if not (record.genotype(sample).phased ):
#				haplo[record.POS] = (record.genotype(sample)['GT'][0],record.genotype(sample)['GT'][2])
#				#homozygous sites are set to the phase set of the last variant position that was phased				
#				if (record.genotype(sample).gt_type == 2 and len(ps) != 0 and phasedcount >= 1):
#					index = -1
#					while not isinstance(ps[index], int):
#						index -= 1
#					if (index >= -len(ps)):
#						ps.append(ps[index])	
#				#unphased sites are set to "unknown"				
#				else:
#					count += 1
#					ps.append('U'+str(count))		
	ends = DefaultOrderedDict(list)
	starts = DefaultOrderedDict(list)
	#The end positions of each phase set are computed and stored as block ending positions
	for i in range(len(ps)):
		phase_set = ps[i]	
		if (i == 0):
			starts[ps[i]].append(i)
		if (i > 0 and phase_set != blocks[i-1]):
			ends[ps[i-1]].append(i-1)
			starts[ps[i]].append(i)
		if (i == len(ps)-1):
			ends[phase_set].append(i)
		blocks.append(phase_set)
	#pairs is a list that contains one list of tuples per phase set ("key")
	for key in ends.keys():
		if (len(ends[key]) == 1):
			pairs.append([(starts[key][0], ends[key][0])])
		else:
			pairlist = []
			for i in range(0,len(ends[key])):
				pairlist.append((starts[key][i], ends[key][i]))
			pairs.append(pairlist)
	#create list E containing block boundaries
	#add first block boundary (if part of an intermediate block, remove later)  
	E.append(pairs[0][-1][1])
	#internal block boundaries are found if present
	for i in range(1,len(pairs)):
		internal = False
		for j in range(0, i):
			if (pairs[i][-1][1] < pairs[j][-1][1]):
				internal = True
			#detect whether any items need to be removed since they are part of nested blocks
			if (len(pairs[i]) > 1):
				if (pairs[i][0][0] > pairs[j][0][0] and pairs[i][0][0] < pairs[j][-1][0] and pairs[i][-1][0] > pairs[j][-1][1]):
					if (pairs[j][-1][1] in E):
						E.remove(pairs[j][-1][1])
						intE.append((pairs[j]))	
		if internal:
			intE.append(pairs[i])
		else:
			E.append(pairs[i][-1][1])
	intE.sort()
	#add ending position of intermediate blocks (for nested blocks: of the last intermediate block) to the boundary set
	E2 = []
	for pairlist in intE:
		E2.append(pairlist[-1][1])
	E_whole = E+E2
	E_whole.sort()
	E_file = open(out_blockends_file,'w')
	for e in E_whole:
		if (e == E_whole[len(E_whole)-1]):
			E_file.write(str(e))
		else:
			E_file.write(str(e))
			E_file.write(',')
	return(E,haplo, intE)
	
#reads a VCF file and extracts the haplotype information, both haplotypes are written into a file as strings
def compute_haplotypes(newfile, resultfilename, variant_set):
	targetset = set()
	vcf_target = VCF(newfile)
	haplo1 = ""
	haplo2 = ""
	for v in vcf_target:
		if (v.POS in variant_set):
			targetset.add(v.POS)
			comp1 = 0
			comp2 = 0
			if (len(v.gt_bases[0].split('|')) > 1):
				leftbase = v.gt_bases[0].split('|')[0]
				rightbase = v.gt_bases[0].split('|')[1]
			elif(len(v.gt_bases[0].split('/')) > 1):
				leftbase = v.gt_bases[0].split('/')[0]
				rightbase = v.gt_bases[0].split('/')[1]
			if (leftbase == v.REF):
				comp1 = 0
			elif (leftbase == v.ALT[0]):
				comp1 = 1
			if (rightbase == v.REF):
				comp2 = 0
			elif (rightbase == v.ALT[0]):
				comp2 = 1
			haplo1 += str(comp1)	
			haplo2 += str(comp2)
	originalhaplofile = open(resultfilename,'w')
	originalhaplofile.write(haplo1)
	originalhaplofile.write('\n')	
	originalhaplofile.write(haplo2)
	return((haplo1,haplo2))

#tests the nested blocks
def update_haplotypes(haplofile, E, pathfile, intE):
	#read the haplotypes from file and store them in strings
	H_A, H_B = "",""
	with open(haplofile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				H_A = line.strip('\n')
			if (i==1):
				H_B = line.strip('\n')	
	#read the paths from file and store them in lists		
	helppath1 = []
	helppath2 = []
	with open(pathfile) as filestream:
		for i,line in enumerate(filestream):
			if (i==0):
				for num in line.strip().split(' '):				
					helppath1.append(num)
			elif (i==1):
				for num in line.strip().split(' '):
					helppath2.append(num)
	helppath1.append(helppath1[len(helppath1)-1])
	helppath2.append(helppath2[len(helppath2)-1])
	path1, path2 = [], []
	for i in range(0, len(helppath1)-1):
		path1.append(int(float(helppath1[i])))
	for i in range(0, len(helppath2)-1):
		path2.append(int(float(helppath2[i])))	
	
	#resolve any switches present in the paths and update the haplotypes accordingly
	(newhaplo1, newhaplo2) = improve_paths(list(path1),list(path2),E,intE, H_A,H_B)
	return((newhaplo1, newhaplo2))

#tests nested blocks
def improve_paths(path1, path2, E, intE, H_A, H_B):
	newpath1 = path1[:]
	newpath2 = path2[:]
	#store original paths for later use to not lose this information when updating the paths
	oldpath1 = newpath1[:]
	oldpath2 = newpath2[:]
	
	newhaplo1 = list(H_A)[:]
	newhaplo2 = list(H_B)[:]
	oldhaplo1 = newhaplo1[:]
	oldhaplo2 = newhaplo2[:]

	
	#at every block ending positions, it is checked whether the paths cross around that position
	for i in E:
		if (i != len(path1)-1 and i < len(path1)):
			#when a crossing between paths is found at position i, the paths are swapped from position 0 to i, the haplotypes accordingly
			if (path1[i] != path1[i+1] and path2[i] != path2[i+1] and (path1[i] == path2[i+1] or path2[i] == path1[i+1])):
				newpath1[0:i+1] = oldpath2[0:i+1]
				newpath2[0:i+1] = oldpath1[0:i+1]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				newhaplo1[0:i+1] = oldhaplo2[0:i+1]
				newhaplo2[0:i+1] = oldhaplo1[0:i+1]
				oldhaplo1 = newhaplo1[:]
				oldhaplo2 = newhaplo2[:]
	for pairlist in intE:
		switch = check_switching(pairlist,newpath1, newpath2)
		if switch:
			for pair in pairlist:				
				i = pair[0]
				j = pair[1]
				newpath1[i:j+1] = oldpath2[i:j+1]
				newpath2[i:j+1] = oldpath1[i:j+1]
				oldpath1 = newpath1[:]
				oldpath2 = newpath2[:]
				newhaplo1[i:j+1] = oldhaplo2[i:j+1]
				newhaplo2[i:j+1] = oldhaplo1[i:j+1]
				oldhaplo1 = newhaplo1[:]
				oldhaplo2 = newhaplo2[:]
	return("".join(newhaplo1),"".join(newhaplo2))

#tests nested blocks
def check_switching(pairlist, path1,path2):
	switch = 0
	to_switch = False
	print("path1: ", path1)
	for pair in pairlist:
		i = pair[0] 
		j = pair[1]
		if (path1[j] != path1[j+1] and path2[j] != path2[j+1] and (path1[j] == path2[j+1] or path2[j]==path1[j+1])):
			switch += 1
		if (path1[j] == path1[j+1] or path2[j]==path2[j+1]):
			switch -= 1
		if (i != 0):
			if (path1[i] != path1[i-1] and path2[i] != path2[i-1] and (path1[i] == path2[i-1] or path2[i]==path1[i-1])):
				switch += 1
			if (path1[i] == path1[i-1] or path2[i] == path2[i-1]):
				switch -= 1
	if (switch > 0):
		to_switch = True
	print("pairlist: ", pairlist, " to switch: ", to_switch)	
	return(to_switch)
	
def compute_referencepanel(ref_file, target_file, variant_set, ref_output):

	#compute size of ref_matrix and initialize it with zeros
	#width=number of variants between first and last relevant variant
	#height=twice the number of samples, as each sample offers two haplotypes
	
	vcf_ref = VCF(ref_file)
	height = 2*sum(1 for sample in vcf_ref.samples)
	width = 0
	z = 0
	start = 0
	end = 0


	#fill the matrix

	strlist = []
	string = ''
	for i in range(0,height):
		strlist.append(string)
	vcf_ref = VCF(ref_file, lazy=True)
	doublelist = []
	doubleset = set()
	for variant in vcf_ref:
		if (variant.POS in variant_set and variant.POS not in doubleset):
			doubleset.add(variant.POS)	
			counter = 0
			for base in variant.gt_bases:			
				if (base[0] == variant.REF):
					strlist[counter]+="0"
				else:		
					strlist[counter]+="1"
				if (base[(2)] == variant.REF):
					strlist[counter+1]+="0"
				else:
					strlist[counter+1]+="1"
				counter += 2
	output = open(ref_output, 'w')	
	for string in strlist:
		output.write(string+'\n')
	output.close()
	return (output)