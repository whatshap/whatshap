import math
import logging
from collections import defaultdict
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp cimport bool
import numpy as np
import itertools as it
from cython.operator import dereference, postincrement
cimport cython
cimport cpp

cdef subsetting(num_vars, clustering, coverage, positions, cov_map, ploidy, genotypes, consensus):

	cdef int num_clusters = len(positions)
	cdef unordered_map[int, pair[int,int]] column
	cdef unordered_map[int, int] pred
	cdef vector[vector[pair[int,int]]] scoring
	cdef vector[pair[int,int]] newcolumn
	
	#initialize first column
	c_ids = cov_map[0]
	#compute the list of all <ploidy>-tuples that cover position 0
	c_tups = list(it.product(c_ids,repeat = ploidy))
	assert(len(c_tups)==len(c_ids)**ploidy)
	#restrict the set of tuples by those whose genotypes coincide. Tuples with the wrong genotype at position 0 are discarded
	geno_c_tups = [tup for tup in c_tups if (compute_tuple_genotype(consensus,tup, 0) == genotypes[0])]
#	geno_c_tups = [tup for tup in c_tups if (compute_tuple_genotype(consensus,tup, 0, genotypes[0]) <= 1)]
	if (len(geno_c_tups) == 0):
		geno_c_tups = c_tups
	#for the first column, only coverage costs are computed, the predecessor is set to -1
	for tup in geno_c_tups:
		newcolumn.push_back((cov_costs(tup, 0, coverage), -1))
	scoring.push_back(newcolumn)
	
	cdef int var
	cdef float mininum
	cdef int minimum_index
	cdef bool min_exists = False
	for var in range(1,num_vars):
		print("computing variant %s of %s variants in total " % (var, num_vars))
		#every column is a vector containing a pair of the costs at this position and the position of the predecessor (used for backtracing) 
		newcolumn.clear()
		#find the suitable clusters that cover position var and compute the list of cluster tuples
		c_ids = cov_map[var]
		c_tups = list(it.product(c_ids,repeat=ploidy))
		#create the genotype conform cluster tuples: compute the genotype expected by the tuple and the true genotype. If they are not equal, this tuple is omitted
		conform_tups = [tup for tup in c_tups if compute_tuple_genotype(consensus,tup, var) == genotypes[var]]
	#	conform_tups = [tup for tup in c_tups if compute_tuple_genotype(consensus,tup, var, genotypes[var]) <= 1]
		if (len(conform_tups) == 0):
			conform_tups = c_tups
		for tup in conform_tups:
			min_exists = False
			pred.clear()
			minimum = 1000000
			minimum_index = 0
			#for the previous column, compute the list of all possible tuples of clusters that appear at position var-1
			pred_tups = list(it.product(cov_map[var-1],repeat = ploidy))
			#compute the list of conform tuples, TODO: find an easier way
			conf_pred_tups = [predtup for predtup in pred_tups if (compute_tuple_genotype(consensus,predtup, var-1) == genotypes[var-1])]
	#		conf_pred_tups = [predtup for predtup in pred_tups if (compute_tuple_genotype(consensus,predtup, var-1, genotypes[var-1]) <=1)]
			if (len(conf_pred_tups) == 0):
				conf_pred_tups = pred_tups
			#compute the minimum of the previous column plus costs for switching:
			for pred_tup_i in range(len(conf_pred_tups)):	
				#the costs for the previous column are added to costs for switching from previous tuples to the current one
				pred[pred_tup_i] = (scoring[var-1][pred_tup_i].first+switch_costs(conf_pred_tups[pred_tup_i], tup,positions,var-1, ploidy))
				#find the  minimum in the previous column
				if pred[pred_tup_i] < minimum:
					min_exists = True
					minimum = pred[pred_tup_i]
					minimum_index = pred_tup_i
			#fill the matrix position with the computed costs and the index that was used from the column before (for simplifying backtracing)
			if min_exists:	
				newcolumn.push_back(((cov_costs(tup, var, coverage) + minimum), minimum_index))
			else:
				#no minimum exists: this is not expected to occur
				newcolumn.push_back((cov_costs(tup, var, coverage), 0))			
		scoring.push_back(newcolumn)	
	#convert cython into python data structure
	scoring_res = []	
	for i in range(num_vars): 
		newcolumn_res = []
		for j in range(len(scoring[i])):
			newcolumn_res.append(scoring[i][j])		
		scoring_res.append(newcolumn_res)
	return(scoring_res)

#computes the genotype that would belong to the given tuple tup at position var
def compute_tuple_genotype(consensus,tup, var):
	genotype = 0
	for i in tup:
		#allele = consensus[i][var]
		allele = consensus[var][i]
		genotype += allele
	return(genotype)
	
#
#def compute_tuple_genotype(consensus,tup, var, geno):
#	genotype = 0
#	for i in tup:
#		allele = consensus[i][var]
#		genotype += allele
#	res = max((geno-genotype),(genotype-geno))
#	return(res)

#computes costs for switching between two cluster tuples c_tuple1 and c_tuple2 at position var
cdef switch_costs(c_tuple1, c_tuple2, positions, var, ploidy):
	cdef int costs = 0
	#switch costs depend on the position: if var is the end of c_tuple1 or var+1 is the beginning of c_tuple2, switching is free
	cdef int i
	cdef int dissim = 0
	for i in range(0,ploidy):
		start = positions[c_tuple2[i]][0]
		end = positions[c_tuple1[i]][1]
	#	if (var != end and var+1 != start and (c_tuple1[i] != c_tuple2[i])):
		if (c_tuple1[i] != c_tuple2[i]):
			costs += 8
#			dissim += 1
#	if (dissim == 1):
#		costs == 32
#	elif (dissim == 2):
#		costs == 24
#	elif (dissim == 3):
#		costs == 16
#	elif (dissim == 4):
#		costs == 8

	#	if (var != end and (c_tuple1[i] != c_tuple2[i])):
	#		costs += 2
	#		dissim += 1
	#if (dissim > 1):
	#	costs += dissim
	#elif dissim == 1:
	#	costs += 5
	return(costs)

#computes the costs for differences between expected copy number (due to coverage) and the real copy number
#TODO: change 'hard' cutoffs to probability function
#TODO: does not work for a general <ploidy> yet
def cov_costs(c_tuple, var, coverage):
	costs = 0
	exp_cn = 0
	#compute copy numbers for every cluster in c_tuple
	for i in range(0,4):	
#	for i in range(0,2):
		cov = coverage[c_tuple[i]][var]
		#if cluster does not cover the position var:
		if (cov == 0):
			return (1000000)
		#else compare the expected copy number to the real one
		else:
			if (cov > 0 and cov < 0.125):
#			if (cov > 0 and cov < 0.2):
				exp_cn = 0
			if (cov >= 0.125 and cov < 0.375):
#			if (cov >= 0.2 and cov < 0.4):				
				exp_cn = 1
			if (cov >= 0.375 and cov < 0.625):
#			if (cov >= 0.4 and cov < 0.6):
				exp_cn = 2
			if (cov >= 0.625 and cov < 0.875):
#			if (cov >= 0.6 and cov < 0.8):
				exp_cn = 3
			if (cov >= 0.875 and cov <= 1):
#			if (cov >= 0.8 and cov <= 1):
				exp_cn = 4
#			if (cov > 0 and cov < 0.33):
#				exp_cn = 0
#			if (cov >= 0.33 and cov < 0.66):
#				exp_cn = 1
#			if (cov >= 0.66 and cov <= 1):
#				exp_cn = 2
		cn = c_tuple.count(c_tuple[i])
		if (exp_cn != cn):
			costs+= 1
	return(costs)

#def cov_costs_general(c_tuple, var, coverage):
#	costs = 0
#	exp_cns = [i for i in range(ploidy+1)]
#	for i in range(ploidy):
#		cov = coverage[c_tuple[i]][var]
#		if (cov == 0):
#			return(1000000)
#		else: 
#			for j in exp_cns:
#				if (cov > j/(ploidy+1) and cov < (j+1)/(ploidy+1)):
#					exp_cn = j
#		cn = c_tuple.count(c_tuple[i])
#		if (cn != exp_cn):
#			costs +=1
#		return(costs)
			

def compute_index(tup, ploidy, num_clusters):
	index = 0
	for i in range(ploidy-1,-1,-1):
		index += tup[(ploidy-1)-i]*(num_clusters**i)
	return(index)
	
def clustering_DP(num_vars,clustering,coverage,positions, cov_map, ploidy, genotypes, consensus):
	scoring_matrix = subsetting(num_vars, clustering, coverage,positions, cov_map, ploidy, genotypes, consensus)
	return(scoring_matrix)

