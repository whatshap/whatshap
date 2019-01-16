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

cdef subsetting(num_vars, all_tuples, coverage, positions, cov_map, ploidy):
	cdef int num_clusters = len(positions)
	cdef unordered_map[int, pair[int,int]] column
	cdef unordered_map[int, int] pred
	cdef vector[unordered_map[int, pair[int,int]]] scoring

	#initialize first column
	c_ids = cov_map[0]
	#compute the list of all <ploidy>-tuples that cover position 0
	c_tups = list(it.product(c_ids,repeat = ploidy))
	assert(len(c_tups)==len(c_ids)**ploidy)
	for tup in c_tups:
		#compute the index of that specific cluster in the list of all clusters
		tup_ind = compute_index(tup, ploidy, num_clusters)
		#in the first column, only costs for coverage (none for switching) must be computed
		column[tup_ind] = (cov_costs(tup, 0, coverage), -1)
	scoring.push_back(column)
	
	cdef int var
	cdef float mininum
	cdef int minimum_index
	cdef bool min_exists = False
	cdef int counter = 0
	for var in range(1,num_vars):
		print("computing variant %s of %s variants in total " % (var, num_vars))
		#every column is a map, mapping the index of the cluster tuple (with respect to all existing tuples) 
		#to a pair of the costs at this position and the predecessor (used for backtracing) 
		column.clear()
		#find the suitable clusters that cover position var and compute the list of cluster tuples
		c_ids = cov_map[var]
		c_tups = list(it.product(c_ids,repeat=ploidy))
		for tup in c_tups:
			#compute_index returns the index of this <ploidy>-tuple in the list of all_tuples. This index is consistently used as key in the maps representing the columns
			tup_ind = compute_index(tup, ploidy, num_clusters)
			min_exists = False
			pred.clear()
			minimum = 1000000
			minimum_index = 0
	
			#for the previous column, compute the list of all possible tuples of clusters that appear at position var-1
			pred_tups = list(it.product(cov_map[var-1],repeat = ploidy))
			#compute the minimum of the column before plus costs for switching:
			for pred_tup in pred_tups:			
				#compute_index returns the index of this <ploidy>-tuple in the list of all_tuples
				pred_tup_ind = compute_index(pred_tup, ploidy, num_clusters)
				#the costs for the previous column are added to costs for switching from previous tuples to the current one
				pred[pred_tup_ind] = (scoring[var-1][pred_tup_ind].first+switch_costs(pred_tup, tup,positions,var-1, ploidy))
				#find the  minimum in the previous column
				if pred[pred_tup_ind] < minimum:
					min_exists = True
					minimum = pred[pred_tup_ind]
					minimum_index = pred_tup_ind
				#fill the matrix position with the computed costs and the index that was used from the column before (for simplifying backtracing)
			if min_exists:				
				column[tup_ind]= ((cov_costs(tup, var, coverage) + minimum), minimum_index)
			else:
				#no minimum exists: this is not expected to occur
				column[tup_ind] = ((cov_costs(tup, var, coverage), 0))
				counter+= 0
		scoring.push_back(column)
	
	#convert cython into python data structure
	scoring_res = []	
	for i in range(num_vars):
		column_res = {} 
		for j in scoring[i]:
			j_key = j.first
			column_res[j_key] = scoring[i][j_key]
		scoring_res.append(column_res)

	return(scoring_res)

#computes costs for switching between two cluster tuples c_tuple1 and c_tuple2 at position var
cdef switch_costs(c_tuple1, c_tuple2, positions, var, ploidy):
	cdef int costs = 0
	cdef vector[int] starts
	cdef vector[int] ends
	#switch costs depend on the position: if var is the end of c_tuple1 or var+1 is the beginning of c_tuple2, switching is free
	cdef int i
	cdef list starts_list = [j[0] for j in positions[c_tuple2[0]]]
	cdef list ends_list = [j[1] for j in positions[c_tuple1[0]]]
	for i in range(0,ploidy):
		starts_list = [j[0] for j in positions[c_tuple2[i]]]
		ends_list = [j[1] for j in positions[c_tuple1[i]]]
		if (var not in ends_list and var+1 not in starts_list):
			if (c_tuple1[i] != c_tuple2[i]):
				costs += 1
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
				exp_cn = 0
			if (cov >= 0.125 and cov < 0.375):
				exp_cn = 1
			if (cov >= 0.375 and cov < 0.625):
				exp_cn = 2
			if (cov >= 0.625 and cov < 0.875):
				exp_cn = 3
			if (cov >= 0.875 and cov <= 1):
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

def compute_index(tup, ploidy, num_clusters):
	index = 0
	for i in range(ploidy-1,-1,-1):
		index += tup[(ploidy-1)-i]*(num_clusters**i)
	return(index)
	
def clustering_DP(num_vars,all_tuples,coverage,positions, cov_map, ploidy):
	scoring_matrix = subsetting(num_vars, all_tuples, coverage,positions, cov_map, ploidy)
	return(scoring_matrix)

