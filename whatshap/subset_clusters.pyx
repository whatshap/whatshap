import math
import logging
from collections import defaultdict
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool
import numpy as np
cimport cython
cimport cpp

cdef subsetting(num_vars, cluster_tuples, coverage, positions):
	print("in subset_clustering")
#	cdef float** Score = new float*[num_vars]
#	scoring = np.zeros([len(cluster_tuples), num_vars],dtype='f,f')
	cdef vector[vector[pair[int,int]]] scoring
	cdef vector[pair[int,int]] column
	cdef vector[int] pred
	uncovered = defaultdict(set)
	print("positions: ", positions)
	print("positions keys: ", positions.keys())
	#initialize first column
	for c_tuple in range(len(cluster_tuples)):
		column.push_back((cov_costs(cluster_tuples[c_tuple], 0, coverage), -1))
	scoring.push_back(column)
	cdef int var
	cdef float mininum
	cdef int minimum_index
	cdef vector[pair[pair[int,int],int]] test
	for var in range(1,num_vars):
		print("var in subset clusters: ", var)
		column.clear()				
		for c_tuple in xrange(len(cluster_tuples)):	
			if (cov_costs(cluster_tuples[c_tuple], var, coverage) == 1000000):
				column.push_back((1000000, c_tuple))
				uncovered[var].add(c_tuple)
			else:
				pred.clear()
			#	minimum = float('inf')
				minimum = 1000000		
				minimum_index = 0
				for i in range(len(cluster_tuples)):
					if (i in uncovered[var-1]):
						pred.push_back(1000000)
					else:
						pred.push_back(scoring[var-1][i].first+switch_costs(cluster_tuples[c_tuple], cluster_tuples[i],positions,var-1))
					if pred[i] < minimum:
						minimum = pred[i]
						minimum_index = i
						#pred.push_back(scoring[var-1][i].first + cov_costs(cluster_tuples[c_tuple], var, coverage) )
			#	minimum, minimum_index = min((val, idx) for (idx,val) in enumerate(pred))
				#if (var == 1766):
				#	print("pred: ", pred)
				column.push_back(((cov_costs(cluster_tuples[c_tuple], var, coverage) + minimum), minimum_index))
				if (var == 2):
					test.push_back(((cov_costs(cluster_tuples[c_tuple], var, coverage),minimum), c_tuple))
		scoring.push_back(column)
#	print("scoring last column: ", scoring[1770])
	print("in subset clustering: length scoring", scoring.size())
	print("in subset clustering: length scoring[0]", scoring[0].size())
#	print("last col: ", scoring[1770])
#	print("not 1000: ", [i for i in scoring[1770] if i[0] != 1000000])
	print("test: ", test)
	scoring_res = [[0]*len(cluster_tuples) for i in range(num_vars)]	
	for i in range(num_vars):
		for j in range(len(cluster_tuples)):
			scoring_res[i][j] = scoring[i][j]

	return(scoring_res)

cdef switch_costs(c_tuple1, c_tuple2, positions, var):
	cdef int costs = 0
	cdef vector[int] starts
	cdef vector[int] ends
	#switch costs depend on the position: if var is the end of c_tuple1 or var+1 is the beginning of c_tuple2, switching is free
	cdef int i
	cdef list starts_list = [j[0] for j in positions[c_tuple2[0]]]
	cdef list ends_list = [j[1] for j in positions[c_tuple1[0]]]
	for i in range(0,4):
		starts_list = [j[0] for j in positions[c_tuple2[i]]]
		ends_list = [j[1] for j in positions[c_tuple1[i]]]
	#	starts.clear()
	#	ends.clear()
	#	for j in positions[c_tuple2[i]]:
	#		starts.push_back(j[0])
	#	for j in positions[c_tuple1[i]]:
	#		ends.push_back(j[1])
		if (var not in ends_list and var+1 not in starts_list):
	#	if (find(ends.begin(), ends.end(), var) != ends.end() and find(starts.begin(), starts.end(), var+1) != starts.end()):	
			if (c_tuple1[i] != c_tuple2[i]):
				costs += 1
	return(costs)

def cov_costs(c_tuple, var, coverage):
	costs = 0
	exp_cn = 0
	#compute copy numbers for every cluster in c_tuple
	for i in range(0,4):
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
		cn = c_tuple.count(c_tuple[i])
		if (exp_cn != cn):
			costs+= 1
	return(costs)

def clustering_DP(num_vars,cluster_tuples,coverage,positions):
	scoring_res = subsetting(num_vars, cluster_tuples, coverage,positions)
	return(scoring_res)

