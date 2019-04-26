import math
import logging
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint32_t
from cython.operator import dereference, postincrement
cimport cython
cimport cpp

cdef class HaploThreader:
	def __cinit__(self, ploidy, switchCost, affineSwitchCost, symmetryOptimization, rowLimit):
		self.thisptr = new cpp.HaploThreader(ploidy, switchCost, affineSwitchCost, symmetryOptimization, rowLimit)
		
	def computePathsBlockwise(self, vector[uint32_t]& blockStarts, vector[vector[uint32_t]]& covMap, vector[vector[double]]& coverage, vector[vector[uint32_t]]& consensus, vector[unordered_map[uint32_t, uint32_t]]& genotypes, vector[vector[vector[double]]]& clusterDissim):
		cdef vector[vector[uint32_t]] path
		path = self.thisptr.computePaths(blockStarts, covMap, coverage, consensus, genotypes, clusterDissim)
		
		# convert to python data structure
		py_path = []
		path_len = path.size()
		if path_len > 0:
			ploidy = path[0].size()
			for i in range(path_len):
				pos = []
				for j in range(ploidy):
					pos.append(path[i][j])
				py_path.append(pos)
		
		return py_path

	def computePaths(self, uint32_t start, uint32_t end, vector[vector[uint32_t]]& covMap, vector[vector[double]]& coverage, vector[vector[uint32_t]]& consensus, vector[unordered_map[uint32_t, uint32_t]]& genotypes, vector[vector[vector[double]]]& clusterDissim):
		cdef vector[vector[uint32_t]] path
		path = self.thisptr.computePaths(start, end, covMap, coverage, consensus, genotypes, clusterDissim)
		
		# convert to python data structure
		py_path = []
		path_len = path.size()
		if path_len > 0:
			ploidy = path[0].size()
			for i in range(path_len):
				pos = []
				for j in range(ploidy):
					pos.append(path[i][j])
				py_path.append(pos)
		
		return py_path