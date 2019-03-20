import math
import logging
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint32_t
from cython.operator import dereference, postincrement
cimport cython
cimport cpp

cdef class HaploThreader:
	def __cinit__(self, ploidy, switchCost, affineSwitchCost):
		self.thisptr = new cpp.HaploThreader(ploidy, switchCost, affineSwitchCost)

	def computePaths(self, uint32_t num_vars, vector[vector[uint32_t]] covMap, vector[vector[double]] coverage, vector[vector[uint32_t]] consensus, vector[uint32_t] genotypes):
		cdef vector[vector[uint32_t]] path
		path = self.thisptr.computePaths(num_vars, covMap, coverage, consensus, genotypes)
		
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