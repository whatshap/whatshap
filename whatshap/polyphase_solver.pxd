# cython: language_level=3

from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp

cdef class ClusterEditingSolver:
	cdef cpp.ClusterEditingSolver *thisptr
	cdef TriangleSparseMatrix m


cdef class TriangleSparseMatrix:
	cdef cpp.TriangleSparseMatrix *thisptr


cdef class ReadScoring:
	cdef cpp.ReadScoring *thisptr


cdef class HaploThreader:
	cdef cpp.HaploThreader *thisptr


cdef class SwitchFlipCalculator:
	cdef cpp.SwitchFlipCalculator *thisptr
	cdef uint32_t ploidy
