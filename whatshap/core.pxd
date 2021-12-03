# cython: language_level=3

from libcpp cimport bool
from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp


cdef class NumericSampleIds:
	cdef dict mapping
	cdef bool frozen


cdef class Read:
	cdef cpp.Read *thisptr
	cdef bool ownsptr


cdef class ReadSet:
	cdef cpp.ReadSet *thisptr


cdef class Pedigree:
	cdef cpp.Pedigree *thisptr
	cdef NumericSampleIds numeric_sample_ids


# cdef class PedigreeDPTable:
# 	cdef cpp.PedigreeDPTable *thisptr
# 	cdef Pedigree pedigree


cdef class PhredGenotypeLikelihoods:
	cdef cpp.PhredGenotypeLikelihoods *thisptr
	
	
cdef class Genotype:
	cdef cpp.Genotype *thisptr
	cdef uint64_t index
	cdef uint32_t ploidy


cdef class GenotypeHMM:
	cdef cpp.GenotypeHMM *thisptr
	cdef Pedigree pedigree
	cdef NumericSampleIds numeric_sample_ids

# cdef class HapChatCore:
# 	cdef cpp.HapChatCore *thisptr


# cdef class ClusterEditingSolver:
# 	cdef cpp.ClusterEditingSolver *thisptr
# 	cdef TriangleSparseMatrix m


# cdef class TriangleSparseMatrix:
# 	cdef cpp.TriangleSparseMatrix *thisptr


# cdef class ReadScoring:
# 	cdef cpp.ReadScoring *thisptr


# cdef class HaploThreader:
# 	cdef cpp.HaploThreader *thisptr
	
# cdef class SwitchFlipCalculator:
# 	cdef cpp.SwitchFlipCalculator *thisptr
# 	cdef uint32_t ploidy
