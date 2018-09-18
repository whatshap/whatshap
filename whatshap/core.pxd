from libcpp cimport bool
cimport cpp

cdef class NumericSampleIds:
	cdef dict mapping


cdef class Read:
	cdef cpp.Read *thisptr
	cdef bool ownsptr


cdef class ReadSet:
	cdef cpp.ReadSet *thisptr


cdef class Pedigree:
	cdef cpp.Pedigree *thisptr
	cdef NumericSampleIds numeric_sample_ids


cdef class PedigreeDPTable:
	cdef cpp.PedigreeDPTable *thisptr
	cdef Pedigree pedigree


cdef class GenotypeLikelihoods:
	cdef cpp.GenotypeLikelihoods *thisptr

cdef class Genotype:
	cdef cpp.Genotype *thisptr

cdef class GenotypeDPTable:
	cdef cpp.GenotypeDPTable *thisptr
	cdef Pedigree pedigree
	cdef NumericSampleIds numeric_sample_ids
