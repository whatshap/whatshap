from libcpp cimport bool
cimport cpp

cdef class NumericSampleIds:
	cdef dict mapping


cdef class Read:
	cdef cpp.Read *thisptr
	cdef bool ownsptr


cdef class ReadSet:
	cdef cpp.ReadSet *thisptr


cdef class DPTable:
	cdef cpp.DPTable *thisptr


cdef class Pedigree:
	cdef cpp.Pedigree *thisptr
	cdef NumericSampleIds numeric_sample_ids


cdef class PedigreeDPTable:
	cdef cpp.PedigreeDPTable *thisptr
	cdef Pedigree pedigree
