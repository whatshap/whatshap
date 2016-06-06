# kate: syntax python
"""
Declarations for all C++ classes that are wrapped from Cython.
"""
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int, int, int) except +
		Read(Read) except +
		string toString() except +
		void addVariant(int, int, int) except +
		string getName() except +
		vector[int] getMapqs() except +
		void addMapq(int) except +
		int getPosition(int) except +
		void setPosition(int, int)  except +
		int getAllele(int) except +
		void setAllele(int, int) except +
		int getVariantQuality(int) except +
		void setVariantQuality(int, int) except +
		int getVariantCount() except +
		void sortVariants() except +
		bool isSorted() except +
		int getSourceID() except +
		int getSampleID() except +


cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		bool contains(int) except +
		void add(int) except +
		int size() except +
		string toString() except +


cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*) except +
		string toString() except +
		int size() except +
		void sort() except +
		Read* get(int) except +
		Read* getByName(string, int) except +
		ReadSet* subset(IndexSet*) except +
		# TODO: Check why adding "except +" here doesn't compile
		vector[unsigned int]* get_positions()


cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(ReadSet*, bool) except +
		void get_super_reads(ReadSet*) except +
		int get_optimal_score() except +
		vector[bool]* get_optimal_partitioning()


cdef extern from "../src/pedigree.h":
	cdef cppclass Pedigree:
		Pedigree() except +
		void addIndividual(unsigned int id, vector[unsigned int] genotypes) except +
		void addRelationship(unsigned int m, unsigned int f, unsigned int c) except +
		unsigned int size()
		string toString() except +


cdef extern from "../src/pedigreedptable.h":
	cdef cppclass PedigreeDPTable:
		PedigreeDPTable(ReadSet*, vector[unsigned int], Pedigree* pedigree) except +
		void get_super_reads(vector[ReadSet*]*, vector[unsigned int]* transmission_vector) except +
		int get_optimal_score() except +
		vector[bool]* get_optimal_partitioning()
