from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int, int) except +
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

cdef class PyRead:
	cdef Read *thisptr
	cdef bool ownsptr

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

cdef class PyReadSet:
	cdef ReadSet *thisptr

cdef extern from "../src/columniterator.h":
	cdef cppclass ColumnIterator:
		ColumnIterator(ReadSet) except +

cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(ReadSet*, bool) except +
		void get_super_reads(ReadSet*) except +
		int get_optimal_score() except +
		vector[bool]* get_optimal_partitioning()

cdef class PyDPTable:
	cdef DPTable *thisptr

cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		bool contains(int) except +
		void add(int) except +
		int size() except +
		string toString() except +

cdef class PyIndexSet:
	cdef IndexSet *thisptr
