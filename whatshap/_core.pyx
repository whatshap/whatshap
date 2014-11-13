# distutils: language = c++
# distutils: sources = src/activelist.cpp src/activelistdelegator.cpp src/columncostcomputer.cpp src/columnindexingiterator.cpp src/columnindexingscheme.cpp src/dp.cpp src/dptable.cpp src/entry.cpp src/graycodes.cpp src/read.cpp src/readset.cpp src/columniterator.cpp

from libcpp.string cimport string

# ====== Read ======
cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int) except +
		Read(Read) except +
		string toString()
		void addVariant(int, char, int, int)

cdef class PyRead:
	cdef Read *thisptr
	def __cinit__(self, str name, int mapq):
		# TODO: Is this the best way to handle string arguments?
		cdef string _name = name.encode('UTF-8')
		self.thisptr = new Read(_name, mapq)
	def __dealloc__(self):
		del self.thisptr
	def __str__(self):
		return str(self.thisptr.toString())
	def addVariant(self, int position, str base, int allele, int quality):
		assert len(base) == 1
		self.thisptr.addVariant(position, ord(base[0]), allele, quality)
	# TODO: add getters

# ====== ReadSet ======
cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*)
		string toString()
		void finalize()

cdef class PyReadSet:
	cdef ReadSet *thisptr
	def __cinit__(self):
		self.thisptr = new ReadSet()
	def __dealloc__(self):
		del self.thisptr
	def add(self, PyRead read):
		self.thisptr.add(new Read(read.thisptr[0]))
	def __str__(self):
		return self.thisptr.toString().decode('utf-8')
	def finalize(self):
		self.thisptr.finalize()
	# TODO: add getters

# ====== ColumnIterator ======
cdef extern from "../src/columniterator.h":
	cdef cppclass ColumnIterator:
		ColumnIterator(ReadSet) except +

# ====== DPTable ======
cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(bool) except +
		void compute_table(ColumnIterator*);

cdef class PyDPTable:
	cdef DPTable *thisptr
	def __cinit__(self, all_heterozygous, PyReadSet readset):
		self.thisptr = new DPTable(all_heterozygous)
		cdef ColumnIterator* iterator = new ColumnIterator(readset.thisptr[0])
		print('CALLING compute_table BEGIN')
		self.thisptr.compute_table(iterator)
		print('CALLING compute_table END')
		del iterator
	def __dealloc__(self):
		del self.thisptr
