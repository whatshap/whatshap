# Do not use the distutils directives here, but configure everything in
# setup.py, such as language and sources. It would work during development, but
# during a regular installation, the module will be compiled from the
# pre-generated .c file and the .pyx file is not read.

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

# ====== Read ======
cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int) except +
		Read(Read) except +
		string toString()
		void addVariant(int, char, int, int)
		string getName()
		vector[int] getMapqs()
		void addMapq(int)
		int getPosition(int)
		char getBase(int)
		int getAllele(int)
		int getBaseQuality(int)
		int getVariantCount()
		void sortVariants()


cdef class PyFrozenRead:
	cdef Read *thisptr
	cdef bool ownsptr

	def __cinit__(self):
		self.thisptr = NULL
		self.ownsptr = False

	def __dealloc__(self):
		if self.ownsptr:
			assert self.thisptr != NULL
			del self.thisptr

	def __str__(self):
		assert self.thisptr != NULL
		return self.thisptr.toString().decode('utf-8')

	def getMapqs(self):
		assert self.thisptr != NULL
		return tuple(self.thisptr.getMapqs())

	def getName(self):
		assert self.thisptr != NULL
		return self.thisptr.getName().decode('utf-8')

	def __iter__(self):
		assert self.thisptr != NULL
		for i in range(len(self)): 
			yield self[i]

	def __len__(self):
		assert self.thisptr != NULL
		return self.thisptr.getVariantCount()

	def __getitem__(self, key):
		assert self.thisptr != NULL
		if isinstance(key, slice):
			raise NotImplementedError("Read does not support slices")
		assert isinstance(key, int)
		if not (0 <= key < self.thisptr.getVariantCount()):
			raise IndexError('Index out of bounds: {}'.format(key))
		return (self.thisptr.getPosition(key), chr(self.thisptr.getBase(key)), self.thisptr.getAllele(key), self.thisptr.getBaseQuality(key))


cdef class PyRead(PyFrozenRead):
	def __cinit__(self, str name, int mapq):
		# TODO: Is this the best way to handle string arguments?
		cdef string _name = name.encode('UTF-8')
		self.thisptr = new Read(_name, mapq)
		self.ownsptr = True

	def addVariant(self, int position, str base, int allele, int quality):
		assert len(base) == 1
		self.thisptr.addVariant(position, ord(base[0]), allele, quality)

	def addMapq(self, int mapq):
		self.thisptr.addMapq(mapq)

	def sort(self):
		self.thisptr.sortVariants()


# ====== ReadSet ======
cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*)
		string toString()
		int size()
		void finalize()
		bool isFinalized()
		Read* get(int)
		Read* getByName(string)
		ReadSet* subset(IndexSet*)
		vector[unsigned int]* get_positions()


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

	def __iter__(self):
		for i in range(self.thisptr.size()):
			yield self[i]

	def __len__(self):
		return self.thisptr.size()

	def __getitem__(self, key):
		if isinstance(key, slice):
			raise NotImplementedError('ReadSet does not support slices')
		assert isinstance(key, int)
		cdef PyFrozenRead read = None
		if self.thisptr.isFinalized():
			read = PyFrozenRead()
			read.thisptr = self.thisptr.get(key)
		else:
			read = PyRead('', 0)
			read.thisptr = self.thisptr.get(key)
			read.ownsptr = False
		return read

	def getByName(self, name):
		cdef string _name = name.encode('UTF-8')
		cdef Read* cread = self.thisptr.getByName(_name)
		cdef PyFrozenRead read = None
		if cread == NULL:
			raise KeyError(name)
		else:
			if self.thisptr.isFinalized():
				read = PyFrozenRead()
			else:
				read = PyRead('', 0)
			read.thisptr = cread
			read.ownsptr = False
			return read

	def finalize(self):
		self.thisptr.finalize()

	def subset(self, PyIndexSet index_set):
		# TODO: is there a way of avoiding to unecessarily creating/destroying a ReadSet object?
		result = PyReadSet()
		del result.thisptr
		result.thisptr = self.thisptr.subset(index_set.thisptr)
		return result

	def getPositions(self):
		return self.thisptr.get_positions()[0]


# ====== ColumnIterator ======
cdef extern from "../src/columniterator.h":
	cdef cppclass ColumnIterator:
		ColumnIterator(ReadSet) except +


# ====== DPTable ======
cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(ReadSet*, bool) except +
		void compute_table()
		void get_super_reads(ReadSet*)


cdef class PyDPTable:
	cdef DPTable *thisptr

	def __cinit__(self, PyReadSet readset, all_heterozygous):
		self.thisptr = new DPTable(readset.thisptr, all_heterozygous)

	def __dealloc__(self):
		del self.thisptr

	def getSuperReads(self):
		result = PyReadSet()
		self.thisptr.get_super_reads(result.thisptr)
		return result


# ====== IndexSet ======
cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		bool contains(int)
		void add(int)
		int size()
		string toString()


cdef class PyIndexSet:
	cdef IndexSet *thisptr

	def __cinit__(self):
		self.thisptr = new IndexSet()

	def __dealloc__(self):
		del self.thisptr

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')

	def __len__(self):
		return self.thisptr.size()

	def contains(self, index):
		return self.thisptr.contains(index)

	def add(self, index):
		self.thisptr.add(index)
