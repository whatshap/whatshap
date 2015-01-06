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
		string toString() except +
		void addVariant(int, char, int, int) except +
		string getName() except +
		vector[int] getMapqs() except +
		void addMapq(int) except +
		int getPosition(int) except +
		char getBase(int) except +
		int getAllele(int) except +
		int getBaseQuality(int) except +
		int getVariantCount() except +
		void sortVariants() except +

cdef class PyRead:
	cdef Read *thisptr
	cdef bool ownsptr

	def __cinit__(self, str name = None, int mapq = 0):
		cdef string _name = ''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			self.thisptr = new Read(_name, mapq)
			self.ownsptr = True

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

	def addVariant(self, int position, str base, int allele, int quality):
		assert self.thisptr != NULL
		assert len(base) == 1
		self.thisptr.addVariant(position, ord(base[0]), allele, quality)

	def addMapq(self, int mapq):
		assert self.thisptr != NULL
		self.thisptr.addMapq(mapq)

	def sort(self):
		assert self.thisptr != NULL
		self.thisptr.sortVariants()


# ====== ReadSet ======
cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*) except +
		string toString() except +
		int size() except +
		void sort() except +
		Read* get(int) except +
		Read* getByName(string) except +
		ReadSet* subset(IndexSet*) except +
		# TODO: Check why adding "except +" here doesn't compile
		vector[unsigned int]* get_positions()


cdef class PyReadSet:
	cdef ReadSet *thisptr

	def __cinit__(self):
		self.thisptr = new ReadSet()

	def __dealloc__(self):
		del self.thisptr

	def add(self, PyRead read):
		"""Adds a read to the set. 
		WARNING: this will internally create a copy of the wrapped C++ Read object,
		so that subsequent changes to the Read don't affect the 
		newly created copy that is added to the ReadSet."""
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
		cdef PyRead read = PyRead()
		read.thisptr = self.thisptr.get(key)
		return read

	def getByName(self, name):
		cdef string _name = name.encode('UTF-8')
		cdef Read* cread = self.thisptr.getByName(_name)
		cdef PyRead read = PyRead()
		if cread == NULL:
			raise KeyError(name)
		else:
			read.thisptr = cread
			return read

	def sort(self):
		"""Sort contained reads by the position of the first variant they contain. Note that 
		this is not necessarily the variant with the lowest position, unless sort() has been 
		called on all contained reads. Ties are resolved by comparing the read name."""
		self.thisptr.sort()

	def subset(self, PyIndexSet index_set):
		# TODO: is there a way of avoiding to unecessarily creating/destroying a ReadSet object?
		result = PyReadSet()
		del result.thisptr
		result.thisptr = self.thisptr.subset(index_set.thisptr)
		return result

	def getPositions(self):
		cdef vector[unsigned int]* v = self.thisptr.get_positions()
		result = list(v[0])
		del v
		return result


# ====== ColumnIterator ======
cdef extern from "../src/columniterator.h":
	cdef cppclass ColumnIterator:
		ColumnIterator(ReadSet) except +


# ====== DPTable ======
cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(ReadSet*, bool) except +
		void compute_table() except +
		void get_super_reads(ReadSet*) except +


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
		bool contains(int) except +
		void add(int) except +
		int size() except +
		string toString() except +


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
