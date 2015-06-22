# Do not use the distutils directives here, but configure everything in
# setup.py, such as language and sources. It would work during development, but
# during a regular installation, the module will be compiled from the
# pre-generated .c file and the .pyx file is not read.

from collections import namedtuple

# A single variant on a read.
Variant = namedtuple('Variant', 'position allele quality')

# ====== Read ======
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

	def __cinit__(self, str name = None, int mapq = 0, int source_id = 0):
		cdef string _name = ''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			self.thisptr = new Read(_name, mapq, source_id)
			self.ownsptr = True

	def __dealloc__(self):
		if self.ownsptr:
			assert self.thisptr != NULL
			del self.thisptr

	def __str__(self):
		assert self.thisptr != NULL
		return self.thisptr.toString().decode('utf-8')

	property mapqs:
		def __get__(self):
			assert self.thisptr != NULL
			return tuple(self.thisptr.getMapqs())

	property name:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getName().decode('utf-8')

	property source_id:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getSourceID()

	def __iter__(self):
		"""Iterate over all variants in this read"""
		assert self.thisptr != NULL
		for i in range(len(self)): 
			yield self[i]

	def __len__(self):
		"""Return number of variants in this read"""
		assert self.thisptr != NULL
		return self.thisptr.getVariantCount()

	def __getitem__(self, key):
		"""Return Variant object at the given integer index"""
		assert self.thisptr != NULL
		if isinstance(key, slice):
			raise NotImplementedError("Read does not support slices")
		assert isinstance(key, int)
		cdef int n = self.thisptr.getVariantCount()
		if not (-n <= key < n):
			raise IndexError('Index out of bounds: {}'.format(key))
		if key < 0:
			key = n + key
		return Variant(
			position=self.thisptr.getPosition(key),
			allele=self.thisptr.getAllele(key),
			quality=self.thisptr.getVariantQuality(key)
		)

	def __setitem__(self, index, variant):
		assert self.thisptr != NULL
		cdef int n = self.thisptr.getVariantCount()
		if not (-n <= index < n):
			raise IndexError('Index out of bounds: {}'.format(index))
		if index < 0:
			index = n + index
		if not isinstance(variant, Variant):
			raise ValueError('Expected instance of Variant, but found {}'.format(type(variant)))
		self.thisptr.setPosition(index, variant.position)
		self.thisptr.setAllele(index, variant.allele)
		self.thisptr.setVariantQuality(index, variant.quality)

	def __contains__(self, position):
		"""Return whether this read contains a variant at the given position.
		A linear search is used.
		"""
		assert self.thisptr != NULL
		assert isinstance(position, int)
		for variant in self:
			if variant.position == position:
				return True
		return False

	def add_variant(self, int position, int allele, int quality):
		assert self.thisptr != NULL
		self.thisptr.addVariant(position, allele, quality)

	def add_mapq(self, int mapq):
		assert self.thisptr != NULL
		self.thisptr.addMapq(mapq)

	def sort(self):
		assert self.thisptr != NULL
		self.thisptr.sortVariants()

	property is_sorted:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.isSorted()


# ====== IndexSet ======
cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		bool contains(int) except +
		void add(int) except +
		int size() except +
		string toString() except +


# ====== ReadSet ======
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
		cdef string name = ''
		cdef Read* cread = NULL
		cdef PyRead read = PyRead()
		if isinstance(key, int):
			read.thisptr = self.thisptr.get(key)
		elif isinstance(key, str):
			raise NotImplementedError('Querying a ReadSet by read name is deprecated, please query by (source_id, name) instead')
		elif isinstance(key, tuple) and (len(key) == 2) and (isinstance(key[0],int) and isinstance(key[1],str)):
			source_id = key[0]
			name = key[1].encode('UTF-8')
			cread = self.thisptr.getByName(name, source_id)
			if cread == NULL:
				raise KeyError(key)
			else:
				read.thisptr = cread
		else:
			assert False, 'Invalid key: {}'.format(key)
		return read

	#def get_by_name(self, name):
		#cdef string _name = name.encode('UTF-8')
		#cdef Read* cread = self.thisptr.getByName(_name)
		#cdef PyRead read = PyRead()
		#if cread == NULL:
			#raise KeyError(name)
		#else:
			#read.thisptr = cread
			#return read

	def sort(self):
		"""Sort contained reads by the position of the first variant they contain. Note that 
		this is not necessarily the variant with the lowest position, unless sort() has been 
		called on all contained reads. Ties are resolved by comparing the read name."""
		self.thisptr.sort()

	def subset(self, reads_to_select):
		# TODO: is there a way of avoiding to unecessarily creating/destroying a ReadSet object?
		cdef IndexSet* index_set = new IndexSet()
		cdef int i
		for i in reads_to_select:
			index_set.add(i)
		result = PyReadSet()
		del result.thisptr
		result.thisptr = self.thisptr.subset(index_set)
		del index_set
		return result

	def get_positions(self):
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
		void get_super_reads(ReadSet*) except +
		int get_optimal_score() except +
		vector[bool]* get_optimal_partitioning()


cdef class PyDPTable:
	cdef DPTable *thisptr

	def __cinit__(self, PyReadSet readset, all_heterozygous):
		"""Build the DP table from the given read set which is assumed to be sorted;
		that is, the variants in each read must be sorted by position and the reads
		in the read set must also be sorted (by position of their left-most variant).
		"""
		self.thisptr = new DPTable(readset.thisptr, all_heterozygous)

	def __dealloc__(self):
		del self.thisptr

	def get_super_reads(self):
		"""Obtain optimal-score haplotypes.
		IMPORTANT: The ReadSet given at construction time must not have been altered.
		DPTable retained a pointer to this set and will access it again. If it has
		been altered, behavior is undefined.
		TODO: Change that.
		"""
		result = PyReadSet()
		self.thisptr.get_super_reads(result.thisptr)
		return result

	def get_optimal_cost(self):
		"""Returns the cost resulting from solving the Minimum Error Correction (MEC) problem."""
		return self.thisptr.get_optimal_score()

	def get_optimal_partitioning(self):
		"""Returns a list of the same size as the read set, where each entry is either 0 or 1,
		telling whether the corresponding read is in partition 0 or in partition 1,"""
		cdef vector[bool]* p = self.thisptr.get_optimal_partitioning()
		result = [0 if x else 1 for x in p[0]]
		del p
		return result

