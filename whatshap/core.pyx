# cython: language_level=3

"""
Wrappers for core C++ classes.
"""
# Do not use the distutils directives here, but configure everything in
# setup.py, such as language and sources. It would work during development, but
# during a regular installation, the module will be compiled from the
# pre-generated .cpp file and the .pyx file is not read.

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp

from .variant import Variant
from collections import namedtuple
from cython.operator cimport dereference as deref
from libcpp.deque cimport deque
from libcpp.pair cimport pair


cdef class NumericSampleIds:
	"""
	Mapping of sample names (strings) to numeric ids.
	"""
	def __cinit__(self):
		self.mapping = {}
		self.frozen = False

	def __getitem__(self, sample):
		if not self.frozen and sample not in self.mapping:
			self.mapping[sample] = len(self.mapping)
		return self.mapping[sample]

	def __len__(self):
		return len(self.mapping)

	def __str__(self):
		return str(self.mapping)

	def freeze(self):
		"""No longer allow modifications"""
		# TODO We should try to have a second class (FrozenNumericSampleIds) for
		# this or try to get rid of NumericSampleIds
		self.frozen = True

	def inverse_mapping(self):
		"""Returns a dict mapping numeric ids to sample names."""
		return {numeric_id:name for name, numeric_id in self.mapping.items()}
	
	def __getstate__(self):
		return (self.mapping, self.frozen)
	
	def __setstate__(self, state):
		mapping, frozen = state
		self.mapping = mapping
		self.frozen = frozen


cdef class Read:
	def __cinit__(self, str name = None, int mapq = 0, int source_id = 0, int sample_id = 0, int reference_start = -1, str BX_tag = None):
		cdef string _name = b''
		cdef string _BX_tag = b''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			if BX_tag is not '' and BX_tag is not None:
				_BX_tag = BX_tag.encode('UTF-8')
			self.thisptr = new cpp.Read(_name, mapq, source_id, sample_id, reference_start, _BX_tag)
			self.ownsptr = True

	def __dealloc__(self):
		if self.ownsptr:
			assert self.thisptr != NULL
			del self.thisptr

	def __repr__(self):
		assert self.thisptr != NULL
		return 'Read(name={!r}, mapq={}, source_id={}, sample_id={}, reference_start={},  BX_tag={}, variants={})'.format(
			self.name, self.mapqs, self.source_id, self.sample_id, self.reference_start, self.BX_tag, list(self))

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

	property sample_id:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getSampleID()
	
	property reference_start:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getReferenceStart()

	property BX_tag:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getBXTag().decode('utf-8')

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
	
	def __getstate__(self):
		mapqs = [mapq for mapq in self.mapqs]
		variants = [(var.position, var.allele, var.quality) for var in self]
		return (mapqs, self.name, self.source_id, self.sample_id, self.reference_start, self.BX_tag, variants)
	
	def __setstate__(self, state):
		mapqs, name, source_id, sample_id, reference_start, BX_tag, variants = state
		
		# TODO: Duplicated code from __cinit__ is ugly, but cinit cannot be used here directly
		cdef string _name = b''
		cdef string _BX_tag = b''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			if BX_tag is not b'' and BX_tag is not None:
				_BX_tag = BX_tag.encode('UTF-8')
			self.thisptr = new cpp.Read(_name, mapqs[0] if len(mapqs) > 0 else 0, source_id, sample_id, reference_start, _BX_tag)
			self.ownsptr = True

		for mapq in mapqs[1:]:
			self.add_mapq(mapq)
		for (pos, allele, quality) in variants:
			self.add_variant(pos, allele, quality)

	def add_variant(self, int position, int allele, int quality):
		assert self.thisptr != NULL
		self.thisptr.addVariant(position, allele, quality)

	def add_mapq(self, int mapq):
		assert self.thisptr != NULL
		self.thisptr.addMapq(mapq)

	def sort(self):
		assert self.thisptr != NULL
		self.thisptr.sortVariants()

	def is_sorted(self):
		assert self.thisptr != NULL
		return self.thisptr.isSorted()

	def has_BX_tag(self):
		assert self.thisptr != NULL
		return self.thisptr.hasBXTag()


cdef class ReadSet:
	def __cinit__(self):
		self.thisptr = new cpp.ReadSet()

	def __dealloc__(self):
		del self.thisptr

	def add(self, Read read):
		"""Adds a read to the set.
		WARNING: this will internally create a copy of the wrapped C++ Read object,
		so that subsequent changes to the Read don't affect the
		newly created copy that is added to the ReadSet."""
		self.thisptr.add(new cpp.Read(read.thisptr[0]))

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
		cdef string name = b''
		cdef cpp.Read* cread = NULL
		cdef Read read = Read()
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
	
	def __getstate__(self):
		return ([read for read in self])
	
	def __setstate__(self, state):
		self.thisptr = new cpp.ReadSet()
		for read in state:
			self.add(read)

	#def get_by_name(self, name):
		#cdef string _name = name.encode('UTF-8')
		#cdef cpp.Read* cread = self.thisptr.getByName(_name)
		#cdef Read read = Read()
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
		cdef cpp.IndexSet* index_set = new cpp.IndexSet()
		cdef int i
		for i in reads_to_select:
			index_set.add(i)
		result = ReadSet()
		del result.thisptr
		result.thisptr = self.thisptr.subset(index_set)
		del index_set
		return result

	def get_positions(self):
		cdef vector[unsigned int]* v = self.thisptr.get_positions()
		result = list(v[0])
		del v
		return result


cdef class PedigreeDPTable:
	def __cinit__(self, ReadSet readset, recombcost, Pedigree pedigree, bool distrust_genotypes = False, positions = None):
		"""Build the DP table from the given read set which is assumed to be sorted;
		that is, the variants in each read must be sorted by position and the reads
		in the read set must also be sorted (by position of their left-most variant).
		"""
		cdef vector[unsigned int]* c_positions = NULL
		if positions is not None:
			c_positions = new vector[unsigned int]()
			for pos in positions:
				c_positions.push_back(pos)
		self.thisptr = new cpp.PedigreeDPTable(readset.thisptr, recombcost, pedigree.thisptr, distrust_genotypes, c_positions)
		self.pedigree = pedigree

	def __dealloc__(self):
		del self.thisptr

	def get_super_reads(self):
		"""Obtain optimal-score haplotypes. Returns a triple (mother, father, child)
		IMPORTANT: The ReadSet given at construction time must not have been altered.
		DPTable retained a pointer to this set and will access it again. If it has
		been altered, behavior is undefined.
		TODO: Change that.
		"""
		cdef vector[cpp.ReadSet*]* read_sets = new vector[cpp.ReadSet*]()

		for i in range(len(self.pedigree)):
			read_sets.push_back(new cpp.ReadSet())
		transmission_vector_ptr = new vector[unsigned int]()
		self.thisptr.get_super_reads(read_sets, transmission_vector_ptr)
		
		results = []
		for i in range(read_sets.size()):
			rs = ReadSet()
			del rs.thisptr
			rs.thisptr = deref(read_sets)[i]
			results.append(rs)

		python_transmission_vector = list(transmission_vector_ptr[0])
		del transmission_vector_ptr
		return results, python_transmission_vector

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


cdef class Pedigree:
	def __cinit__(self, numeric_sample_ids):
		self.thisptr = new cpp.Pedigree()
		self.numeric_sample_ids = numeric_sample_ids

	def __dealloc__(self):
		del self.thisptr

	def add_individual(self, id, genotypes, genotype_likelihoods=None):
		cdef vector[cpp.Genotype*] gt_vector
		for gt in genotypes:
			gt_vector.push_back(new cpp.Genotype((<Genotype?>gt).thisptr[0]))
		cdef vector[cpp.PhredGenotypeLikelihoods*] gl_vector
		if genotype_likelihoods:
			for gl in genotype_likelihoods:
				if gl is None:
					gl_vector.push_back(NULL)
				else:
					gl_vector.push_back(new cpp.PhredGenotypeLikelihoods((<PhredGenotypeLikelihoods?>gl).thisptr[0]) )
		else:
			for _ in genotypes:
				gl_vector.push_back(NULL)
		self.thisptr.addIndividual(self.numeric_sample_ids[id], gt_vector, gl_vector)

	def add_relationship(self, father_id, mother_id, child_id):
		self.thisptr.addRelationship(self.numeric_sample_ids[father_id], self.numeric_sample_ids[mother_id], self.numeric_sample_ids[child_id])

	property variant_count:
		"""Number of variants stored for each individual."""
		def __get__(self):
			return self.thisptr.get_variant_count()

	def genotype(self, sample_id, unsigned int variant_index):
		cdef const cpp.Genotype* gt = self.thisptr.get_genotype_by_id(self.numeric_sample_ids[sample_id], variant_index)
		return Genotype(gt[0].as_vector())

	def genotype_likelihoods(self, sample_id, unsigned int variant_index):
		cdef const cpp.PhredGenotypeLikelihoods* gl = self.thisptr.get_genotype_likelihoods_by_id(self.numeric_sample_ids[sample_id], variant_index)
		if gl == NULL:
			return None
		else:
			return PhredGenotypeLikelihoods(gl[0].as_vector(), gl[0].get_ploidy(), gl[0].get_nr_alleles())

	def __len__(self):
		return self.thisptr.size()

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')


cdef class PhredGenotypeLikelihoods:
	def __cinit__(self, vector[double] gl, unsigned int ploidy=2, unsigned int nr_alleles=2):
		self.thisptr = new cpp.PhredGenotypeLikelihoods(gl, ploidy, nr_alleles)

	def __dealloc__(self):
		del self.thisptr

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')

	def __getitem__(self, Genotype genotype):
		assert self.thisptr != NULL
		assert genotype.is_diploid_and_biallelic()
		return self.thisptr.get(genotype.thisptr[0])

	def __len__(self):
		return self.thisptr.size()

	def __iter__(self):
		for genotype in self.genotypes():
			yield self[genotype]
			
	def __eq__(self, PhredGenotypeLikelihoods other):
		if self.genotypes() != other.genotypes():
			return False
		for genotype in self.genotypes():
			if self[genotype] != other[genotype]:
				return False
		return True
		
	def genotypes(self):
		cdef vector[cpp.Genotype]* genotypes = new vector[cpp.Genotype]()
		self.thisptr.get_genotypes(deref(genotypes))
		result = [Genotype(genotype.as_vector()) for genotype in genotypes[0]]
		del genotypes
		return result
	
	
def binomial_coefficient(int n, int k):
	return cpp.binomial_coefficient(n, k)

			
cdef class Genotype:
	
	def __cinit__(self, vector[uint32_t] alleles):
		self.thisptr = new cpp.Genotype(alleles)

	def __dealloc__(self):
		del self.thisptr

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')
	
	def __repr__(self):
		return self.thisptr.toString().decode('utf-8')

	def is_none(self):
		return self.thisptr.is_none()

	def get_index(self):
		return self.thisptr.get_index()

	def as_vector(self):
		result = []
		cdef vector[uint32_t] alleles = self.thisptr.as_vector()
		for allele in alleles:
			result.append(allele)
		return alleles

	def is_homozygous(self):
		return self.thisptr.is_homozygous()
	
	def is_diploid_and_biallelic(self):
		return self.thisptr.is_diploid_and_biallelic()
	
	def get_ploidy(self):
		return self.thisptr.get_ploidy()

	def __eq__(self, Genotype g):
		return self.thisptr[0] == g.thisptr[0]

	def __ne__(self, Genotype g):
		return self.thisptr[0] != g.thisptr[0]
	
	def __lt__(self, Genotype g):
		return self.thisptr[0] < g.thisptr[0]

	def __getstate__(self):
		return (self.thisptr.get_index(), self.thisptr.get_ploidy())

	def __setstate__(self, state):
		index, ploidy = state
		cdef vector[uint32_t] alleles = cpp.convert_index_to_alleles(index, ploidy)
		if self.thisptr != NULL:
			del self.thisptr
		self.thisptr = new cpp.Genotype(alleles)

	def __deepcopy__(self, memo):
		return Genotype.__new__(Genotype, self.as_vector())

	def __hash__(self):
		return hash(self.thisptr.get_index())
	
	
def get_max_genotype_ploidy():
	return cpp.get_max_genotype_ploidy()


def get_max_genotype_alleles():
	return cpp.get_max_genotype_alleles()


cdef class GenotypeDPTable:
	def __cinit__(self, numeric_sample_ids, ReadSet readset, recombcost, Pedigree pedigree, positions = None):
		"""Build the DP table from the given read set which is assumed to be sorted;
		that is, the variants in each read must be sorted by position and the reads
		in the read set must also be sorted (by position of their left-most variant).
		"""
		cdef vector[unsigned int]* c_positions = NULL
		if positions is not None:
			c_positions = new vector[unsigned int]()
			for pos in positions:
				c_positions.push_back(pos)
		self.thisptr = new cpp.GenotypeDPTable(readset.thisptr, recombcost, pedigree.thisptr, c_positions)
		self.pedigree = pedigree
		self.numeric_sample_ids = numeric_sample_ids

	def __dealloc__(self):
		del self.thisptr

	def get_genotype_likelihoods(self, sample_id, unsigned int pos):
		return PhredGenotypeLikelihoods(self.thisptr.get_genotype_likelihoods(self.numeric_sample_ids[sample_id],pos))


def compute_genotypes(ReadSet readset, positions = None):
	cdef vector[cpp.Genotype]* genotypes_vector = new vector[cpp.Genotype]()
	cdef vector[cpp.GenotypeDistribution]* gl_vector = new vector[cpp.GenotypeDistribution]()
	cdef vector[unsigned int]* c_positions = NULL
	if positions is not None:
		c_positions = new vector[unsigned int]()
		for pos in positions:
			c_positions.push_back(pos)
	cpp.compute_genotypes(readset.thisptr[0], genotypes_vector, gl_vector, c_positions)
	# TODO: Inefficient
	genotypes = [Genotype(genotype.as_vector()) for genotype in genotypes_vector[0]]
	#genotypes = list(genotypes_vector[0])
	gls = [(gl_vector[0][i].probabilityOf(0), gl_vector[0][i].probabilityOf(1), gl_vector[0][i].probabilityOf(2)) for i in range(gl_vector[0].size())]
	del genotypes_vector
	del gl_vector
	return genotypes, gls


cdef class HapChatCore:
	def __cinit__(self, ReadSet readset):
		self.thisptr = new cpp.HapChatCore(readset.thisptr)
	def __dealloc__(self):
		del self.thisptr
	def get_length(self):
		return self.thisptr.get_length()
	def get_super_reads(self):
		cdef vector[cpp.ReadSet*]* read_sets = new vector[cpp.ReadSet*]()
		leng=self.thisptr.get_length()
		for i in range(leng):
			read_sets.push_back(new cpp.ReadSet())
		self.thisptr.get_super_reads(read_sets)
		
		results = []
		for i in range(read_sets.size()):
			rs = ReadSet()
			del rs.thisptr
			rs.thisptr = deref(read_sets)[i]
			results.append(rs)
		
		return results, None
	def get_optimal_cost(self):
		return self.thisptr.get_optimal_cost()
	def get_optimal_partitioning(self):
		cdef vector[bool]* p = self.thisptr.get_optimal_partitioning()
		result = ['*' for x in p[0]]
		del p
		return result
		
cdef class Caller:
	def __cinit__(self, string reference, int k, int window ):
		self.thisptr = new cpp.Caller(reference, k, window)
		
	def __dealloc__(self):
		del self.thisptr
		
	def all_variants(self, vector[pair[int,int]] variants_list):
		cdef deque[pair[int,int]] v_list
		for v in variants_list:
			v_list.push_back(v)
		self.thisptr.all_variants(v_list)

	def add_read(self, int bam_alignment_pos, vector[vector[int]]  bam_alignment_cigartuples, string bam_alignment_query, string outfile):
		self.thisptr.add_read(bam_alignment_pos, bam_alignment_cigartuples, bam_alignment_query, outfile)
	
	def final_pop(self, string outfile):
		self.thisptr.final_pop(outfile)
		
	def finish(self):
		pass
		
