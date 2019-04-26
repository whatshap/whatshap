# kate: syntax python
"""
Declarations for all C++ classes that are wrapped from Cython.
"""
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t
from libcpp.unordered_map cimport unordered_map

cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int, int, int, int, string) except +
		Read(Read) except +
		string toString() except +
		void addVariant(int, int, vector[unsigned int]) except +
		string getName() except +
		vector[int] getMapqs() except +
		void addMapq(int) except +
		int getPosition(int) except +
		void setPosition(int, int)  except +
		int getAllele(int) except +
		void setAllele(int, int) except +
		vector[unsigned int] getVariantQuality(int) except +
		void setVariantQuality(int, vector[unsigned int]) except +
		int getVariantCount() except +
		void sortVariants() except +
		bool isSorted() except +
		int getSourceID() except +
		int getSampleID() except +
		int getReferenceStart() except +
		string getBXTag() except +
		bool hasBXTag() except +


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


cdef extern from "../src/pedigree.h":
	cdef cppclass Pedigree:
		Pedigree(unsigned int ploidy) except +
		void addIndividual(unsigned int id, vector[Genotype*] genotypes, vector[GenotypeLikelihoods*]) except +
		void addRelationship(unsigned int f, unsigned int m, unsigned int c) except +
		unsigned int size()
		string toString() except +
		const Genotype* get_genotype_by_id(unsigned int, unsigned int) except +
		const GenotypeLikelihoods* get_genotype_likelihoods_by_id(unsigned int, unsigned int) except +
		unsigned int get_variant_count() except +
		unsigned int triple_count() except +


cdef extern from "../src/pedigreedptable.h":
	cdef cppclass PedigreeDPTable:
		PedigreeDPTable(ReadSet*, vector[unsigned int], Pedigree* pedigree, unsigned int ploidy, bool distrust_genotypes, vector[unsigned int]* allele_counts, vector[unsigned int]* positions, vector[unsigned int]* partitioning) except +
		void get_super_reads(vector[ReadSet*]*, vector[unsigned int]* transmission_vector) except +
		int get_optimal_score() except +
#		vector[unsigned int]* get_block_boundaries() except +
		vector[unsigned int]* get_optimal_partitioning()

cdef extern from "../src/genotypedptable.h":
	cdef cppclass GenotypeDPTable:
		GenotypeDPTable(ReadSet*, vector[unsigned int], Pedigree* pedigree, unsigned int ploidy, vector[unsigned int]* positions) except +
		vector[long double] get_genotype_likelihoods(unsigned int individual, unsigned int position) except +

cdef extern from "../src/genotypelikelihoods.h":
	cdef cppclass GenotypeLikelihoods:
		GenotypeLikelihoods() except +
		GenotypeLikelihoods(unsigned int ploidy, unsigned int n_alleles, vector[double], bool is_phred) except +
		GenotypeLikelihoods(GenotypeLikelihoods) except +
		double get(Genotype) except +
		unsigned int genotype_count() except +
		vector[double] as_vector() except +
		string toString() except +
		unsigned int get_ploidy() except +
		unsigned int get_n_alleles() except +
		Genotype get_likeliest_genotype(double threshold_prob) except +
		bool is_phred() except +

cdef extern from "../src/genotype.h":
	cdef cppclass Genotype:
		Genotype() except +
		Genotype(vector[unsigned int]) except +
		Genotype(Genotype) except +
		void add_allele(unsigned int) except +
		vector[unsigned int] as_vector() except +
		bool is_none() except +
		unsigned int get_index(unsigned int, unsigned int) except +
		string toString() except +
		bool is_homozygous() except +
	cdef bool operator==(Genotype,Genotype) except +
	cdef bool operator!=(Genotype,Genotype) except +

cdef extern from "../src/genotypedistribution.h":
	cdef cppclass GenotypeDistribution:
		GenotypeDistribution(double hom_ref_prob, double het_prob, double hom_alt_prob) except +
		double probabilityOf(Genotype genotype) except +

cdef extern from "../src/genotyper.h":
	void compute_genotypes(ReadSet, vector[Genotype]* genotypes, vector[GenotypeDistribution]* genotype_likelihoods, vector[unsigned int]* positions)  except +

cdef extern from "../src/clusterediting/DynamicSparseGraph.h":
	cdef cppclass DynamicSparseGraph:
		DynamicSparseGraph(uint32_t) except +
		void addEdge(uint32_t, uint32_t, double) except +
		void setWeight(uint32_t, uint32_t, double) except +
		void clearAndResize(uint32_t) except +
		
cdef extern from "../src/clusterediting/CoreAlgorithm.h":
	cdef cppclass CoreAlgorithm:
		CoreAlgorithm(DynamicSparseGraph graph, bool bundleEdges) except +
		ClusterEditingSolutionLight run() except +

cdef extern from "../src/clusterediting/ClusterEditingSolutionLight.h":
	cdef cppclass ClusterEditingSolutionLight:
		ClusterEditingSolutionLight() except +
		ClusterEditingSolutionLight(ClusterEditingSolutionLight) except +
		ClusterEditingSolutionLight(double pTotalCost, vector[vector[int]] pClusters) except +
		vector[unsigned int] getCluster(int i) except +
		double getTotalCost() except +
		int getNumClusters() except +

cdef extern from "../src/clusterediting/TriangleSparseMatrix.h":
	cdef cppclass TriangleSparseMatrix:
		TriangleSparseMatrix() except +
		unsigned int entryToIndex(unsigned int i, unsigned int j) except +
		unsigned int size() except +
		float get(unsigned int i, unsigned int j) except +
		void set(unsigned int i, unsigned int j, float v) except +
		vector[pair[uint32_t, uint32_t]] getEntries() except +

cdef extern from "../src/clusterediting/ReadScoring.h":
	cdef cppclass ReadScoring:
		ReadScoring() except +
		void scoreReadsetGlobal(TriangleSparseMatrix* result, ReadSet* readset, uint32_t minOverlap,uint32_t ploidy) except +
		void scoreReadsetLocal(TriangleSparseMatrix* result, ReadSet* readset, uint32_t minOverlap, uint32_t ploidy) except +
		void scoreReadsetPatterns(TriangleSparseMatrix* result, ReadSet* readset, uint32_t minOverlap, uint32_t ploidy, double errorrate, uint32_t windowSize) except +
		
cdef extern from "../src/threading/HaploThreader.h":
	cdef cppclass HaploThreader:
		HaploThreader(uint32_t ploidy, double switchCost, double affineSwitchCost, bool symmetryOptimization, uint32_t rowLimit) except +
		vector[vector[uint32_t]] computePaths(uint32_t start, uint32_t end,
					vector[vector[uint32_t]]& covMap,
                    vector[vector[double]]& coverage, 
                    vector[vector[uint32_t]]& consensus,
                    vector[unordered_map[uint32_t, uint32_t]]& genotypes,
					vector[vector[vector[double]]]& clusterDissim) except +
		vector[vector[uint32_t]] computePaths(vector[uint32_t]& blockStarts,
					vector[vector[uint32_t]]& covMap,
                    vector[vector[double]]& coverage, 
                    vector[vector[uint32_t]]& consensus,
                    vector[unordered_map[uint32_t, uint32_t]]& genotypes,
					vector[vector[vector[double]]]& clusterDissim) except +