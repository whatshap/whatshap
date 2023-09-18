# cython: language_level=3

"""
Declarations for all C++ classes that are wrapped from Cython.
"""
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport int8_t, uint32_t, uint64_t
from libcpp.unordered_map cimport unordered_map
from libcpp.deque cimport deque


cdef extern from "read.h":
    cdef cppclass Read:
        Read(string, int, int, int, int, string) except +
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
        int getReferenceStart() except +
        string getBXTag() except +
        bool hasBXTag() except +


cdef extern from "indexset.h":
    cdef cppclass IndexSet:
        IndexSet() except +
        bool contains(int) except +
        void add(int) except +
        int size() except +
        string toString() except +


cdef extern from "readset.h":
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


cdef extern from "pedigree.h":
    cdef cppclass Pedigree:
        Pedigree() except +
        void addIndividual(unsigned int id, vector[Genotype*] genotypes, vector[PhredGenotypeLikelihoods*]) except +
        void addRelationship(unsigned int f, unsigned int m, unsigned int c) except +
        unsigned int size()
        string toString() except +
        const Genotype* get_genotype_by_id(unsigned int, unsigned int) except +
        const PhredGenotypeLikelihoods* get_genotype_likelihoods_by_id(unsigned int, unsigned int) except +
        unsigned int get_variant_count() except +
        unsigned int triple_count() except +


cdef extern from "pedigreedptable.h":
    cdef cppclass PedigreeDPTable:
        PedigreeDPTable(ReadSet*, vector[unsigned int], Pedigree* pedigree, bool distrust_genotypes, vector[unsigned int]* positions) except +
        void get_super_reads(vector[ReadSet*]*, vector[unsigned int]* transmission_vector) except +
        int get_optimal_score() except +
        vector[bool]* get_optimal_partitioning()
        
        
cdef extern from "binomial.h":
    cdef int binomial_coefficient(int n, int k) except +
        
        
cdef extern from "genotype.h":
    cdef cppclass Genotype:
        Genotype() except +
        Genotype(vector[uint32_t]) except +
        Genotype(Genotype) except +
        vector[uint32_t] as_vector() except +
        bool is_none() except +
        uint64_t get_index() except +
        string toString() except +
        bool is_homozygous() except +
        bool is_diploid_and_biallelic() except +
        uint32_t get_ploidy() except +
    cdef bool operator==(Genotype,Genotype) except +
    cdef bool operator!=(Genotype,Genotype) except +
    cdef bool operator<(Genotype,Genotype) except +
    cdef vector[uint32_t] convert_index_to_alleles(uint64_t index, uint32_t ploidy) except +
    cdef uint32_t get_max_genotype_ploidy() except +
    cdef uint32_t get_max_genotype_alleles() except +


cdef extern from "genotypedptable.h":
    cdef cppclass GenotypeDPTable:
        GenotypeDPTable(ReadSet*, vector[unsigned int], Pedigree* pedigree, vector[unsigned int]* positions) except +
        vector[long double] get_genotype_likelihoods(unsigned int individual, unsigned int position) except +

cdef extern from "phredgenotypelikelihoods.h":
    cdef cppclass PhredGenotypeLikelihoods:
        PhredGenotypeLikelihoods(vector[double], unsigned int, unsigned int) except +
        PhredGenotypeLikelihoods(PhredGenotypeLikelihoods) except +
        double get(Genotype) except +
        string toString() except +
        unsigned int get_ploidy() except +
        unsigned int get_nr_alleles() except +
        unsigned int size() except +
        vector[double] as_vector() except +
        void get_genotypes(vector[Genotype]&) except +


cdef extern from "genotypedistribution.h":
    cdef cppclass GenotypeDistribution:
        GenotypeDistribution(double hom_ref_prob, double het_prob, double hom_alt_prob) except +
        double probabilityOf(unsigned int genotype) except +


cdef extern from "genotyper.h":
    void compute_genotypes(ReadSet, vector[Genotype]* genotypes, vector[GenotypeDistribution]* genotype_likelihoods, vector[unsigned int]* positions)  except +


cdef extern from "hapchat/hapchatcore.cpp":
    cdef cppclass HapChatCore:
        HapChatCore(ReadSet*)
        void get_super_reads(vector[ReadSet*]*)
        vector[bool]* get_optimal_partitioning()
        int get_length()
        int get_optimal_cost()


cdef extern from "polyphase/clustereditingsolver.h":
    cdef cppclass ClusterEditingSolver:
        ClusterEditingSolver(TriangleSparseMatrix m, bool bundleEdges) except +
        ClusterEditingSolution run() except +


cdef extern from "polyphase/clustereditingsolution.h":
    cdef cppclass ClusterEditingSolution:
        ClusterEditingSolution() except +
        ClusterEditingSolution(ClusterEditingSolution) except +
        ClusterEditingSolution(double pTotalCost, vector[vector[int]] pClusters) except +
        vector[unsigned int] getCluster(int i) except +
        double getTotalCost() except +
        int getNumClusters() except +


cdef extern from "polyphase/allelematrix.h":
    cdef cppclass AlleleMatrix:
        AlleleMatrix(ReadSet* rs) except +
        AlleleMatrix(vector[unordered_map[uint32_t, int8_t]]& readList, vector[uint32_t]& posList, vector[uint32_t]& idList) except +
        uint64_t size() except +
        uint64_t getNumPositions() except +
        vector[uint32_t] getPositions() except +
        int8_t getAllele(uint32_t readId, uint32_t position) except +
        int8_t getAlleleGlobal(uint32_t readId, uint32_t genPosition) except +
        vector[pair[uint32_t, int8_t]] getRead(uint32_t readId) except +
        uint32_t getFirstPos(uint32_t readId) except +
        uint32_t getLastPos(uint32_t readId) except +
        uint32_t getGlobalId(uint32_t readId) except +
        uint32_t globalToLocal(uint32_t genPosition) except +
        uint32_t localToGlobal(uint32_t position) except +
        vector[uint32_t] getAlleleDepths(uint32_t position) except +
        AlleleMatrix* extractInterval(uint32_t start, uint32_t end, bool removeEmpty) except +
        AlleleMatrix* extractSubMatrix(vector[uint32_t] positions, vector[uint32_t] readIds, bool removeEmpty) except +


cdef extern from "polyphase/trianglesparsematrix.h":
    cdef cppclass TriangleSparseMatrix:
        TriangleSparseMatrix() except +
        unsigned int entryToIndex(unsigned int i, unsigned int j) except +
        unsigned int size() except +
        float get(unsigned int i, unsigned int j) except +
        void set(unsigned int i, unsigned int j, float v) except +
        vector[pair[uint32_t, uint32_t]] getEntries() except +


cdef extern from "polyphase/readscoring.h":
    cdef cppclass ReadScoring:
        ReadScoring() except +
        void scoreReadset(TriangleSparseMatrix* result, AlleleMatrix* readset, uint32_t minOverlap, uint32_t ploidy, double err) except +


cdef extern from "polyphase/haplothreader.h":
    cdef cppclass HaploThreader:
        HaploThreader(uint32_t ploidy, double switchCost, double affineSwitchCost, uint32_t maxClusterGap, uint32_t rowLimit) except +
        vector[vector[uint32_t]] computePaths(uint32_t start, uint32_t end,
                    vector[vector[uint32_t]]& covMap,
                    vector[unordered_map[uint32_t, unordered_map[uint32_t, uint32_t]]]& alleleDepths) except +
        vector[vector[uint32_t]] computePaths(vector[uint32_t]& blockStarts,
                    vector[vector[uint32_t]]& covMap,
                    vector[unordered_map[uint32_t, unordered_map[uint32_t, uint32_t]]]& alleleDepths) except +


cdef extern from "polyphase/switchflipcalculator.h":
    cdef cppclass SwitchFlipCalculator:
        SwitchFlipCalculator(uint32_t ploidy, double switchCost, double flipCost) except +
        pair[double, double] compare(vector[vector[uint32_t]]& phasing0,
                    vector[vector[uint32_t]]& phasing1,
                    vector[uint32_t]& switchesInColumn,
                    vector[vector[uint32_t]]& flippedHapsInColumn,
                    vector[vector[uint32_t]]& permInColumn) except +


cdef extern from "polyphase/progenygenotypelikelihoods.h":
    cdef cppclass ProgenyGenotypeLikelihoods:
        ProgenyGenotypeLikelihoods(uint32_t ploidy, uint32_t numSamples, uint32_t numPositions) except +
        double getGl(uint32_t pos, uint32_t sampleId, uint32_t genotype) except +
        vector[double] getGlv(uint32_t pos, uint32_t sampleId) except +
        void setGl(uint32_t pos, uint32_t sampleId, uint32_t genotype, double l) except +
        void setGlv(uint32_t pos, uint32_t sampleId, vector[double] l) except +
        uint32_t getPloidy() except +
        uint32_t getNumSamples() except +
        uint32_t getNumPositions() except +
        double getSimplexNulliplexScore(uint32_t pos1, uint32_t pos2) except +
        double getSimplexSimplexScore(uint32_t pos1, uint32_t pos2) except +
        double getDuplexNulliplexScore(uint32_t pos1, uint32_t pos2) except +
        
    
cdef extern from "../src/caller.h":
    ctypedef pair[int,int] cpp_pair
    ctypedef unordered_map[int,int] cpp_map
    ctypedef pair[cpp_pair,cpp_map] cpp_Pair
    cdef cppclass Caller:
        Caller(string&, int, int) except +
        void all_variants(deque[pair[int,int]]) except +
        void add_read(int,vector[vector[int]], string, string) except +
        void finish() except +
        pair[int,int] get_column(int) except +
        void pop_column() except +
        void process_complete_columns(int, string) except +
        void advance_to(int) except +
        void enumerate_reference_kmers(string&, int) except +
        void final_pop(string) except +
        void enumerate_kmers(int, string, int, vector[vector[int]]) except +
