# cython: language_level=3

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport int8_t, uint32_t, uint64_t
from libcpp.unordered_map cimport unordered_map

from whatshap cimport cpp
from whatshap.core cimport ReadSet

cdef class ClusterEditingSolver:
    def __cinit__(self, TriangleSparseMatrix m, bundleEdges):
        self.thisptr = new cpp.ClusterEditingSolver(m.thisptr[0], bundleEdges)
        self.m = m

    def run(self):
        cdef cpp.ClusterEditingSolution solution = self.thisptr.run()
        clusters = []
        n_clusters = solution.getNumClusters()
        for i in range(n_clusters):
            clusters.append(solution.getCluster(i))
        return clusters


cdef class AlleleMatrix:
    def __cinit__(self, ReadSet rs=None):
        if rs:
            self.thisptr = new cpp.AlleleMatrix(rs.thisptr)
        else:
            self.thisptr = NULL

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
        
    def getNumPositions(self):
        return self.thisptr.getNumPositions()
        
    def getPositions(self):
        return self.thisptr.getPositions()

    def getAllele(self, uint32_t readId, uint32_t position):
        return self.thisptr.getAllele(readId, position)

    def getAlleleGlobal(self, uint32_t readId, uint32_t position):
        return self.thisptr.getAlleleGlobal(readId, position)

    def getRead(self, uint32_t readId):
        return self.thisptr.getRead(readId)

    def getFirstPos(self, uint32_t readId):
        return self.thisptr.getFirstPos(readId)

    def getLastPos(self, uint32_t readId):
        return self.thisptr.getLastPos(readId)
    
    def getGlobalId(self, uint32_t readId):
        return self.thisptr.getGlobalId(readId)

    def globalToLocal(self, uint32_t position):
        return self.thisptr.globalToLocal(position)
        
    def localToGlobal(self, uint32_t position):
        return self.thisptr.localToGlobal(position)
    
    def getAlleleDepths(self, uint32_t position):
        return self.thisptr.getAlleleDepths(position)
        
    def extractInterval(self, uint32_t start, uint32_t end, bool removeEmpty=True):
        cdef AlleleMatrix mx = AlleleMatrix.__new__(AlleleMatrix)
        mx.thisptr = self.thisptr.extractInterval(start, end, removeEmpty)
        return mx
    
    def extractSubMatrix(self, vector[uint32_t] positions, vector[uint32_t] readIds, bool removeEmpty=True):
        cdef AlleleMatrix mx = AlleleMatrix.__new__(AlleleMatrix)
        mx.thisptr = self.thisptr.extractSubMatrix(positions, readIds, removeEmpty)
        return mx

    def __iter__(self):
        for i in range(self.thisptr.size()):
            yield self.getRead(i)

    def __len__(self):
        return self.thisptr.size()

    def __getstate__(self):
        read_list = [{pos: allele for pos, allele in read} for read in self]
        pos_list = self.getPositions()
        id_list = [self.getGlobalId(i) for i in range(len(self))]
        return read_list, pos_list, id_list

    def __setstate__(self, state):
        read_list, pos_list, id_list = state
        self.thisptr = new cpp.AlleleMatrix(read_list, pos_list, id_list)


cdef class TriangleSparseMatrix:
    def __cinit__(self):
        self.thisptr = new cpp.TriangleSparseMatrix()

    def __dealloc__(self):
        del self.thisptr

    def get(self, int i, int j):
        return self.thisptr.get(i, j)
    
    def set(self, int i, int j, float v):
        return self.thisptr.set(i, j, v)

    def size(self):
        return self.thisptr.size()

    def getEntries(self):
        return self.thisptr.getEntries()

    def __iter__(self):
        pairs = self.thisptr.getEntries()
        for i in range(self.thisptr.size()):
            yield pairs[i]

    def __len__(self):
        return self.thisptr.size()


cdef class ReadScoring:
    def __cinit__(self):
        self.thisptr = new cpp.ReadScoring()
    
    def scoreReadset(self, AlleleMatrix am, uint32_t minOverlap, uint32_t ploidy, double err):
        sim = TriangleSparseMatrix()
        self.thisptr.scoreReadset(sim.thisptr, am.thisptr, minOverlap, ploidy, err)
        return sim


def scoreReadset(am, minOverlap, ploidy, err=0.0):
    readscoring = ReadScoring()
    sim = readscoring.scoreReadset(am, minOverlap, ploidy, err)
    del readscoring
    return sim
    
    
cdef class HaploThreader:
    def __cinit__(self, ploidy, switchCost, affineSwitchCost, maxClusterGap, rowLimit):
        self.thisptr = new cpp.HaploThreader(ploidy, switchCost, affineSwitchCost, maxClusterGap, rowLimit)
        
    def computePathsBlockwise(self, vector[uint32_t]& blockStarts, vector[vector[uint32_t]]& covMap, vector[unordered_map[uint32_t, unordered_map[uint32_t, uint32_t]]]& alleleDepths):
        cdef vector[vector[uint32_t]] path
        path = self.thisptr.computePaths(blockStarts, covMap, alleleDepths)
        
        # convert to python data structure
        py_path = []
        path_len = path.size()
        if path_len > 0:
            ploidy = path[0].size()
            for i in range(path_len):
                pos = []
                for j in range(ploidy):
                    pos.append(path[i][j])
                py_path.append(pos)
        
        return py_path

    def computePaths(self, uint32_t start, uint32_t end, vector[vector[uint32_t]]& covMap, vector[unordered_map[uint32_t, unordered_map[uint32_t, uint32_t]]]& alleleDepths):
        cdef vector[vector[uint32_t]] path
        path = self.thisptr.computePaths(start, end, covMap, alleleDepths)
        
        # convert to python data structure
        py_path = []
        path_len = path.size()
        if path_len > 0:
            ploidy = path[0].size()
            for i in range(path_len):
                pos = []
                for j in range(ploidy):
                    pos.append(path[i][j])
                py_path.append(pos)
        
        return py_path

cdef class SwitchFlipCalculator:
    def __cinit__(self, ploidy, switch_cost=1, flip_cost=1):
        self.thisptr = new cpp.SwitchFlipCalculator(ploidy, switch_cost, flip_cost)
        self.ploidy = ploidy
    
    def compute_switch_flips_poly(self, phasing0, phasing1):
        assert len(phasing0) == len(phasing1) == self.ploidy
        assert self.ploidy >= 2
        assert len(phasing0[0]) > 0

        # convert haplotype-wise phasing over strings to position-wise phasing over int
        num_vars = len(phasing0[0])
        input0 = [[int(phasing0[k][i]) for k in range(self.ploidy)] for i in range(num_vars)]
        input1 = [[int(phasing1[k][i]) for k in range(self.ploidy)] for i in range(num_vars)]
        
        # create result space for error counts and detailed configuration
        cdef vector[uint32_t] switches_in_column
        cdef vector[vector[uint32_t]] flips_in_column
        cdef vector[vector[uint32_t]] perm_in_column
        
        # run compare algorithm
        cdef pair[double, double] result
        result = self.thisptr.compare(input0, input1, switches_in_column, flips_in_column, perm_in_column)
        
        backtracking_info = True
        if switches_in_column.size() == 0 or flips_in_column.size() == 0 or perm_in_column.size() == 0:
            backtracking_info = False
        else:
            assert num_vars == switches_in_column.size() or switches_in_column.size() == 0
            assert num_vars == flips_in_column.size() or flips_in_column.size() == 0
            assert num_vars == perm_in_column.size() or perm_in_column.size() == 0
        
        # convert to python datastructures
        py_switches_in_column = []
        py_flips_in_column = []
        py_perm_in_column = []
        if backtracking_info:
            for i in range(num_vars):
                py_switches_in_column.append(switches_in_column[i])
                flip = []
                for k in range(len(flips_in_column[i])):
                    flip.append(flips_in_column[i][k])
                py_flips_in_column.append(flip)
                perm = []
                for k in range(self.ploidy):
                    perm.append(perm_in_column[i][k])
                py_perm_in_column.append(perm)

        return result.first, result.second, py_switches_in_column, py_flips_in_column, py_perm_in_column


cdef class ProgenyGenotypeLikelihoods:
    def __cinit__(self, ploidy, numSamples, numPositions):
        self.thisptr = new cpp.ProgenyGenotypeLikelihoods(ploidy, numSamples, numPositions)

    def __dealloc__(self):
        del self.thisptr

    def getGl(self, uint32_t pos, uint32_t sample, uint32_t genotype):
        return self.thisptr.getGl(pos, sample, genotype)
    
    def getGlv(self, uint32_t pos, uint32_t sample):
        return self.thisptr.getGlv(pos, sample)
    
    def setGl(self, uint32_t pos, uint32_t sample, uint32_t genotype, double l):
        self.thisptr.setGl(pos, sample, genotype, l)
    
    def setGlv(self, uint32_t pos, uint32_t sample, vector[double] l):
        self.thisptr.setGlv(pos, sample, l)

    def getPloidy(self):
        return self.thisptr.getPloidy()
    
    def getNumSamples(self):
        return self.thisptr.getNumSamples()
    
    def getNumPositions(self):
        return self.thisptr.getNumPositions()
    
    def getSimplexNulliplexScore(self, uint32_t pos1, uint32_t pos2):
        return self.thisptr.getSimplexNulliplexScore(pos1, pos2)
    
    def getSimplexSimplexScore(self, uint32_t pos1, uint32_t pos2):
        return self.thisptr.getSimplexSimplexScore(pos1, pos2)
    
    def getDuplexNulliplexScore(self, uint32_t pos1, uint32_t pos2):
        return self.thisptr.getDuplexNulliplexScore(pos1, pos2)

    def __len__(self):
        return self.thisptr.getNumPositions()


