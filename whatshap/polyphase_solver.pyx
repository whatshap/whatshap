# cython: language_level=3

from collections import defaultdict
import itertools as it

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.unordered_map cimport unordered_map
from libc.stdint cimport uint32_t
from . cimport cpp
cimport cython


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
    
    def scoreReadset(self, ReadSet readset, uint32_t minOverlap, uint32_t ploidy, double err):
        sim = TriangleSparseMatrix()
        self.thisptr.scoreReadset(sim.thisptr, readset.thisptr, minOverlap, ploidy, err)
        return sim


def scoreReadset(readset, minOverlap, ploidy, err=0.0):
    readscoring = ReadScoring()
    sim = readscoring.scoreReadset(readset, minOverlap, ploidy, err)
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
