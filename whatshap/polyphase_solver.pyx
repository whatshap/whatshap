from collections import defaultdict
import itertools as it
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.unordered_map cimport unordered_map
from libc.stdint cimport uint32_t
cimport cpp
cimport cython


cdef class DynamicSparseGraph:
	def __cinit__(self, uint32_t initialNumNodes):
		self.thisptr = new cpp.DynamicSparseGraph(initialNumNodes)
	def addEdge(self, uint32_t node_id1, uint32_t node_id2, double weight):
		self.thisptr.addEdge(node_id1, node_id2, weight)
	def setWeight(self, uint32_t node_id1, uint32_t node_id2, double weight):
		self.thisptr.setWeight(node_id1, node_id2, weight)
	def clearAndResize(self, uint32_t initialNumNodes):
		self.thisptr.clearAndResize(initialNumNodes)


cdef class ClusterEditingSolver:
	def __cinit__(self, DynamicSparseGraph graph, bundleEdges):
		self.thisptr = new cpp.ClusterEditingSolver(graph.thisptr[0], bundleEdges)
		self.graph = graph
	def run(self):
		cdef cpp.ClusterEditingSolutionLight solution = self.thisptr.run()
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

	def scoreReadsetGlobal(self, ReadSet readset, uint32_t minOverlap, uint32_t ploidy):
		sim = TriangleSparseMatrix()
		self.thisptr.scoreReadsetGlobal(sim.thisptr, readset.thisptr, minOverlap, ploidy)
		return sim
	
	def scoreReadsetLocal(self, ReadSet readset, vector[vector[uint32_t]] refHaplotypes, uint32_t minOverlap, uint32_t ploidy):
		sim = TriangleSparseMatrix()
		self.thisptr.scoreReadsetLocal(sim.thisptr, readset.thisptr, refHaplotypes, minOverlap, ploidy)
		return sim
	
	
cdef class HaploThreader:
	def __cinit__(self, ploidy, switchCost, affineSwitchCost, symmetryOptimization, rowLimit):
		self.thisptr = new cpp.HaploThreader(ploidy, switchCost, affineSwitchCost, symmetryOptimization, rowLimit)
		
	def computePathsBlockwise(self, vector[uint32_t]& blockStarts, vector[vector[uint32_t]]& covMap, vector[vector[double]]& coverage, vector[vector[uint32_t]]& consensus, vector[unordered_map[uint32_t, uint32_t]]& genotypes, vector[vector[vector[double]]]& clusterDissim):
		cdef vector[vector[uint32_t]] path
		path = self.thisptr.computePaths(blockStarts, covMap, coverage, consensus, genotypes, clusterDissim)
		
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

	def computePaths(self, uint32_t start, uint32_t end, vector[vector[uint32_t]]& covMap, vector[vector[double]]& coverage, vector[vector[uint32_t]]& consensus, vector[unordered_map[uint32_t, uint32_t]]& genotypes, vector[vector[vector[double]]]& clusterDissim):
		cdef vector[vector[uint32_t]] path
		path = self.thisptr.computePaths(start, end, covMap, coverage, consensus, genotypes, clusterDissim)
		
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
	
	
# Cython implementation of the DP threading
cdef subsetting(num_vars, clustering, coverage, positions, cov_map, ploidy, consensus, geno_map):

	cdef int num_clusters = len(positions)
	cdef unordered_map[int, pair[int,int]] column
	cdef unordered_map[int, int] pred
	cdef vector[vector[pair[int,int]]] scoring
	cdef vector[pair[int,int]] newcolumn
	
	#initialize first column
	c_tups = geno_map[0]
	conf_tups = list(it.chain.from_iterable(sorted(list(set(it.permutations(x)))) for x in c_tups))
	for tup in conf_tups:
		#for the first column, only coverage costs are computed, the predecessor is set to -1
		newcolumn.push_back((cov_costs(tup, 0, coverage), -1))
	scoring.push_back(newcolumn)
	
	cdef int var
	cdef float mininum
	cdef int minimum_index
	cdef bool min_exists = False
	cdef int avg_length_before = 0
	cdef int avg_length_after = 0
	for var in range(1,num_vars):
		print("computing variant %s of %s variants in total " % (var, num_vars))
		#every column is a vector containing a pair of the costs at this position and the position of the predecessor (used for backtracing) 
		newcolumn.clear()
		#find the suitable clusters that cover position var and compute the list of cluster tuples
		c_tuples = geno_map[var]
		conform_tups = list(it.chain.from_iterable(sorted(list(set(it.permutations(x)))) for x in c_tuples))		
		for tup in conform_tups:
			min_exists = False
			pred.clear()
			minimum = 1000000000
			minimum_index = 0
			#for the previous column, compute the list of all possible tuples of clusters that appear at position var-1
			pred_tups = geno_map[var-1]
			conf_pred_tups = list(it.chain.from_iterable(sorted(list(set(it.permutations(x)))) for x in pred_tups))	
			#compute the minimum of the previous column plus costs for switching:
			for pred_tup_i in range(len(conf_pred_tups)):	
				#the costs for the previous column are added to costs for switching from previous tuples to the current one
				pred[pred_tup_i] = (scoring[var-1][pred_tup_i].first+switch_costs(conf_pred_tups[pred_tup_i], tup,positions,var-1, ploidy))
				#find the  minimum in the previous column
				if pred[pred_tup_i] < minimum:
					min_exists = True
					minimum = pred[pred_tup_i]
					minimum_index = pred_tup_i
			#fill the matrix position with the computed costs and the index that was used from the column before (for simplifying backtracing)
			if min_exists:	
				newcolumn.push_back(((cov_costs(tup, var, coverage) + minimum), minimum_index))
			else:
				#no minimum exists: this is not expected to occur
				newcolumn.push_back((cov_costs(tup, var, coverage), 0))			
		scoring.push_back(newcolumn)
	#convert cython into python data structure
	scoring_res = []	
	for i in range(num_vars): 
		newcolumn_res = []
		for j in range(len(scoring[i])):
			newcolumn_res.append(scoring[i][j])		
		scoring_res.append(newcolumn_res)
	return(scoring_res)

#computes the genotype that would belong to the given tuple tup at position var

cdef int compute_tuple_genotype(consensus,tup,int var):
	cdef int genotype = 0
	cdef int i
	for i in tup:
		#allele = consensus[i][var]
		allele = consensus[var][i]
		genotype += allele
	return(genotype)
	

def compute_tuple_genotype_soft(consensus,tup, var, geno):
	genotype = 0
	for i in tup:
		allele = consensus[var][i]
		#allele = consensus[i][var]
		genotype += allele
	res = max((geno-genotype),(genotype-geno))
	return(res)

#computes costs for switching between two cluster tuples c_tuple1 and c_tuple2 at position var
cdef int switch_costs(c_tuple1, c_tuple2, positions,int var,int ploidy):
	cdef int costs = 0
	#switch costs depend on the position: if var is the end of c_tuple1 or var+1 is the beginning of c_tuple2, switching is free
	cdef int i
	cdef int errors = 0
	for i in range(0,ploidy):
	#	start = positions[c_tuple2[i]][0]
	#	end = positions[c_tuple1[i]][1]
	#	if (var != end and var+1 != start and (c_tuple1[i] != c_tuple2[i])):
		if (c_tuple1[i] != c_tuple2[i]):
			costs += 32
			#errors += 1
	#costs = errors + 8
	return(costs)

#computes the costs for differences between expected copy number (due to coverage) and the real copy number
#TODO: change 'hard' cutoffs to probability function
#TODO: does not work for a general <ploidy> yet
cdef int cov_costs(c_tuple, int var, coverage):
	#print(str(c_tuple)+" "+str(var))
	cdef int costs = 0
	cdef int exp_cn = 0
	cdef int i = 0
	#compute copy numbers for every cluster in c_tuple
	for i in range(0,4):	
#	for i in range(0,2):
		cov = coverage[var][c_tuple[i]] if c_tuple[i] in coverage[var] else 0
		#if cluster does not cover the position var:
		if (cov == 0):
			return (1000000000)
		#else compare the expected copy number to the real one
		else:
			if (cov > 0 and cov < 0.125):
#			if (cov > 0 and cov < 0.2):
				exp_cn = 0
			if (cov >= 0.125 and cov < 0.375):
#			if (cov >= 0.2 and cov < 0.4):				
				exp_cn = 1
			if (cov >= 0.375 and cov < 0.625):
#			if (cov >= 0.4 and cov < 0.6):
				exp_cn = 2
			if (cov >= 0.625 and cov < 0.875):
#			if (cov >= 0.6 and cov < 0.8):
				exp_cn = 3
			if (cov >= 0.875 and cov <= 1):
#			if (cov >= 0.8 and cov <= 1):
				exp_cn = 4
#			if (cov > 0 and cov < 0.33):
#				exp_cn = 0
#			if (cov >= 0.33 and cov < 0.66):
#				exp_cn = 1
#			if (cov >= 0.66 and cov <= 1):
#				exp_cn = 2
		cn = c_tuple.count(c_tuple[i])
		if (exp_cn != cn):
			costs+= 1
	return(costs)

#def cov_costs_general(c_tuple, var, coverage):
#	costs = 0
#	exp_cns = [i for i in range(ploidy+1)]
#	for i in range(ploidy):
#		cov = coverage[c_tuple[i]][var]
#		if (cov == 0):
#			return(1000000)
#		else: 
#			for j in exp_cns:
#				if (cov > j/(ploidy+1) and cov < (j+1)/(ploidy+1)):
#					exp_cn = j
#		cn = c_tuple.count(c_tuple[i])
#		if (cn != exp_cn):
#			costs +=1
#		return(costs)
			

def compute_index(tup, ploidy, num_clusters):
	index = 0
	for i in range(ploidy-1,-1,-1):
		index += tup[(ploidy-1)-i]*(num_clusters**i)
	return(index)
	
def clustering_DP(num_vars,clustering,coverage,positions, cov_map, ploidy, consensus, geno_map):
	scoring_matrix = subsetting(num_vars, clustering, coverage,positions, cov_map, ploidy, consensus, geno_map)
	return(scoring_matrix)