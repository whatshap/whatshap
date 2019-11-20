from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t
cimport cpp

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