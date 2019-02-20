from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
cimport cpp

cdef class StaticSparseGraph:
	def __cinit__(self, int numNodes):
		self.thisptr = new cpp.StaticSparseGraph(numNodes)
	def addEdge(self, int node_id1, int node_id2, double weight):
		self.thisptr.addEdge(node_id1, node_id2, weight)

cdef class CoreAlgorithm:
	def __cinit__(self, StaticSparseGraph graph, bundleEdges):
		self.thisptr = new cpp.CoreAlgorithm(graph.thisptr[0], bundleEdges)
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

	def scoreReadsetGlobal(self, ReadSet readset, unsigned int minOverlap, unsigned int ploidy, double errorrate):
		sim = TriangleSparseMatrix()
		self.thisptr.scoreReadsetGlobal(sim.thisptr, readset.thisptr, minOverlap, ploidy, errorrate)
		return sim
	
	def scoreReadsetLocal(self, ReadSet readset, unsigned int minOverlap, unsigned int ploidy):
		sim = TriangleSparseMatrix()
		self.thisptr.scoreReadsetLocal(sim.thisptr, readset.thisptr, minOverlap, ploidy)
		return sim
	
	def scoreReadsetPatterns(self, ReadSet readset, unsigned int minOverlap, unsigned int ploidy, double errorrate, unsigned int windowSize):
		sim = TriangleSparseMatrix()
		self.thisptr.scoreReadsetPatterns(sim.thisptr, readset.thisptr, minOverlap, ploidy, errorrate, windowSize)
		return sim