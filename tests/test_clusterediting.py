from whatshap.core import ReadSet, LightCompleteGraph, CoreAlgorithm
from whatshap.testhelpers import string_to_readset, brute_force_phase
from whatshap.readscoring import score
import itertools

def test_clusterediting1():

	reads = """
		110000010111
		1100000101  
		 1000 01    
		 00 0 0 000 
		 1000001 11 
		  1111101   
		  0 10010 1 
		   0000 010 
		   1110     
		   0000 011 
		    000  00 
		    0001011 
		    0  10110
		    00010111
		    000 0000
		"""

	# construct a ReadSet
	readset = string_to_readset(reads)

	# compute similarities
	similarities = score(readset, 4, 0.1, 5)

	# create read graph
	n_reads = len(readset)
	graph = LightCompleteGraph(n_reads, True)
	
	# insert edges
	for id1 in range(n_reads):
		for id2 in range(id1+1, n_reads):
			graph.setWeight(id1, id2, similarities[id1][id2 - id1 - 1])

	# run cluster editing
	clusterediting = CoreAlgorithm(graph)	
	readpartitioning = clusterediting.run()

	print('computed clusters: ', readpartitioning)

	# make sure each read occurs only once
	read_ids = list(itertools.chain.from_iterable(readpartitioning))
	duplicates = set([ r for r in read_ids if read_ids.count(r) > 1 ])
	print('duplicates:', duplicates)
	assert(len(duplicates)  == 0)


