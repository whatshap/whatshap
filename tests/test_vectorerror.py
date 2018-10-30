"""
Test VectorError
"""

from whatshap.vectorerror import vector_error, num_switches

def test_num_switches():
	assert num_switches([0,1,2,3,4], [0,1,2,3,4]) == 0
	assert num_switches([0,1,3,2,4], [0,1,2,3,4]) == 2
	assert num_switches([0,4,3,2,1], [0,1,2,3,4]) == 4
	assert num_switches([0,4,3,2,1], [0,1,2,3,4]) == num_switches([0,1,2,3,4], [0,4,3,2,1])
	assert num_switches([2,3,1,4], [2,3,4,1]) == 2

def test_vector_error():
	phasing = [[0,0,0,0,0,0,0,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [0,0,0,0,1,1,1,1]]
	assert vector_error(phasing, truth) == 0
	
	phasing = [[0,0,0,0,0,0,0,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,1,1,1,1], [0,0,0,0,0,0,0,0]]
	assert vector_error(phasing, truth) == 0
	
	phasing = [[0,0,0,0,0,0,0,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0]]
	assert vector_error(phasing, truth) == 4
	
	phasing = [[1,1,1,1,0,0,0,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth) == 2
	
	phasing = [[1,1,1,1,0,0,1,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth) == 3
	
	phasing = [[1,1,1,1,0,0,1,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth, 5, 1) == 7
	
	phasing = [[1,1,1,1,0,0,1,0], [0,0,0,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth, 1, 10) == 7
	
	phasing = [[0,0,0,1,0,0,0,0], [1,1,1,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth, 1, 1) == 2
	
	phasing = [[0,0,0,1,0,0,0,0], [1,1,1,0,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth, 5, 1) == 4
	
	phasing = [[0,0,0,1,0,0,0,0], [1,1,1,1,1,1,1,1]]
	truth =   [[0,0,0,0,0,0,0,0], [1,1,1,1,1,1,1,1]]
	assert vector_error(phasing, truth, float("inf"), 1) == float("inf")