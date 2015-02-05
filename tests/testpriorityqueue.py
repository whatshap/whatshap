from whatshap.priorityqueue import PriorityQueue

def test_queue():
	pq = PriorityQueue()
	pq.push(10, 'a')
	pq.push(5, 'b')
	pq.push(12, 'c')
	pq.push(3, 'd')
	assert len(pq) == 4
	assert pq.pop() == (12, 'c')
	assert pq.pop() == (10, 'a')
	assert pq.pop() == (5, 'b')
	assert pq.pop() == (3, 'd')

def test_change_score():
	pq = PriorityQueue()
	pq.push(10, 'a')
	pq.push(5, 'b')
	pq.change_score('a', 2)
	pq.push(12, 'c')
	pq.push(3, 'd')
	pq.change_score('c', 1)
	pq.change_score('d', 15)
	assert len(pq) == 4
	assert pq.pop() == (15, 'd')
	assert pq.pop() == (5, 'b')
	assert pq.pop() == (2, 'a')
	assert pq.pop() == (1, 'c')


def test_construction():
	pq=PriorityQueue()

	#only for testing of the readstructure
	#testarray_for_heap_analysis = {}
	#testarray_for_heap_analysis[0] = [[1, 2, [0, 1]], [2, 2, [0, 1]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]],[7, 4, [0, 1, 2, 3]]]
	#testarray_for_heap_analysis[1] = [[1, 2, [0, 1]], [2, 2, [0, 1]], [3, 2, [1, 2]], [4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]], [7, 4, [0, 1, 2, 3]]]
	#testarray_for_heap_analysis[2] = [[3, 2, [1, 2]], [4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]], [7, 4, [0, 1, 2, 3]]]
	#testarray_for_heap_analysis[3] = [[4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [7, 4, [0, 1, 2, 3]]]
	#testarray_for_heap_analysis[4]= [[8,9,[4]]]
	#testarray_for_heap_analysis[5]=[]

def test_comparison():
	pq=PriorityQueue()

	pq.push(50,'a')
	pq.push(40,'b')


	assert pq.getscore(1)== 50