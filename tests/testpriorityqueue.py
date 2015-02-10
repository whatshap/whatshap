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

def test_queue2():
	pq = PriorityQueue()
	pq.push(1, 'a')
	pq.push(50, 'b')
	pq.push(2, 'c')
	pq.push(30, 'd')
	assert len(pq) == 4
	assert pq.pop() == (50, 'b')
	assert pq.pop() == (30, 'd')
	assert pq.pop() == (2, 'c')
	assert pq.pop() == (1, 'a')


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

def test_change_score_sorting():
	pq=PriorityQueue()
	pq.push(50,'1')
	pq.push(40,'2')
	pq.push(30,'3')
	pq.push(20,'4')
	pq.push(10,'5')
	pq.change_score('5',100)
	pq.change_score('2',45)
	assert len(pq) == 5
	assert pq.pop() == (100,'5')
	#assert pq.pop() == (45,'2')
	pq.push(60,'8')
	assert pq.pop() == (60,'8')
	pq.change_score('2',40)
	assert pq.pop()== (50,'1')
	assert pq.pop() == (40,'2')


def test_read_construction():
	pq=PriorityQueue()
	pq.push(2,(1,(0,1)))
	pq.push(2,(2,(0,1)))
	pq.push(4,(5,(0,1,2,3)))
	pq.push(3,(6,(0,1,2)))
	pq.push(4,(7,(0,1,2,3)))
	pq.push(2,(3,(1,2)))
	pq.push(3,(4,(1,2,3)))
	assert len(pq) == 7
	pq.pop()
	pq.pop()
	assert len(pq) == 5
	pq.change_score((4,(1,2,3)),1)
	assert pq.pop()== (3,(6,(0,1,2)))




'''
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
'''
