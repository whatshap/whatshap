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
