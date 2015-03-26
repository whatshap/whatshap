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

def test_privat_methods():
	pq=PriorityQueue()
	pq.push(10,'A')
	pq.push(9,'B')
	pq.push(8,'C')
	pq.push(7,'D')
	pq.push(6,'E')
	pq.push(5,'F')
	pq.push(4,'G')
	assert len(pq) == 7
	pq._swap(0,6)
	assert pq.pop()== (4,'G')
	pq._sift_up(5)
	assert pq.pop()==(10,'A')
	pq.push(50,'H')
	assert pq._get_score(0) == 50
	assert pq._get_item(0) == 'H'
	assert pq._get_index('H') == 0
	assert pq.get_score_by_item('D')== 7
	pq.pop()
	pq.pop()
	pq.pop()
	pq.pop()
	pq.pop()
	pq.pop()
	assert pq.is_empty()

def test_is_empty():
	pq=PriorityQueue()
	assert pq.is_empty()
	pq.push(10,'A')
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()
	pq.push(9,'B')
	assert not pq.is_empty()
	pq.push(8,'C')
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()
	pq.push(7,'D')
	assert not pq.is_empty()
	pq.push(6,'E')
	assert not pq.is_empty()
	pq.push(5,'F')
	assert not pq.is_empty()
	pq.push(4,'G')
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert not pq.is_empty()
	pq.pop()
	assert pq.is_empty()


def test_tuple_score():
	pq= PriorityQueue()
	assert pq.is_empty()
	pq.push((4,0,2),'A')
	(score,item ) = pq.pop()
	assert score== (4,0,2)
	assert item=='A'
	pq.push((1,0,2),'B')
	pq.push((2,0,2),'C')
	pq.push((3,0,2),'D')
	pq.push((40,0,2),'E')
	pq.push((0,0,2),'F')
	pq.push((50,0,2),'G')
	(score,item )=pq.pop()
	assert score == (50,0,2)
	assert item =='G'
	(score,item )=pq.pop()
	assert score == (40,0,2)
	assert item =='E'
	(score,item )=pq.pop()
	assert score == (3,0,2)
	assert item =='D'
	(score,item )=pq.pop()
	assert score == (2,0,2)
	assert item =='C'
	pq.pop()
	pq.pop()
	assert pq.is_empty()

def test_tuple_score_sorting():
	pq= PriorityQueue()
	pq.push((10,0,0),'B')
	pq.push((10,2,6),'C')
	pq.push((10,3,2),'D')
	pq.push((10,4,3),'E')
	pq.push((10,2,2),'F')
	pq.push((10,0,2),'G')
	(score,item ) = pq.pop()
	assert score== (10,4,3)
	assert item=='E'
	(score,item ) = pq.pop()
	assert score== (10,3,2)
	assert item=='D'
	(score,item ) = pq.pop()
	assert score== (10,2,6)
	assert item=='C'
	(score,item ) = pq.pop()
	assert score== (10,2,2)
	assert item=='F'
	(score,item ) = pq.pop()
	assert score== (10,0,2)
	assert item=='G'
	pq.push((1,10,4),'X')
	pq.push((5,0,6),'Y')
	pq.push((1,8,2),'Z')
	pq.change_score('Y',(100,100,100))
	pq.change_score('Z',(0,0,0))
	(score,item ) = pq.pop()
	assert score== (100,100,100)
	assert item=='Y'
	(score,item ) = pq.pop()
	assert score== (10,0,0)
	assert item=='B'
	(score,item ) = pq.pop()
	assert score== (1,10,4)
	assert item=='X'
	(score,item ) = pq.pop()
	assert score== (0,0,0)
	assert item=='Z'
	assert pq.is_empty()