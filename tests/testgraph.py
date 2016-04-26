from nose.tools import raises
from whatshap.pedigree import Graph, CyclicGraphError


def test_graph():
	graph = Graph()
	tuples = [
		# mother, father, child
		('mmm', 'mmf', 'mm'),
		('mf', 'mm', 'm'),
		('m', 'f', 'c1'),
		('m', 'f', 'c2'),
		('ff', 'fm', 'f'),
	]
	for mother, father, child in tuples:
		graph.add_edge(child, mother)
		graph.add_edge(child, father)
	t = graph.toposorted()
	assert len(t) == len(set(t))
	# parents must be before children
	for mother, father, child in tuples:
		assert t.index(mother) < t.index(child)
		assert t.index(father) < t.index(child)


@raises(CyclicGraphError)
def test_cyclic():
	graph = Graph()
	tuples = [
		# mother, father, child
		('mmm', 'mmf', 'mm'),
		('mf', 'mm', 'm'),
		('m', 'f', 'c1'),
		('m', 'f', 'c2'),
		('ff', 'fm', 'f'),
		('c1', 'c2', 'mmf')  # this makes it cyclic
	]
	for mother, father, child in tuples:
		graph.add_edge(child, mother)
		graph.add_edge(child, father)
	graph.toposorted()
