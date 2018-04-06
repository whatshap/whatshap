from nose.tools import raises
from whatshap.pedigree import Graph, CyclicGraphError


def assert_toposort(tuples):
	graph = Graph()
	for mother, father, child in tuples:
		graph.add_edge(child, mother)
		graph.add_edge(child, father)
	t = graph.toposorted()
	assert len(t) == len(set(t))
	# parents must be before children
	for mother, father, child in tuples:
		assert t.index(mother) < t.index(child)
		assert t.index(father) < t.index(child)


def test_graph():
	tuples = [
		# mother, father, child
		('mmm', 'mmf', 'mm'),
		('mf', 'mm', 'm'),
		('m', 'f', 'c1'),
		('m', 'f', 'c2'),
		('ff', 'fm', 'f'),
	]
	assert_toposort(tuples)


def test_charles_ii():
	# source: https://en.wikipedia.org/wiki/Charles_II_of_Spain
	# fields are: child, father, mother
	pedigree = \
	"""
	Charles II of Spain, Philip IV of Spain, Mariana of Austria

	Mariana of Austria, Ferdinand III, Maria Anna of Spain

	Philip IV of Spain, Philip III of Spain, Margaret of Austria
	Maria Anna of Spain, Philip III of Spain, Margaret of Austria
	Ferdinand III, Ferdinand II, Maria Anna of Bavaria (1574-1616)

	Philip III of Spain, Philip II of Spain, Anna of Austria (1549-80)
	Margaret of Austria, Charles II Archduke of Austria, Maria Anna of Bavaria
	Ferdinand II, Charles II Archduke of Austria, Maria Anna of Bavaria
	Maria Anna of Bavaria (1574-1616), William V Duke of Bavaria, Renata of Lorraine

	Anna of Austria (1549-80), Maximilian II, Maria of Spain
	Maria Anna of Bavaria, Albert V Duke of Bavaria, Anna of Austria
	William V Duke of Bavaria, Albert V Duke of Bavaria, Anna of Austria
	Renata of Lorraine, Francis I Duke of Lorraine, Christina of Denmark

	Philip II of Spain, Charles V, Isabella of Portugal
	Maria of Spain, Charles V, Isabella of Portugal
	Maximilian II, Ferdinand I, Anna of Bohemia and Hungary
	Charles II Archduke of Austria, Ferdinand I, Anna of Bohemia and Hungary
	Anna of Austria, Ferdinand I, Anna of Bohemia and Hungary
	Christina of Denmark, Christian II of Denmark, Isabella of Austria
	Charles V, Philip I of Castile, Joanna of Castille
	Ferdinand I, Philip I of Castile, Joanna of Castille
	Anna of Bohemia and Hungary, Philip I of Castile, Joanna of Castille
	Isabella of Austria, Philip I of Castile, Joanna of Castille
	"""
	individuals = set()
	tuples = []
	for line in pedigree.split('\n'):
		line = line.strip()
		if line == '':
			continue
		child, father, mother = line.split(', ')
		assert child not in individuals, child
		tuples.append((mother, father, child))
	assert_toposort(tuples)


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
