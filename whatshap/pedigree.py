"""
Pedigree-related functions
"""
import math
from collections import namedtuple, Counter, OrderedDict

RecombinationMapEntry = namedtuple('RecombinationMapEntry', ['position', 'cum_distance'])


def interpolate(point, start_pos, end_pos, start_value, end_value):
	assert start_pos <= point <= end_pos
	if start_pos == point == end_pos:
		assert start_value == end_value
		return start_value
	return start_value + ((point - start_pos) * (end_value - start_value) / (end_pos - start_pos))


def recombination_cost_map(genetic_map, positions):

	assert len(genetic_map) > 0


	# Step 1: compute cumulative genetic distances from start of chromosome
	#         to each position.
	cumulative_distances = []
	# i and j are such that genetic_map[i].position <= position <= genetic_map[j].position
	# i and j are None if no such values exist (because we are at the end of the list)
	i = None
	j = 0

	for position in positions:
		# update i to meet the invariant
		if (i == None) and (genetic_map[0].position <= position):
			i = 0
		while (i != None) and (i+1 < len(genetic_map)) and (genetic_map[i+1].position <= position):
			i += 1

		# update j to meet the invariant
		while (j != None) and (genetic_map[j].position < position):
			if j+1 < len(genetic_map):
				j += 1
			else:
				j = None

		# interpolate
		if i == None:
			assert j != None
			d = interpolate(position, 0, genetic_map[j].position, 0, genetic_map[j].cum_distance)
		elif j == None:
			# Point outside the genetic map --> extrapolating using average recombination rate
			avg_rate = genetic_map[-1].cum_distance / genetic_map[-1].position
			d = genetic_map[-1].cum_distance + (position - genetic_map[-1].position) * avg_rate
		else:
			assert genetic_map[i].position <= position <= genetic_map[j].position
			d = interpolate(position, genetic_map[i].position, genetic_map[j].position, genetic_map[i].cum_distance, genetic_map[j].cum_distance)
		cumulative_distances.append(d)

	# Step 2: compute costs (= phred-scaled recombination probabilities between two positions)
	result = [0]
	for i in range(1, len(cumulative_distances)):
		d = cumulative_distances[i] - cumulative_distances[i-1]
		result.append(round(centimorgen_to_phred(d)))

	return result


def uniform_recombination_map(recombrate, positions):
	"""
	For a list of positions and a constant recombination rate (in cM/Mb),
	return a list "results" of the same length as "positions" such that
	results[i] is the phred-scaled recombination probability between
	positions[i-1] and positions[i].
	"""
	return [0] + [ round(centimorgen_to_phred((positions[i] - positions[i-1]) * 1e-6 * recombrate)) for i in range(1, len(positions)) ]


def centimorgen_to_phred(distance):
	assert distance >= 0
	if distance == 0:
		return large_value
	p = (1.0-math.exp(-(2.0*distance)/100))/2.0
	return -10 * math.log10(p)


def load_genetic_map(filename):
	genetic_map = []

	with open(filename,'r') as fid:

		# read and ignore first line
		fid.readline()

		# for each line only store the first and third value in two seperate list
		for line in fid:
			line_spl = line.strip().split()
			assert len(line_spl) == 3
			genetic_map.append(
				RecombinationMapEntry(position=int(line_spl[0]), cum_distance=float(line_spl[2]))
			)

	return genetic_map


mendelian_conflict_sets = set([
	(0, 0, 1),
	(0, 0, 2),
	(0, 1, 2),
	(0, 2, 0),
	(0, 2, 2),
	(1, 2, 0),
	(2, 2, 0),
	(2, 2, 1)
])


def mendelian_conflict(genotypem, genotypef, genotypec):
	l = [genotypem, genotypef]
	l.sort()
	l.append(genotypec)
	return tuple(l) in mendelian_conflict_sets


def find_recombination(transmission_vector, components, positions, recombcost, recombination_list_filename=None):
	assert len(transmission_vector) == len(positions) == len(recombcost)
	assert set(components.keys()).issubset(set(positions))
	position_to_index = { pos: i for i, pos in enumerate(positions) }
	blocks = defaultdict(list)
	for position, block_id in components.items():
		blocks[block_id].append(position)

	RecombinationEvent = namedtuple('RecombinationEvent', ['position1', 'position2', 'inheritance_value1', 'inheritance_value2' , 'recombination_cost'])
	event_count = 0
	event_list = []
	cum_recomb_cost = 0
	for block_id, block in blocks.items():
		block.sort()
		block_transmission_vector = [ transmission_vector[position_to_index[i]] for i in block ]
		block_recomb_cost = [ recombcost[position_to_index[i]] for i in block ]
		if len(block) <= 2:
			continue
		for i in range(2, len(block)):
			if block_transmission_vector[i-1] != block_transmission_vector[i]:
				event_list.append(RecombinationEvent(block[i-1], block[i], block_transmission_vector[i-1], block_transmission_vector[i], block_recomb_cost[i]))
				cum_recomb_cost += block_recomb_cost[i]
				event_count += 1

	if recombination_list_filename is not None:
		event_list.sort()
		f = open(recombination_list_filename, 'w')
		print('#position1', 'position2', 'inheritance_value1', 'inheritance_value2' , 'recombination_cost', file=f)
		for e in event_list:
			print(e.position1, e.position2, e.inheritance_value1, e.inheritance_value2, e.recombination_cost, file=f)
		f.close()

	logger.info('Cost accounted for by recombination events: %d', cum_recomb_cost)
	return event_count


class ParseError(Exception):
	pass


# TODO rename this to Trio with attributes mother, father, child
Individual = namedtuple('Individual', ['id', 'mother_id', 'father_id'])


class PedReader:
	"""
	A parser for PED/FAM files as used by PLINK and other tools.

	According to <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>:
	The PED file is a white-space (space or tab) delimited file: the first six
	columns are mandatory:
	* Family ID
	* Individual ID
	* Paternal ID
	* Maternal ID
	* Sex (1=male; 2=female; other=unknown)
	* Phenotype

	All fields except the individual, maternal and paternal ID are ignored by
	this class. The entire file is read upon construction.
	"""
	def __init__(self, file, numeric_sample_ids):
		if isinstance(file, str):
			with open(file) as f:
				self.individuals = self._parse(f)
		else:
			self.individuals = self._parse(file)
		# Ensure that all individuals mentioned in the ped have a numeric id
		for individual in self.individuals:
			numeric_sample_ids[individual.id]

	@staticmethod
	def _parse_record(line):
		"""
		Parse a single non-comment line of a PED or FAM file.
		"""
		fields = line.split()
		if len(fields) < 6:
			raise ParseError("Less than six fields found in PED/FAM file")
		individual_id, paternal_id, maternal_id = fields[1:4]
		if paternal_id == '0':
			paternal_id = None
		if maternal_id == '0':
			maternal_id = None
		return Individual(id=individual_id, mother_id=maternal_id, father_id=paternal_id)

	def _parse(self, file):
		individuals = []
		for line in file:
			if line.startswith('#') or line == '\n':
				continue
			individuals.append(self._parse_record(line))
		self._sanity_check(individuals)
		return individuals

	@staticmethod
	def _sanity_check(individuals):
		"""
		Ensure that each individual occurs only once in the file.
		"""
		individual_ids = [ individual.id for individual in individuals ]
		if not individual_ids:
			return
		id, count = Counter(individual_ids).most_common()[0]
		if count > 1:
			raise ParseError('Individual {!r} occurs more than once in PED file'.format(id))

	def __iter__(self):
		return iter(self.individuals)


class CyclicGraphError(Exception):
	pass


class Graph:
	"""Directed graph that can sort topologically"""
	def __init__(self):
		# map node to a list of neighbors
		self._neighbors = OrderedDict()

	def add_edge(self, node1, node2):
		"""The edge is directed from node1 to node2"""
		if node1 not in self._neighbors:
			self._neighbors[node1] = []
		self._neighbors[node1].append(node2)
		if node2 not in self._neighbors:
			self._neighbors[node2] = []

	def toposorted(self):
		"""
		Return nodes of the graph sorted topologically.
		For all edges u -> v that the graph has, node v will appear
		before node u.
		"""
		order = []
		colors = { node: 'white' for node in self._neighbors }
		def visit(node):
			assert colors[node] == 'white'
			colors[node] = 'gray'
			for neighbor in self._neighbors[node]:
				if colors[neighbor] == 'white':
					visit(neighbor)
				elif colors[neighbor] == 'gray':
					raise CyclicGraphError(
						'Cycle involving {!r} and {!r} detected'.format(
						node, neighbor))
			order.append(node)
			colors[node] = 'black'

		for node in self._neighbors:
			if colors[node] == 'white':
				visit(node)
		return order
