"""
Pedigree-related functions
"""
import math
from collections import namedtuple, Counter, OrderedDict
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

RecombinationMapEntry = namedtuple('RecombinationMapEntry', ['position', 'cum_distance'])


def interpolate(point, start_pos, end_pos, start_value, end_value):
	assert start_pos <= point <= end_pos
	if start_pos == point == end_pos:
		assert start_value == end_value
		return start_value
	return start_value + ((point - start_pos) * (end_value - start_value) / (end_pos - start_pos))

MINIMUM_GENETIC_DISTANCE = 1e-10 # cM

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
		d = max(d, MINIMUM_GENETIC_DISTANCE)
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
		raise ValueError('Cannot convert genetic distance of zero to phred.')
	elif distance < 1e-10:
		return -10 * (math.log10(distance) - 2)
	else:
		p = (1.0-math.exp(-(2.0*distance)/100))/2.0
		return -10 * math.log10(p)


def load_genetic_map(filename):
	genetic_map = []

	warned_zero_distance = False
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
			if len(genetic_map) >= 2:
				if not warned_zero_distance and (genetic_map[-2].cum_distance == genetic_map[-1].cum_distance):
					logger.warning('Zero genetic distances encountered in %s', filename)
					warned_zero_distance = True

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


def find_recombination(transmission_vector, components, positions, recombcost):
	assert len(transmission_vector) == len(positions) == len(recombcost)
	assert set(components.keys()).issubset(set(positions))
	position_to_index = { pos: i for i, pos in enumerate(positions) }
	blocks = defaultdict(list)
	for position, block_id in components.items():
		blocks[block_id].append(position)

	RecombinationEvent = namedtuple('RecombinationEvent', ['position1', 'position2', 'transmitted_hap_mother1', 'transmitted_hap_mother2' ,'transmitted_hap_father1', 'transmitted_hap_father2', 'recombination_cost'])
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
				event_list.append(RecombinationEvent(block[i-1], block[i], block_transmission_vector[i-1]%2, block_transmission_vector[i]%2, block_transmission_vector[i-1]//2, block_transmission_vector[i]//2, block_recomb_cost[i]))
				cum_recomb_cost += block_recomb_cost[i]

	logger.info('Cost accounted for by recombination events: %d', cum_recomb_cost)
	event_list.sort()
	return event_list


class ParseError(Exception):
	pass


Trio = namedtuple('Trio', ['child', 'mother', 'father'])
Trio.__doc__ = """
Relationships are modelled as a set of trios (mother, father, child).
"""


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
				self.trios = self._parse(f)
		else:
			self.trios = self._parse(file)
		# Ensure that all mentioned individuals have a numeric id
		for trio in self.trios:
			numeric_sample_ids[trio.child]

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
		return Trio(child=individual_id, mother=maternal_id, father=paternal_id)

	def _parse(self, file):
		trios = []
		for line in file:
			if line.startswith('#') or line == '\n':
				continue
			trios.append(self._parse_record(line))
		self._sanity_check(trios)
		return trios

	@staticmethod
	def _sanity_check(trios):
		"""
		Ensure that each individual occurs only once in the file.
		"""
		children = [ trio.child for trio in trios ]
		if not children:
			return
		id, count = Counter(children).most_common()[0]
		if count > 1:
			raise ParseError('Individual {!r} occurs more than once in PED file'.format(id))

	def __iter__(self):
		return iter(self.trios)


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
