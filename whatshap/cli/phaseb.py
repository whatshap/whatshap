"""
Phase reads mapped to bubble chains

The output is a FASTA file written to standard output that needs to be piped to a file
"""

import logging
import os
import pickle
import sys
from copy import deepcopy
from whatshap.core import ReadSet, Read, Pedigree, PedigreeDPTable, NumericSampleIds, readselection
from whatshap.pedigree import uniform_recombination_map
from whatshap.vg_pb2 import Alignment
from whatshap.graph import ComponentFinder
import stream
import pdb

__author__ = "Fawaz Dabbaghie"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Node:
	def __init__(self, identifier):
		self.id = identifier
		self.seq = ""
		self.seq_len = 0
		self.start = []
		self.end = []
		self.visited = False
		self.which_chain = 0
		self.which_sb = 0
		self.which_b = 0
		self.which_allele = -1


class Bubble:
	def __init__(self, source, sink, alleles):
		self.source = source  # node id
		self.sink = sink  # node id
		self.alleles = alleles  # dictionary of two keys 0 and 1 and node ids


def merge_end(nodes, n, k, none_nodes):
	# print("in merge END with n {} and neighbor {}".format(n, nodes[n].end[0]))
	# print("and n's START are {}".format(nodes[n].start))

	if n != nodes[n].end[0][0]:
		neighbor = nodes[n].end[0]
		# checking if the neighbor is connected at start (so we have + + edge)
		# and that it only have one node from the start which is n
		if (neighbor[1] == 0) and (len(nodes[neighbor[0]].start) == 1):
			# the ends of n becomes the ends of neighbor
			# and the sequence and seq_len get updated
			nodes[n].end = deepcopy(nodes[neighbor[0]].end)
			# Here I need to check the new neighbors at end, and remove the merged node
			# and n to them
			for nn in nodes[neighbor[0]].end:
				# We are connected to it from start
				if nn[1] == 0:
					nodes[nn[0]].start.remove((neighbor[0], 1))
					nodes[nn[0]].start.append((n, 1))
				elif nn[1] == 1:
					# the if else here needed in case of there was a self loop on the end side of neighbor
					# the self loops is added to the merged node
					if nn[0] != neighbor[0]:
						nodes[nn[0]].end.remove((neighbor[0], 1))
						nodes[nn[0]].end.append((n, 1))
					else:
						nodes[n].end.remove((neighbor[0], 1))
						nodes[n].end.append((n, 1))
			nodes[n].seq += nodes[neighbor[0]].seq[k - 1:]
			nodes[n].seq_len = len(nodes[n].seq)
			nodes[neighbor[0]] = None
			none_nodes.append(neighbor[0])
			if len(nodes[n].end) == 1:
				merge_end(nodes, n, k, none_nodes)

		elif (neighbor[1] == 1) and (len(nodes[neighbor[0]].end) == 1):
			# the ends of n becomes the start of neighbor (because it's flipped)
			# and the sequence and seq_len get updated
			nodes[n].end = deepcopy(nodes[neighbor[0]].start)
			# Here I need to check the new neighbors at end of n, and remove the merged node
			# and add n to them
			for nn in nodes[neighbor[0]].start:
				# We are connected to it from start
				if nn[1] == 0:
					# the if else here needed in case of there was a self loop on the end side of neighbor
					# the self loops is added to the merged node
					if nn[0] != neighbor[0]:
						nodes[nn[0]].start.remove((neighbor[0], 0))
						nodes[nn[0]].start.append((nodes[n].id, 1))
					else:
						nodes[n].end.remove((neighbor[0], 0))
						nodes[n].end.append((n, 1))
				elif nn[1] == 1:
					nodes[nn[0]].end.remove((neighbor[0], 0))
					nodes[nn[0]].end.append((nodes[n].id, 1))

			reverse = reverse_complement(nodes[neighbor[0]].seq)

			nodes[n].seq += reverse[k - 1:]
			nodes[n].seq_len = len(nodes[n].seq)
			nodes[neighbor[0]] = None
			none_nodes.append(neighbor[0])

			if len(nodes[n].end) == 1:
				merge_end(nodes, n, k, none_nodes)


def merge_start(nodes, n, k, none_nodes):
	# print("in merge START with n {} and neighbor {}".format(n, nodes[n].start[0]))
	# print("and n's END are {}".format(nodes[n].end))
	if n != nodes[n].start[0][0]:  # no self loop
		neighbor = nodes[n].start[0]
		# checking if the neighbor is connected at start (so we have - + edge)
		# and that it only have one node from the start which is n
		if (neighbor[1] == 0) and (len(nodes[neighbor[0]].start) == 1):
			# the start of n becomes the ends of neighbor
			# and the sequence and seq_len get updated
			nodes[n].start = deepcopy(nodes[neighbor[0]].end)
			# Here I need to check the new neighbors at end, and remove the merged node
			# and n to them
			for nn in nodes[neighbor[0]].end:
				# We are connected to it from start
				if nn[1] == 0:
					nodes[nn[0]].start.remove((neighbor[0], 1))
					nodes[nn[0]].start.append((nodes[n].id, 0))

				elif nn[1] == 1:
					if nn[0] != neighbor[0]:
						nodes[nn[0]].end.remove((neighbor[0], 1))
						nodes[nn[0]].end.append((nodes[n].id, 0))
					else:
						nodes[n].start.remove((neighbor[0], 1))
						nodes[n].start.append((n, 0))

			reverse = reverse_complement(nodes[neighbor[0]].seq)

			nodes[n].seq = reverse[:len(reverse) - (k - 1)] + nodes[n].seq
			nodes[n].seq_len = len(nodes[n].seq)
			nodes[neighbor[0]] = None
			none_nodes.append(neighbor[0])

			if len(nodes[n].start) == 1:
				merge_start(nodes, n, k, none_nodes)

		elif (neighbor[1] == 1) and (len(nodes[neighbor[0]].end) == 1):
			merge_end(nodes, neighbor[0], k, none_nodes)


def compact_graph(nodes, k):
	# keeping the nodes that got merged to remove later
	none_nodes = []
	for n in nodes.keys():

		# checking that it's not a node that already got merged in the while loop
		if nodes[n] is not None:
			# checking if it has one neighbor and it's not a self loop
			if len(nodes[n].end) == 1:
				merge_end(nodes, n, k, none_nodes)
			if len(nodes[n].start) == 1:
				merge_start(nodes, n, k, none_nodes)

	# removing merged nodes
	for n in none_nodes:
		del nodes[n]

	return nodes


def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	return ''.join([complement[base] for base in seq[::-1]])


def add_arguments(parser):
	add = parser.add_argument
	add('gfa_file', metavar='GFA', type=str, help='Give the modified GFA file path')
	add('k', type=int, help='Give the value of K in the graph to remove the overlaps in the output')
	add('gam_file', metavar='GAM', type=str, help='Give the alignment GAM file path')


def validate(args, parser):

	if not os.path.exists(args.gfa_file):
		os.remove("phaseb.log")
		parser.error("The GFA file was not found")

	if not os.path.exists(args.gam_file):
		os.remove("phaseb.log")
		parser.error("The GAM file was not found")


def read_gfa(gfa_file_path, modified=False):
	"""
	:param gfa_file_path: gfa graph file.
	:param modified: if I'm reading my modified GFA with extra information for the nodes
	:return: Dictionary of node ids and Node objects.
	"""
	if not os.path.exists(gfa_file_path):
		print("the gfa file path you gave does not exists, please try again!")
		sys.exit()

	bubble_chains = dict()
	nodes = dict()
	edges = []
	with open(gfa_file_path, "r") as lines:
		for line in lines:
			if line.startswith("S"):
				if modified:
					line = line.split("\t")
					n_id = int(line[1])
					nodes[n_id] = Node(n_id)
					nodes[n_id].seq_len = len(line[2])
					nodes[n_id].seq = str(line[2])
					# this is my extra columns
					# constructed as which_chain:which_sb:which_b:which_allele
					# e.g. 5:0:3:1 (chain 5, bubble 3, allele 1)
					specifications = str(line[3])
					specifications = specifications.split(":")
					nodes[n_id].which_chain = int(specifications[0])
					nodes[n_id].which_sb = int(specifications[1])
					nodes[n_id].which_b = int(specifications[2])
					nodes[n_id].which_allele = int(specifications[3])

					if nodes[n_id].which_chain not in bubble_chains:
						bubble_chains[nodes[n_id].which_chain] = {nodes[n_id].which_b: [n_id]}
					else:
						if nodes[n_id].which_b not in bubble_chains[nodes[n_id].which_chain]:
							bubble_chains[nodes[n_id].which_chain][nodes[n_id].which_b] = [n_id]
						else:
							bubble_chains[nodes[n_id].which_chain][nodes[n_id].which_b].append(n_id)

				else:

					line = line.split()
					n_id = int(line[1])
					n_len = len(line[2])
					nodes[n_id] = Node(n_id)
					nodes[n_id].seq_len = n_len
					nodes[n_id].seq = str(line[2])

			elif line.startswith("L"):
				edges.append(line)

	for e in edges:
		line = e.split()

		k = int(line[1])
		neighbor = int(line[3])
		if (k not in nodes) or (neighbor not in nodes):
			continue

		if line[2] == "-":
			from_start = True
		else:
			from_start = False

		if line[4] == "-":
			to_end = True
		else:
			to_end = False

		if from_start is True and to_end is True:  # from start to end L x - y -
			if (neighbor, 1) not in nodes[k].start:
				nodes[k].start.append((neighbor, 1))
			if (k, 0) not in nodes[neighbor].end:
				nodes[neighbor].end.append((k, 0))

		elif from_start is True and to_end is False:  # from start to start L x - y +

			if (neighbor, 0) not in nodes[k].start:
				nodes[k].start.append((neighbor, 0))

			if (k, 0) not in nodes[neighbor].start:
				nodes[neighbor].start.append((k, 0))

		elif from_start is False and to_end is False:  # from end to start L x + y +
			if (neighbor, 0) not in nodes[k].end:
				nodes[k].end.append((neighbor, 0))

			if (k, 1) not in nodes[neighbor].start:
				nodes[neighbor].start.append((k, 1))

		elif from_start is False and to_end is True:  # from end to end L x + y -
			if (neighbor, 1) not in nodes[k].end:
				nodes[k].end.append((neighbor, 1))

			if (k, 1) not in nodes[neighbor].end:
				nodes[neighbor].end.append((k, 1))

	return nodes, bubble_chains


def build_readsets(nodes, gam_file_path):
	"""
	:param nodes: dictionary of nodes objects
	:param gam_file_path: gam file path
	:return: returns a dictionary of readsets, where each read set object represent a bubble chain
	"""
	all_readsets = {}

	# problem_read = "m140730_040106_42161_c100690211270000001823141603241592_s1_p0/33844/13172_22572"
	# debugging_nodes = [225, 1982, 2122, 2145, 2162, 2312, 2400, 2507, 2592, 3053, 3484, 3702, 3975,
	# 				   4066, 4215, 4533, 4728, 4924, 5280, 5360, 5552, 5606, 2018, 2325, 2482, 3042, 3626, 4335]
	#
	# debugging_dict = {}
	# for x in debugging_nodes:
	# 	debugging_dict[x] = 0

	with stream.open(str(gam_file_path), "rb") as instream:
		counter = 0
		for data in instream:
			counter += 1
			if (counter % 1000000) == 0:
				logger.info("Processed {} reads and we have {} readsets".format(counter, len(all_readsets)))

			g = Alignment()
			g.ParseFromString(data)

			# construct read object
			read = Read(g.name)  # leave other values as default

			for m in g.path.mapping:

				n_id = int(m.position.node_id)

				if n_id not in nodes:  # in case one of the nodes wasn't in the bubble chain
					continue

				# checking of there's a readset for the bubble chain
				if nodes[n_id].which_chain not in all_readsets:
					all_readsets[nodes[n_id].which_chain] = ReadSet()

				# the which allele value is either 1 or 0, otherwise the node doesn't belong to a bubble
				if (int(nodes[n_id].which_allele) == 0) or (int(nodes[n_id].which_allele) == 1):
					# if g.name == "m141230_065853_42225_c100724012550000001823140404301562_s1_p0/71345/0_14242":
					# 	n = nodes[n_id]
					# 	print("node {} adding variant {} in allele {} in chain {}".format(n_id, n.which_b,
					# 																	  n.which_allele, n.which_chain))

					read.add_variant(nodes[n_id].which_b, int(nodes[n_id].which_allele), 30)

			# only reads that cover 2 bubbles or more
			if len(read) >= 2:
				# some reads are present more than once and map to the same chain
				# maybe due to some big gaps that didn't map
				if n_id not in nodes:
					continue
				try:
					all_readsets[nodes[n_id].which_chain].add(read)
				except RuntimeError:
					if (int(nodes[n_id].which_allele) == 0) or (int(nodes[n_id].which_allele) == 1):
						read.add_variant(nodes[n_id].which_b, nodes[n_id].which_allele, 30)
				except KeyError:
					all_readsets[nodes[n_id].which_chain] = ReadSet()
					if (int(nodes[n_id].which_allele) == 0) or (int(nodes[n_id].which_allele) == 1):
						read.add_variant(nodes[n_id].which_b, nodes[n_id].which_allele, 30)

	return all_readsets


def build_readsets_from_dict(nodes, chains):
	"""
	:param nodes: dictionary of nodes objects
	:param chains: dictionary of read sets, inside dictionary of reads and variants
	:return: returns a dictionary of readsets, where each read set object represent a bubble chain
	"""
	all_readsets = {}

	counter = 0
	nodes_not_in_graph = 0

	for chain_n, chain in chains.items():
		if (counter % 1000000) == 0:
			logger.info("Processed {} reads and we have {} readsets".format(counter, len(all_readsets)))

		all_readsets[chain_n] = ReadSet()
		for read_name, read in chain.items():
			counter += 1
			read_obj = Read(read_name)
			# variants_to_add = list(set(read))
			# if len(variants_to_add) < 2:
			# 	continue
			for node in read:
				if node in nodes:
					read_obj.add_variant(nodes[node].which_b, int(nodes[node].which_allele), 30)
				else:
					nodes_not_in_graph += 1

			all_readsets[chain_n].add(read_obj)

	logger.info("Nodes in reads not in the graph are {}".format(nodes_not_in_graph))
	return all_readsets


def build_readsets_from_dict_debug_v(nodes, chains, list_of_chains):
	"""
	:param nodes: dictionary of nodes objects
	:param chains: dictionary of read sets, inside dictionary of reads and variants
	:param list_of_chains: for debugging
	:return: returns a dictionary of readsets, where each read set object represent a bubble chain
	"""
	all_readsets = {}

	counter = 0
	nodes_not_in_graph = 0

	for chain_n in list_of_chains:
		chain = chains[chain_n]
		if (counter % 1000000) == 0:
			logger.info("Processed {} reads and we have {} readsets".format(counter, len(all_readsets)))

		all_readsets[chain_n] = ReadSet()
		for read_name, read in chain.items():
			counter += 1
			read_obj = Read(read_name)
			# variants_to_add = list(set(read))
			# if len(variants_to_add) < 2:
			# 	continue
			for node in read:
				if node in nodes:
					read_obj.add_variant(nodes[node].which_b, int(nodes[node].which_allele), 30)
				else:
					nodes_not_in_graph += 1

			all_readsets[chain_n].add(read_obj)

	logger.info("Nodes in reads not in the graph are {}".format(nodes_not_in_graph))
	return all_readsets


def return_seq(nodes, haplotig, k, chain_n, bubble_haplotigs):
	"""
	:param nodes: dictionary of node objects
	:param haplotig: list of nodes in the haplotig to make a subgraph from and run
	:param k: overlap length
	:param chain_n: the chain key
	:param bubble_haplotigs: dict of bubble_haplotigs
	:return:
	"""
	path = dict()
	for n in haplotig:
		path[n] = Node(n)
		path[n].seq = nodes[n].seq
		for neighbor in nodes[n].start:
			if neighbor[0] in haplotig:
				path[n].start.append(neighbor)
		for neighbor in nodes[n].end:
			if neighbor[0] in haplotig:
				path[n].end.append(neighbor)

	# pdb.set_trace()
	path = compact_graph(path, k)
	try:
		assert len(path) == 1
		for value in path.values():
			return [value.seq]

	except AssertionError:
		logger.info("The path has length {} and its nodes are {} error happened on chain {} "
					"and the bubble haplotigs are {} "
					"and the haplotig is {}".format(len(path), list(path.keys()), chain_n, bubble_haplotigs, haplotig))
		seqs = []
		for value in path.values():
			seqs.append(value.seq)

		return seqs


def output_fasta(nodes, bubble_chains, bubble_membership, k, split_components, bubble_haplotigs, overall_components,
				 chain_n):
	"""
	:param nodes: dictionary of node objects
	:param bubble_chains: dict of bubble chains, each is a dict of bubble ids and nodes belonging to that bubble.
	:param k: overlap between nodes
	:param split_components: dicionary of readsets, and inside each one is split based on block
	:param bubble_haplotigs: dictionary mapping each bubble to the haplotigs
	:param overall_components: from whatshap
	:param chain_n: for debugging
	"""
	# I had to add this dictionary because some read set can bridge between two chains
	# I only know the original chain, but not the one it bridges too, but I know which bubbles it covered when
	# it bridged, so I can use the bubble_id to know which chain
	# I can probably change this to have bubble objects in a dictionary, then I don't care about the bubble_chain dict

	for block_n, block in split_components.items():
		for sub_block_n, sub_block in block.items():
			haplotig_1 = []
			haplotig_2 = []
			for bubble in sub_block:
				for node in bubble_chains[bubble_membership[bubble]][bubble]:
					try:
						if nodes[node].which_allele == -1:
							haplotig_1.append(node)
							haplotig_2.append(node)

						elif nodes[node].which_allele == bubble_haplotigs[bubble][0]:
							haplotig_1.append(node)
						elif nodes[node].which_allele == bubble_haplotigs[bubble][1]:
							haplotig_2.append(node)
					except KeyError:
						logger.info("Key error for bubble {}, this is weird".format(bubble))

			try:
				assert len(haplotig_1) == len(haplotig_2)
			except AssertionError:
				logger.info("haplotig_1 and 2 are not the same size for {}".format(split_components))
				continue

			seq_1 = return_seq(nodes, haplotig_1, k, chain_n, bubble_haplotigs)
			seq_2 = return_seq(nodes, haplotig_2, k, chain_n, bubble_haplotigs)
			try:
				assert len(seq_1) == len(seq_2)
			except AssertionError:
				logger.info("seq_1 and seq_2 were not equal in {}".format(overall_components))
				continue

			for idx, seq in enumerate(seq_1):
				if seq is not None:
					print(">" + "chain:" + str(chain_n) + "|block" + str(block_n) +
						  "|sub_block" + str(sub_block_n) + "_" + str(idx) + "|haplotig1")
					print(seq)
			for idx, seq in enumerate(seq_2):
				# todo implement the reverse thing better. However, the aligner should handle both version
				if seq[0:10] == reverse_complement(seq_1[idx])[0:10]:

					seq = reverse_complement(seq)
				if seq is not None:
					print(">" + "chain:" + str(chain_n) + "|block" + str(block_n) +
						  "|sub_block" + str(sub_block_n) + "_" + str(idx) + "|haplotig2")
					print(seq)


def split_blocks(overall_components):
	"""
	:param overall_components: the dictionary output from find_components
	:return: components: dictionary of blocks
	"""
	components = {}
	for pos, comp in overall_components.items():
		if comp not in components:
			components[comp] = [pos]
		else:
			components[comp].append(pos)

	# components now is a dictionary according to blocks
	# superbubbles can break the blocks further, so we need to account for that and have sub-blocks
	for comp_idx, comp in components.items():
		comp.sort()
		breaking_point = []
		new_comp = {}
		# pdb.set_trace()
		for i in range(1, len(comp)):
			# if there's a superbubble then the numbers for bubbles won't be consecutive
			if comp[i] != comp[i - 1] + 1:
				breaking_point.append(comp[i - 1] + 1)

		if breaking_point:
			for idx, i in enumerate(breaking_point):
				if idx == 0:
					block = [x for x in comp if x < i]
					if block:
						new_comp[len(new_comp) + 1] = block
				if 0 < idx < len(breaking_point) - 1:
					block = [x for x in comp if breaking_point[idx - 1] < x < i]
					if block:
						new_comp[len(new_comp) + 1] = [x for x in comp if breaking_point[idx - 1] < x < i]
				if idx == len(breaking_point) - 1:
					block1 = [x for x in comp if breaking_point[idx - 1] < x < i]
					block2 = [x for x in comp if x > i]
					if block1:
						new_comp[len(new_comp) + 1] = block1
					if block2:
						new_comp[len(new_comp) + 1] = block2

				components[comp_idx] = new_comp
		else:
			components[comp_idx] = {1: comp}

	# overall_components_keys = list(overall_components.keys())
	# overall_components_keys.sort()
	# split_components = []
	# for block_n, block in components.items():
	# 	for sub_block_n, sub_block in block.items():
	# 		split_components += sub_block
	# try:
	# 	print(overall_components_keys)
	# 	print(split_components)
	# 	assert len(overall_components_keys) == len(split_components)
	# except AssertionError:
	# 	pdb.set_trace()

	return components


def find_components(phased_positions, reads, master_block=None, heterozygous_positions=None):
	"""
	:param phased_positions: variatnts positions from the readset
	:param reads: a readset object
	:param master_block: Sure
	:param heterozygous_positions:
	:return: components dictionary

	Return a dict that maps each variant position to the component it is in.
	Variants are considered to be in the same component if a read exists that
	covers both. A component is identified by the position of its leftmost
	variant.
	master_block -- List of positions in a "master block", i.e. all blocks containing
					any of these positions are merged into one block.
	heterozygous_positions -- A dictionary mapping numeric sample ids to sets of
							  positions. Component building is then restricted to variants
							  at these positions. If none, all variants are used.
	"""
	# logger.debug('Finding connected components ...')
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		if heterozygous_positions is None:
			positions = [variant.position for variant in read if variant.position in phased_positions]
		else:
			positions = [variant.position for variant in read
						 if (variant.position in phased_positions) and (
								 variant.position in heterozygous_positions[read.sample_id])
						 ]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	if master_block is not None:
		for position in master_block[1:]:
			component_finder.merge(master_block[0], position)
	components = {position: component_finder.find(position) for position in phased_positions}
	return components


def read_pickled_dict(file_name):
	with open(file_name, 'rb') as handle:
		reads = pickle.load(handle)
	return reads


def run_phaseb(gfa_file, k, gam_file):
	"""
	:param gfa_file: path to the graph file
	:param k: value of overlaps between nodes in De Bruijn graphs
	:param gam_file: path to gam file
	:return:
	"""

#	debugging_chains = [48765  # 4 bubbles, had the same SE chr10:66669806
#		, 214543  # 4 bubbles, different SE, chr13:111364901
#		, 103916  # 11 bubbles, different SE chr13:111108634
#					   ]

	logger.info("reading graph file now...")
	nodes, bubble_chains = read_gfa(gfa_file, modified=True)
	logger.info("finished reading graph file.")
	bubble_membership = dict()
	for chain_k, chain in bubble_chains.items():
		for bubble in list(chain.keys()):
			bubble_membership[bubble] = chain_k

	logger.info("reading alignments and building readsets")
	if gam_file.endswith("pickle"):
		chains = read_pickled_dict(gam_file)

		# all_readsets = build_readsets_from_dict(nodes, chains)
		all_readsets = build_readsets_from_dict_debug_v(nodes, chains, debugging_chains)
	elif gam_file.endswith("gam"):
		all_readsets = build_readsets(nodes, gam_file)
	else:
		logger.info("It's neither a gam file nor a pickled dictionary")
		sys.exit(1)

	logger.info("finished reading alignments and building readset, there are {} readsets".format(len(all_readsets)))

	counter = 0

#	for chain_n in debugging_chains:
#		readset = all_readsets[chain_n]
#		pdb.set_trace()
	for chain_n, readset in all_readsets.items():

		logger.info("In chain {} with len {}".format(chain_n, len(bubble_chains[chain_n])))
		counter += 1
		if (counter % 10000) == 0:
			logger.info("processed {} readsets".format(counter))
		if len(readset) == 0:
			continue
		# if chain_n != 1:
		# 	continue
		# sort by starting positions of the reads
		all_read_pos = set()
		problematic_reads = set()
		for idx, read in enumerate(readset):
			all_read_pos.add(idx)
			try:
				read.sort()
			except RuntimeError:
				# this is happening when the same snp is mapped on both sides
				# because of loop I think, need th through these reads out I guess
				problematic_reads.add(idx)

		if len(problematic_reads) > 0:
			try:
				logger.info("Chain number {} with length {} "
							"had {} problematic reads".format(chain_n, len(bubble_chains[chain_n]), len(problematic_reads)))
			except KeyError:
				pdb.set_trace()

		# removing problematic reads
		readset = readset.subset(all_read_pos - problematic_reads)

		readset.sort()
		# print(len(readset))
		selected_indices = readselection(readset, 15)
		selected_reads = readset.subset(selected_indices)

		positions = selected_reads.get_positions()
		# construct pedigree
		numeric_sample_ids = NumericSampleIds()
		# print(numeric_sample_ids)

		recombination_cost = uniform_recombination_map(1.26, positions)
		pedigree = Pedigree(numeric_sample_ids)
		# print(recombination_cost)
		try:
			assert len(recombination_cost) == len(positions)
		except AssertionError:
			logger.info("For some reason recombination_cost didn't have the same length as positions")
			continue
			# pdb.set_trace()
		# assume all genotypes are heterozygous
		distrust_genotypes = False
		genotypes = [1] * len(positions)
		genotype_likelihoods = None
		sample_name = gfa_file.split(".")[0]
		pedigree.add_individual(sample_name, genotypes, genotype_likelihoods)
		# solve MEC
		# check positions
		dp_table = PedigreeDPTable(selected_reads, recombination_cost, pedigree, distrust_genotypes, positions)

		# cost = dp_table.get_optimal_cost()
		# vector indicating for each read in which bartition (0 or 1) it ended up (ordered according to ReadSet)
		read_partitioning = dp_table.get_optimal_partitioning()

		overall_components = find_components(positions, readset)
		if len(overall_components) < 2:  # in case there was only one bubble covered
			continue

		split_components = split_blocks(overall_components)
		# pdb.set_trace()

		# how to get the order of reads in readset:
		ordered_readnames = [read.name for read in selected_reads]
		try:
			assert len(read_partitioning) == len(ordered_readnames)
		except AssertionError:
			logger.info("For sosme reason read_partitioning and ordered_readnames are not equal")
			continue

		superreads_list, transmission_vector = dp_table.get_super_reads()
		two_haplotigs = superreads_list[0]
		# print("for the bubble_chain {} the two haplotigs are and the number of reads is {}".format(chain_n,len(readset)))
		haplotig = []
		# todo need to change this maybe to take both and make sure they are complementary
		for idx, read in enumerate(two_haplotigs):
			if idx == 0:  # just taking the first path as the second is just the opposite positions
				for x in list(read):
					haplotig.append(x.allele)

		# I need to make a dictionary of bubble_id haplotig
		try:
			assert len(haplotig) == len(overall_components.keys())
		except AssertionError:
			logger.info("haplotig and overall_components are not the same, something weird went wrong {}\t{}".format(haplotig, overall_components.keys()))

		bubble_haplotigs = {}
		for idx, b in enumerate(list(overall_components.keys())):
			if haplotig[idx] == 0:
				bubble_haplotigs[b] = [0, 1]
			elif haplotig[idx] == 1:
				bubble_haplotigs[b] = [1, 0]
			else:
				# this is in case I had equal scores, I'll change the seq in the branch nodes to N
				# and I randomly assign 0 and 1
				# I need the bubble_membership dict because the overall component can give back bubbles not in the
				# current branch (some reads can branch out to other chains), so I can check where does the bubble come
				# from, then check the nodes inside it
				bubble_haplotigs[b] = [0, 1]
				for node in bubble_chains[bubble_membership[b]][b]:
					if (nodes[node].which_allele == 1) or (nodes[node].which_allele == 0):
						nodes[node].seq = "N"*nodes[node].seq_len

		try:
			assert len(bubble_haplotigs) == len(overall_components)
		except AssertionError:
			logger.info("Split components are {} and haplotig vector is {} and "
						"the superreads are {}".format(split_components, haplotig, two_haplotigs))
			continue
		pdb.set_trace()
		if chain_n == 25:
			print("These are the selected reads")
			print(selected_reads)
			print("These are the bubble haplotigs")
			print(bubble_haplotigs)
			print("These are the superreads")
			print(two_haplotigs)
			print("This is the order of reads")
			for idx, r in enumerate(ordered_readnames):
				print("idx\t{} belongs to haplotype {}".format(str(r), read_partitioning[idx]))
			
		# print(split_components)
		# pdb.set_trace()
		# print(overall_components)
		# pdb.set_trace()
		# print("##################################################################")
		#
		# output_fasta(nodes, bubble_chains, bubble_membership, k, split_components, bubble_haplotigs, overall_components,
		# 			 chain_n)


def main(args):
	# here all the main stuff going to happen
	run_phaseb(**vars(args))
