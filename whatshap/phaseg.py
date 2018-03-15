"""
create association between reads and bubbles.
"""
# Dependencies to be installed: pysam, pyfaidx,xopen,pyvcf,protobuf

import pyfaidx

from xopen import xopen
import stream
import logging
from . import vg_pb2
from collections import Counter
import sys
from collections import defaultdict
from .core import ReadSet, Read
from functools import reduce
import operator as op
from itertools import groupby

import random
import collections
from collections import OrderedDict, namedtuple


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
from .core import PyIndexSet as IndexSet
from .core import ReadSet, readselection, Pedigree, PedigreeDPTable, NumericSampleIds, PhredGenotypeLikelihoods
from .graph import ComponentFinder
from .pedigree import (PedReader, mendelian_conflict, recombination_cost_map,
                       load_genetic_map, uniform_recombination_map, find_recombination)
from .bam import BamIndexingError, SampleNotFoundError
from .timer import StageTimer
from .variants import ReadSetReader, ReadSetError
from heapq import heappush, heappop
from itertools import count


__author__ = "Shilpa Garg, Tobias Marschall"

logger = logging.getLogger(__name__)


class CoverageMonitor:
        '''TODO: This is a most simple, naive implementation. Could do this smarter.'''
        def __init__(self, length):
                self.coverage = [0] * length

        def max_coverage_in_range(self, begin, end):
                return max(self.coverage[begin:end])

        def add_read(self, begin, end):
                for i in range(begin, end):
                        self.coverage[i] += 1



def find_components(phased_positions, reads, master_block=None):
	"""
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
	logger.debug('Finding connected components ...')
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [ variant.position for variant in read if variant.position in phased_positions ]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	if not master_block is None:
		for position in master_block[1:]:
			component_finder.merge(master_block[0], position)
	components = { position : component_finder.find(position) for position in phased_positions }
	return components


def find_largest_component(components):
	"""
	Determine the largest component and return a sorted list of positions
	contained in it.
	components -- dictionary mapping positin to block_id as returned by find_components.
	"""
	blocks = defaultdict(list)
	for position, block_id in components.items():
		blocks[block_id].append(position)
	largest = []
	for block in blocks.values():
		if len(block) > len(largest):
			largest = block
	largest.sort()
	return largest
      

"""
output the possible allele-pairs for a bubble.
"""
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer//denom
  

"""
build node-sequence list for vg graph
"""
def vg_graph_reader(vg_file):
	node_seq_list= defaultdict()
	edge_connections = defaultdict(list)
	with stream.open(str(vg_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.Graph()
			l.ParseFromString(data)
			for i in range(len(l.node)):
				index = l.node[i].id
				seq = l.node[i].sequence
				node_seq_list[index]=seq
			for j in range(len(l.edge)):
				from_edge = getattr(l.edge[j], "from")
				edge_connections[from_edge].append(l.edge[j].to)
	return node_seq_list, edge_connections


"""
Output phased SnarlTraversal using superreads_list. 
Then take input vg graph and phased SnarlTraversal to reconstruct underlying two sequences.
"""
"""
Input: Phase variants from Locus file and aligned reads from GAM file.

It creates an association between het variants and read alignments. 

Output: The optimal partitioning is written to standard output.
"""
def vg_reader(locus_file, gam_file):
	"""
	input: sorted locus and sorted GAM file output from vg.
	output: sorted readset for core DP.
	assumptions: 
	1. locus file consists of linear ordering of simple bubbles only and hence sorted. Each locus file does not contain start and end vertex.
	2. paths in the locus should be covered by atleast one pacbio read.
	2. GAM file is sorted and restricted to locus file.
	3. files consists of all DAG connected components.
	4. add variant only when it identifies the branch uniquely.
	"""

	locus_count=0
	prev_startsnarl = 0
	prev_endsnarl = 0
	locus_branch_mapping=OrderedDict()
	locus_count=0
	prev_startsnarl = 0
	prev_startsnarl_orientation = -1
	prev_endsnarl = 0
	prev_endsnarl_orientation = -1
	start_end_bubblenods = set()
	with stream.open(str(locus_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.SnarlTraversal()
			l.ParseFromString(data)
			#TODO: make ordered doctionary locus_branch_mapping
			# handle forward and backward case of nodes
			current_startsnarl = l.snarl.start.node_id
			current_startsnarl_orientation = l.snarl.start.backward
			current_endsnarl = l.snarl.end.node_id
			current_endsnarl_orientation = l.snarl.end.backward
			path_in_bubble =[]
			start_end_bubblenods.add(l.snarl.end.node_id)
			start_end_bubblenods.add(l.snarl.start.node_id)
			if len(l.visits) ==0:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.snarl.end.node_id)))
			else:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id, l.visits[-1].node_id)))
					for i in range(0,len(l.visits)-1):
						path_in_bubble.append(tuple((l.visits[i+1].node_id, l.visits[i].node_id)))
					path_in_bubble.append(tuple ((l.visits[0].node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
					for i in range(0,len(l.visits)-1):
						path_in_bubble.append(tuple((l.visits[i].node_id, l.visits[i+1].node_id)))
					path_in_bubble.append(tuple ((l.visits[-1].node_id, l.snarl.end.node_id))) 

			if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
				per_locus.append(path_in_bubble)
			else:
				locus_count=locus_count+1
				per_locus = []
				per_locus.append(path_in_bubble)
			prev_startsnarl = current_startsnarl
			prev_startsnarl_orientation = current_startsnarl_orientation
			prev_endsnarl = current_endsnarl
			prev_endsnarl_orientation = current_endsnarl_orientation
			locus_branch_mapping[locus_count]=per_locus

	print('The number of hets:')
	het_count= 0
	for k,v in locus_branch_mapping.items():
		if len(v) >1:
			het_count = het_count +1
	print(het_count)
	# keep branch of paths in each bubble.
	alleles_per_pos= defaultdict()
	for k,v in locus_branch_mapping.items():
		alleles_per_pos[k]=len(v)

	# both simple and complex bubbles: key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles)
	reverse_mapping= defaultdict(list)
	for k,v in locus_branch_mapping.items():
		if len(v) > 1: # more than one branch
			for i,b in enumerate(v):
				if len(b) > 0:
					for p,j in enumerate(b):
						reverse_mapping[j].append([k,i, len(v)]) # in complex bubbles, a node can map to multiple branches.


	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	readset=ReadSet()
	count =0
	duplicated = 0
	#TODO: consider reads with only positive score.
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			# hard-coded source id, mapping quality and other values.
			val1 = True
			val2 = False

			count1 =0
			count2=0
			score = g.score/len(g.sequence)

			if score > 0.2:
				continue
			read=Read(g.name, 0, 0, 0) # create read for each read alignment

			prev_tmp=[]
			prev_locus= -1
			locus = -1

			for i in range(0,len(g.path.mapping)-1):
			#for i in g.path.mapping: # go over the mapping in a read
				# TODO: check for forward or reverse strand, we may not need it for DAG.
				edge1 = tuple((g.path.mapping[i].position.node_id, g.path.mapping[i+1].position.node_id)) # go over nodes in a mapping
				edge2 = tuple((g.path.mapping[i+1].position.node_id, g.path.mapping[i].position.node_id)) # go over nodes in a mapping
				if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
					if edge1 in reverse_mapping:
						qualities = [10]* reverse_mapping[edge1][0][2]
						node_inf= [tuple(k[0:2]) for k in reverse_mapping[edge1]] # consider (locus, branch)
					else:
						qualities = [10]* reverse_mapping[edge2][0][2]
						node_inf= [tuple(k[0:2]) for k in reverse_mapping[edge2]]
					tmp = [x for x in node_inf]
					if prev_locus != tmp[0][0]:
						prev_tmp = tmp
						prev_locus = tmp[0][0]
						
					interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
					if len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
						qualities[interset_tmp[0][1]] = 0
						if i== len(g.path.mapping)-1:
							read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
						else:
							if g.path.mapping[i+1].position.node_id in start_end_bubblenods:
								read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)    

						locus= interset_tmp[0][0]
						

			readset.add(read)

	readset1=ReadSet()
	tmp_duplicated=set()
	for read in readset:
		if read.sort() ==1:
			duplicated = duplicated +1
			tmp=[]
			for variant in read:
				tmp.append(variant.position)
			#print("duplicated variant")
			x = [item for item, count in collections.Counter(tmp).items() if count > 1]
			for a in x:
				tmp_duplicated.add(a)
			continue
		else:
			if len(read) >=2:
			      readset1.add(read)
	print("length of duplicated bubbles")
	print(tmp_duplicated)


	readset1.sort()
	
	print("reads considered before read-selection")
	print(len(readset1))
	return readset1, alleles_per_pos, locus_branch_mapping


def reverse_complement(seq):
	seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def generate_haplotigs(sample_superreads, components, node_seq_list, locus_branch_mapping, canu_alignments, vg_file, pred_haplotigs, locus_file):
	sample = 0
	pred_haplotigs_file = open(pred_haplotigs, 'w')

	# This holds a dict from (node ID, orientation) pair, where true is reverse
	# (leftward) and false is forward (rightward) to a set of (node ID,
	# orientation) pairs of the nodes you reach, and their orientations when you get
	# there, reading off of the node in the specified orientation.
	# We will call these pairs "traversals".
	traversals_after = defaultdict(set)


	with stream.open(str(vg_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.Graph()
			l.ParseFromString(data)
			for j in range(len(l.edge)):
				from_traversal = (getattr(l.edge[j], "from"), l.edge[j].from_start)
				to_traversal = (l.edge[j].to, l.edge[j].to_end)

				# Put the edge in the way it was read
				traversals_after[from_traversal].add(to_traversal)
				# Also store it in the other orientation, so you can follow it backward
				traversals_after[(to_traversal[0], not to_traversal[1])].add((from_traversal[0], not from_traversal[1]))

	for haptype in range(2):
		# for second haplotype
		prev_comp = -1
		hap1 =''
		hapseq1= defaultdict(list)
		hapseq2= defaultdict(list)
		haplotype_over_bubbles = defaultdict(list)
		start_node_to_bubble = defaultdict(list)

		for sample, superreads in sample_superreads.items():
			print(superreads)
			for v1, v2 in zip(*superreads):
				v = v1 if haptype == 0 else v2
				b = locus_branch_mapping[v.position][v.allele]
				# tmp stores the nodes over the haplotype path in a bubble
				tmp =list()
				tmp.append(b[0][0])
				for p,j in enumerate(b):
					tmp.append(j[-1])

				def dfs_path(start, goal, tmp):
					stack = [((start, True), [(start, True)]),((start, False), [(start, False)])]
					visited = set()
					visited.add(start)
					count = 0
					while stack:
						(traversal, path) = stack.pop()
						for next in traversals_after[traversal]:
							if count > 5000:
								break
							if next[0] in tmp and next not in visited:
								#if "{}_{}".format(vertex, next) in edge_connections_sign:
								if next[0] == goal:
									if len(path) == len(tmp) -1:
										return path + [next]
								else:
									count+=1
									visited.add(next)
									stack.append((next, path + [next]))
					return []

				path = dfs_path(tmp[0], tmp[-1], tmp)
				if len(path) != len(tmp):
					path = dfs_path(tmp[-1], tmp[0], tmp)

				# We need a function to flip a traversal
				def reverse_traversal(trav):
					return (trav[0], not trav[1])

				# We need a function to flip a path and all its traversals
				def reverse_path(to_reverse):
					return [reverse_traversal(t) for t in reversed(path)]

				# store the haplotype path with start or end as key
				if len(path) == len(tmp):
					haplotype_over_bubbles[path[0]] = path # from start
					haplotype_over_bubbles[reverse_traversal(path[-1])] = reverse_path(path)
					start_node_to_bubble[path[0]] = v.position
					start_node_to_bubble[reverse_traversal(path[-1])] = v.position

		# consider underlying graph as bidirected graph
		# start from canu contigs and make break them based on whatshap components
		# In bubbles, consider the haplotype path made up of nodes stored and whether to traverse the path in forward or backward, decide based on canu
		# at non-bubble region, consider path based on canu by considering the underlying graph.

		nodes_list = set()
		dummy_list = ['0']*1000
		orderalignment = defaultdict(list)
		orderalignment = defaultdict(lambda: [-1]*10000, orderalignment)
		with stream.open(str(canu_alignments), "rb") as istream:
			for data in istream:
				g = vg_pb2.Alignment()
				contig_nodes = []
				contig_nodes_blocks = []
				contig_nodes_seq = ''
				g.ParseFromString(data)
				save_nodes = []
				canu_nodes_toseq = defaultdict()
				for i in range(0,len(g.path.mapping)):
					index1 =  g.path.mapping[i].position.node_id
					orientation_canu = g.path.mapping[i].position.is_reverse
					save_nodes.append((index1, orientation_canu))
					canu_nodes_toseq[index1] = g.path.mapping[i].edit[0].sequence

				# What component was the last bubble in, if there was a last bubble
				prev_component = None

				it_val = 0
				already_done = set()
				for i in range(0,len(save_nodes)):
					if i >= it_val:
						index1 =  save_nodes[i][0]
						orientation_canu = save_nodes[i][1]

						# to take care of components, break when the bubbleid of previous and current is not equal
						if (index1, orientation_canu) in start_node_to_bubble:
							bubbleid = start_node_to_bubble[(index1, orientation_canu)]
							component = components[bubbleid]
							if prev_component is not None and component != prev_component:
								# We have moved to a new component of bubbles
								contig_nodes.append(contig_nodes_blocks)
								contig_nodes_blocks = []
								prev_component = component
							elif prev_component is None:
								# Remember the first component
								prev_component = component

						if (index1, orientation_canu) not in haplotype_over_bubbles:
							if orientation_canu == False:
								already_done.add(index1)
								contig_nodes_blocks.append(str(index1)+"_"+str(0))
							else:
								already_done.add(index1)
								contig_nodes_blocks.append(str(index1)+"_"+str(1))

						if (index1, orientation_canu) in haplotype_over_bubbles:
							if  haplotype_over_bubbles[(index1, orientation_canu)][-1] in save_nodes: # taking ordering from graph:
								for traversal in haplotype_over_bubbles[(index1, orientation_canu)][:-1]:
									if traversal[0] not in already_done:
										# Put each traversal that appears in the bubble in the contig node blocks
										# Except for the last one, which will be in the next bubble or in Canu again
										contig_nodes_blocks.append(str(traversal[0])+"_"+ ("1" if traversal[1] else "0"))
										already_done.add(traversal[0])


						if (index1, orientation_canu) in haplotype_over_bubbles and haplotype_over_bubbles[(index1, orientation_canu)][-1] in save_nodes:
							if save_nodes.index(haplotype_over_bubbles[(index1, orientation_canu)][-1]) > save_nodes.index(haplotype_over_bubbles[(index1, orientation_canu)][0]):
								# Skip to the last traversal in the bubble
								# It will also be shared by Canu
								it_val = save_nodes.index(haplotype_over_bubbles[(index1, orientation_canu)][-1]) # end node is not repeated

					else:
						# Don't do this Canu visit, it's part of a bubble we already did.
						continue


				contig_nodes.append(contig_nodes_blocks) # for the last one.
				print(contig_nodes)
				# build the contig sequence taking care of reverse complements for every canu contigs
				for j, contig_blocks in enumerate(contig_nodes):
					contig_nodes_seq = ''
					for i in contig_blocks:
						node = int(i.split("_")[0])
						if i.split("_")[1] == '1':
							contig_nodes_seq = contig_nodes_seq + reverse_complement(str(node_seq_list[node]))
						else:
							contig_nodes_seq = contig_nodes_seq + str(node_seq_list[node])
					pred_haplotigs_file.write(">seq" + str(j) +"_" + str(locus_file) + "_"+ str(haptype+1) + "\n")
					pred_haplotigs_file.write(contig_nodes_seq + '\n')

    

def run_phaseg(locus_file, gam_file, vg_file, canu_alignments, pred_haplotigs):
	recombrate=1.26
	max_coverage = 15
	all_heterozygous = False
	distrust_genotypes = True
	with ExitStack() as stack:
		node_seq_list, edge_connections = vg_graph_reader(vg_file)
		all_reads, alleles_per_pos, locus_branch_mapping = vg_reader(locus_file, gam_file)
		all_positions = sorted(all_reads.get_positions())
		all_components = find_components(all_positions, all_reads)
		blocks = defaultdict(list)
		print("all_components")
		for position, block_id in all_components.items():
			blocks[block_id].append(locus_branch_mapping[position][0][0][0])
		for k,v in blocks.items():
			print(k,v)
		print("all_components")
		

		selected_indices = readselection(all_reads, max_coverage)
		selected_reads = all_reads.subset(selected_indices)

		print("reads after read-selection")
		print(len(selected_reads))
		print("positions covered by atleast one read after read selection")
		print(len(selected_reads.get_positions()))

		accessible_positions = sorted(selected_reads.get_positions())
		
		pedigree = Pedigree(NumericSampleIds())
		# compute the number of alleles at each position.
		alleles_per_accessible_pos =[]
		genotype_likelihoods = []
		for pos in accessible_positions:
			if pos in alleles_per_pos:
				n_alleles = alleles_per_pos[pos]  
				possible_genotypes = n_alleles + ncr(n_alleles, 2)
				genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
		# random input of genotypes, since distrust_genotypes is always ON.
		pedigree.add_individual('individual0', [0]* len(accessible_positions), genotype_likelihoods)
		recombination_costs = uniform_recombination_map(recombrate, accessible_positions)
		# Finally, run phasing algorithm

		if len(selected_reads) == 0:
			print('I cannot phase this canu contig')
			sys.exit()

		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()

		cost = dp_table.get_optimal_cost()

		read_partitions = dp_table.get_optimal_partitioning()

		
		overall_components = find_components(accessible_positions, selected_reads)		
		n_phased_blocks = len(set(overall_components.values()))
		all_phased_blocks = len(set(all_components.values()))
		print('No. of phased blocks: %d', n_phased_blocks)
		largest_component = find_largest_component(overall_components)
		print('No. of blocks from all the reads: %d', all_phased_blocks)
		largest_component_all_reads = find_largest_component(all_components)
		if len(largest_component) > 0:
			print('Largest component contains %d variants',len(largest_component))
		if len(largest_component_all_reads) > 0:
			print('Largest component contains %d variants',len(largest_component_all_reads))

		sample = 0
		superreads, components = dict(), dict()
		superreads[sample] = superreads_list[0]
		components[sample] = overall_components
		generate_haplotigs(superreads, overall_components, node_seq_list, locus_branch_mapping, canu_alignments, vg_file, pred_haplotigs, locus_file)



def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT', help='read alignments in GAM file ')
	arg('vg_file', metavar='GRAPH', help='sequence graph')
	arg('canu_alignments', metavar='CANU_ALNS', help='contigs from canu.')
	arg('pred_haplotigs', metavar='PREED_HAPLOTIGS', help='write predicted haplotigs for every block')

def main(args):
	run_phaseg(**vars(args))
