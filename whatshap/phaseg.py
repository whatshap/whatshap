"""
create association between reads and bubbles.
"""
import pyfaidx
from xopen import xopen
import stream
import logging
from . import vg_pb2
from collections import Counter
import networkx as nx
from collections import defaultdict
from .core import ReadSet, Read
from functools import reduce
import operator as op
import networkx as nx
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

def slice_reads(reads, max_coverage):
	"""
	Iterate over all read in random order and greedily retain those reads whose
	addition does not lead to a local physical coverage exceeding the given threshold.
	Return a ReadSet containing the retained reads.

	max_coverage -- Slicing ensures that the (physical) coverage does not exceed max_coverage anywhere along the chromosome.
	reads -- a ReadSet
	"""
	
	SEED = 448
	random.seed(SEED)
	shuffled_indices = list(range(len(reads)))
	random.shuffle(shuffled_indices)

	position_list = reads.get_positions()
	logger.info('Found %d SNP positions', len(position_list))

	# dictionary to map SNP position to its index
	position_to_index = { position: index for index, position in enumerate(position_list) }

	# List of slices, start with one empty slice ...
	slices = [IndexSet()]
	# ... and the corresponding coverages along each slice
	slice_coverages = [CoverageMonitor(len(position_list))]
	skipped_reads = 0
	accessible_positions = set()
	for index in shuffled_indices:
		read = reads[index]
		# Skip reads that cover only one SNP
		if len(read) < 2:
			skipped_reads += 1
			continue
		for position, base, allele in read:
			accessible_positions.add(position)
		first_position, first_base, first_allele = read[0]
		last_position, last_base, last_allele = read[len(read)-1]
		begin = position_to_index[first_position]
		end = position_to_index[last_position] + 1
		slice_id = 0
		while True:
			# Does current read fit into this slice?
			if slice_coverages[slice_id].max_coverage_in_range(begin, end) < max_coverage:
				slice_coverages[slice_id].add_read(begin, end)
				slices[slice_id].add(index)
				break
			else:
				slice_id += 1
				# do we have to create a new slice?
				if slice_id == len(slices):
					slices.append(IndexSet())
					slice_coverages.append(CoverageMonitor(len(position_list)))
	logger.info('Skipped %d reads that only cover one SNP', skipped_reads)

	unphasable_snps = len(position_list) - len(accessible_positions)
	if position_list:
		logger.info('%d out of %d variant positions (%.1d%%) do not have a read '
			'connecting them to another variant and are thus unphasable',
			unphasable_snps, len(position_list),
			100. * unphasable_snps / len(position_list))

	# Print stats
	for slice_id, index_set in enumerate(slices):
		logger.info('Slice %d contains %d reads', slice_id, len(index_set))

	return reads.subset(slices[0])

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
	# create a dictionary of branches for each locus based on locus file.
	locus_branch_mapping=OrderedDict()
	locus_count=0
	prev_startsnarl = 0
	prev_endsnarl = 0
	locus_branch_mapping=defaultdict()
	locus_count=0
	prev_startsnarl = 0
	prev_startsnarl_orientation = -1
	prev_endsnarl = 0
	prev_endsnarl_orientation = -1
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
			cyclic_bubbles = [102838,102840,102846,102850,52424,52430,52708,52711,54914,54917,60635,60638,60965,60968,61857,61861,61906,61909,65760,65762,67841,67844,67858,67862,70509,70513,73378,73380,83218,83220,83224,83231,83676,83678,86581,86586,92007,92012,92467,92474,97403,97405,99187,99190]
			if l.snarl.end.node_id in cyclic_bubbles or l.snarl.start.node_id in cyclic_bubbles:
				continue
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
	#for i in [1, 2, 131, 132, 509, 6, 3, 646, 10, 12, 13, 269, 143, 16, 17, 657, 659, 407, 280, 667, 31, 672, 169, 301, 687, 560, 48, 691, 563, 693, 694, 569, 572, 317, 573, 574, 319, 577, 701, 579, 580, 325, 582, 583, 584, 201, 330, 586, 588, 589, 585, 590, 592, 593, 594, 337, 339, 597, 85, 599, 87, 345, 601, 67, 607, 608, 609, 482, 612, 614, 360, 65, 632, 581, 494, 371, 500, 501, 629, 120, 506, 380, 381]:
		#del locus_branch_mapping[i]
	#for i in [7, 8, 12, 15, 19, 20, 22, 23, 27, 29, 31, 32, 35, 38, 40, 42, 43, 45, 46, 52, 53, 59, 60, 61, 62, 63, 65, 69, 70, 71, 72, 78, 81, 82, 87, 91, 92, 94, 98, 100, 102, 104, 108, 114, 115, 118, 127, 128, 129, 142, 149, 156, 162, 163, 164, 165, 167, 170, 171, 172, 177, 182, 185, 186, 195, 198, 203, 211, 212, 213, 216, 223, 226, 227, 229, 231, 233, 235, 237, 238, 242, 243, 248, 249, 258, 259, 260, 266, 270, 271, 272, 273, 277, 278, 280, 287, 288, 290, 295, 298, 299, 301, 304, 305, 307, 308, 310, 311, 312, 315, 316, 317, 319, 322, 323, 328, 330, 335, 338, 340, 343, 346, 347, 348, 351, 354, 362, 364, 366, 367, 368, 369, 372, 373, 379, 383, 384, 385, 387, 391, 392, 394, 395, 396, 397, 399, 403, 404, 405, 406, 407, 409, 411, 413, 415, 417, 419, 422, 425, 428, 429, 432, 433, 436, 437, 438, 440, 442, 443, 444, 446, 452, 453, 458, 459, 462, 464, 465, 467]:
		#del locus_branch_mapping[i]
	#for i in [5, 19, 20, 22, 23, 27, 29, 31, 32, 35, 38, 40, 42, 43, 45, 46, 52, 53, 59, 60, 61, 62, 63, 65, 69, 70, 71, 72, 78, 81, 82, 87, 91, 92, 94, 98, 100, 102, 104, 108, 114, 115, 118, 127, 128, 129, 142, 149, 156, 162, 163, 164, 165, 167, 170, 171, 172, 177, 182, 185, 186, 195, 198, 203, 211, 212, 213, 216, 223, 226, 227, 229, 231, 233, 235, 237, 238, 242, 243, 248, 249, 258, 259, 260, 266, 270, 271, 272, 273, 277, 278, 280, 287, 288, 290, 295, 298, 299, 301, 304, 305, 307, 308, 310, 311, 312, 315, 316, 317, 319, 322, 323, 328, 330, 335, 338, 340, 343, 346, 347, 348, 351, 354, 362, 364, 366, 367, 368, 369, 372, 373, 379, 383, 384, 385, 387, 391, 392, 394, 395, 396, 397, 399, 403, 404, 405, 406, 407, 409, 411, 413, 415, 417, 419, 422, 425, 428, 429, 432, 433, 436, 437, 438, 440, 442, 443, 444, 446, 452, 453, 458, 459, 462, 464, 465, 467]:
		#del locus_branch_mapping[i]
	#for i in [7, 8, 12, 15, 19, 20, 22, 23, 27, 29, 31, 32, 35, 38, 40, 42, 43, 45, 46, 52, 53, 59, 60, 61, 62, 63, 65, 69, 70, 71, 72, 78, 81, 82, 87, 91, 92, 94, 98, 100, 102, 104, 108, 114, 115, 118, 127, 128, 129, 142, 149, 156, 162, 163, 164, 165, 167, 170, 171, 172, 177, 182, 185, 186, 195, 198, 203, 211, 212, 213, 216, 223, 226, 227, 229, 231, 233, 235, 237, 238, 242, 243, 248, 249, 258, 259, 260, 266, 270, 271, 272, 273, 277, 278, 280, 287, 288, 290, 295, 298, 299, 301, 304, 305, 307, 308, 310, 311, 312, 315, 316, 317, 319, 322, 323, 328, 330, 335, 338, 340, 343, 346, 347, 348, 351, 354, 362, 364, 366, 367, 368, 369, 372, 373, 379, 383, 384, 385, 387, 391, 392, 394, 395, 396, 397, 398, 399, 400, 401, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 422, 425, 428, 429, 432, 433, 436, 437, 438, 440, 442, 443, 444, 446, 452, 453, 458, 459, 462, 464, 465, 467]:
		#del locus_branch_mapping[i]
	#for i in [7, 8, 12, 15, 19, 20, 22, 23, 27, 29, 31, 32, 35, 38, 40, 42, 43, 45, 46, 52, 53, 59, 60, 61, 62, 63, 65, 69, 70, 71, 72, 78, 81, 82, 87, 91, 92, 94, 98, 100, 102, 104, 108, 114, 115, 118, 127, 128, 129, 142, 149, 156, 162, 163, 164, 165, 167, 170, 171, 172, 177, 182, 185, 186, 195, 198, 203, 211, 212, 213, 216, 223, 226, 227, 229, 231, 233, 235, 237, 238, 242, 243, 248, 249, 258, 259, 260, 266, 270, 271, 272, 273, 277, 278, 280, 287, 288, 290, 295, 298, 299, 301, 304, 305, 307, 308, 310, 311, 312, 315, 316, 317, 319, 322, 323, 328, 330, 335, 338, 340, 343, 346, 347, 348, 351, 354, 362, 364, 366, 367, 368, 369, 372, 373, 379, 383, 384, 385, 387, 391, 392, 394, 395, 396, 397, 399, 403, 404, 405, 406, 407, 409, 411, 413, 415, 417, 419, 422, 425, 428, 429, 432, 433, 436, 437, 438, 440, 442, 443, 444, 446, 452, 453, 458, 459, 462, 464, 465, 467]:
		#del locus_branch_mapping[i]
	#for i in [395, 396]:
		#del locus_branch_mapping[i]
	

	#print(locus_branch_mapping)
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
	#print(reverse_mapping)
	print(locus_branch_mapping)

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
			if score < 0.75:
				continue
			read=Read(g.name, 0, 0, 0) # create read for each read alignment
			#readnames= ["S1_Y12_290","S1_SK1_290","S1_Y12_430","S1_SK1_657","S1_Y12_139","S1_Y12_427","S1_SK1_427","S1_Y12_657","S1_SK1_588","S1_Y12_588","S1_SK1_139","S1_SK1_430","S1_Y12_76","S1_Y12_463","S1_SK1_463","S1_SK1_76"]
			#readnames = ["S1_Y12_259"]
			#if g.name not in readnames:
				#continue
			print(g.name)
			prev_tmp=[]
			prev_locus= -1
			locus = -1
			#for i in range(0,len(g.path.mapping)):
				#if g.path.mapping[i].position.is_reverse != val1:
					#val1 = False
					#break
				#else:
					#count1 = count1 +1
					
			#if count1 == len(g.path.mapping):
				#count = count+1
				##print(g.name)
				
			#for i in range(0,len(g.path.mapping)):
				#if g.path.mapping[i].position.is_reverse != val2:
					#val2 = True
					#break
				#else:
					#count2 = count2 +1
					
			#if count2 == len(g.path.mapping):
				#count = count+1
				##print(g.name)
			#print(val1)
			#print(val2)
			#if val1 ==val2:
			for i in range(0,len(g.path.mapping)-1):
			#for i in g.path.mapping: # go over the mapping in a read
				# TODO: check for forward or reverse strand, we may not need it for DAG.
				edge1 = tuple((g.path.mapping[i].position.node_id, g.path.mapping[i+1].position.node_id)) # go over nodes in a mapping
				edge2 = tuple((g.path.mapping[i+1].position.node_id, g.path.mapping[i].position.node_id)) # go over nodes in a mapping
				if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
					if edge1 in reverse_mapping:
						qualities = [10]* reverse_mapping[edge1][0][2]
						node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge1]] # consider (locus, branch)
					else:
						qualities = [10]* reverse_mapping[edge2][0][2]
						node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge2]]
					tmp = [x for x in node_inf]
					if prev_locus != tmp[0][0]:
						prev_tmp = tmp
						prev_locus = tmp[0][0]
						
					interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
					if len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
						qualities[interset_tmp[0][1]] = 0
						if i== len(g.path.mapping)-2:
							read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
						else:
							next_edge1 = tuple((g.path.mapping[i+1].position.node_id, g.path.mapping[i+2].position.node_id))
							next_edge2 = tuple((g.path.mapping[i+2].position.node_id, g.path.mapping[i+1].position.node_id))

							if next_edge1 not in reverse_mapping and next_edge2 not in reverse_mapping:
								read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)    

						locus= interset_tmp[0][0]
						
					#if prev_locus!=locus:
						#prev_tmp = []
					#else:
						#for i in tmp:
							#prev_tmp.append(i)
					#prev_locus = locus
			print(read)

			if len(read) >= 2:
				readset.add(read)
	print("non-shattered")
	print(count)
	#print(readset)
	readset1=ReadSet()
	tmp_duplicated=set()
	for read in readset:
		if read.sort() ==1:
			duplicated = duplicated +1
			tmp=[]
			for variant in read:
				tmp.append(variant.position)
			print("duplicated variant")
			x = [item for item, count in collections.Counter(tmp).items() if count > 1]
			for a in x:
				tmp_duplicated.add(a)
			continue
		else:
			readset1.add(read)
	print("length of duplicated bubbles")
	print(tmp_duplicated)
	print(len(list(tmp_duplicated)))


	readset1.sort()
	#print("******")
	#for i,read in enumerate(readset1):
		#for j,variant in enumerate(read):
			#print(str(i)+" "+str(variant.position)+" "+str(variant.allele)+ " "+"10")
	#print("******")
	print("duplicated")
	print(duplicated)
	print("reads considered before read-selection")
	print(len(readset1))
	return readset1, alleles_per_pos, locus_branch_mapping


"""
consider only top-k paths from complex bubbles using Yen's algorithm for later.
"""

def k_shortest_paths(G, source, target, k=1, weight='weight'):
	"""Returns the k-shortest paths from source to target in a weighted graph G.
	"""
	if source == target:
		return ([0], [[source]]) 
	   
	length, path = nx.single_source_dijkstra(G, source, target, weight=weight)
	if target not in length:
		raise nx.NetworkXNoPath("node %s not reachable from %s" % (source, target))
		
	lengths = [length[target]]
	paths = [path[target]]
	c = count()		
	B = []						
	G_original = G.copy()	
	
	for i in range(1, k):
		for j in range(len(paths[-1]) - 1):			
			spur_node = paths[-1][j]
			root_path = paths[-1][:j + 1]
			
			edges_removed = []
			for c_path in paths:
				if len(c_path) > j and root_path == c_path[:j + 1]:
					u = c_path[j]
					v = c_path[j + 1]
					if G.has_edge(u, v):
						edge_attr = G.edge[u][v]
						G.remove_edge(u, v)
						edges_removed.append((u, v, edge_attr))
			
			for n in range(len(root_path) - 1):
				node = root_path[n]
				# out-edges
				for u, v, edge_attr in G.copy().edges_iter(node, data=True):
					G.remove_edge(u, v)
					edges_removed.append((u, v, edge_attr))
				
				if G.is_directed():
					# in-edges
					for u, v, edge_attr in G.in_edges_iter(node, data=True):
						G.remove_edge(u, v)
						edges_removed.append((u, v, edge_attr))
			
			spur_path_length, spur_path = nx.single_source_dijkstra(G, spur_node, target, weight=weight)			
			if target in spur_path and spur_path[target]:
				total_path = root_path[:-1] + spur_path[target]
				total_path_length = get_path_length(G_original, root_path, weight) + spur_path_length[target]				
				heappush(B, (total_path_length, next(c), total_path))
				
			for e in edges_removed:
				u, v, edge_attr = e
				G.add_edge(u, v, edge_attr)
					   
		if B:
			(l, _, p) = heappop(B)
			lengths.append(l)
			paths.append(p)
		else:
			break
	
	return (lengths, paths)

def get_path_length(G, path, weight='weight'):
	length = 0
	if len(path) > 1:
		for i in range(len(path) - 1):
			u = path[i]
			v = path[i + 1]
			
			length += G.edge[u][v].get(weight, 1)
	
	return length 

"""
To generate two haplotype sequences.
Assumption: positions in one component occur consecutive
"""
def generate_hap_contigs(sample_superreads, sample_components, node_seq_list, locus_branch_mapping, edge_connections):
	sample = 0
	components = sample_components[sample]

	# TODO: sort components dictionary by value.
	prev_comp = -1
	hap1 =''
	hapseq1= defaultdict(list)
	hapseq2= defaultdict(list)

	for sample, superreads in sample_superreads.items():
		for v1, v2 in zip(*superreads):	
			#TODO: handle ambiguos cases
			b = locus_branch_mapping[v1.position][v1.allele]
			if v1.allele == -2:
				b = locus_branch_mapping[v1.position][1]
			tmp =list()
			if v1.position == 1:
				tmp.append(b[0][0])
			for p,j in enumerate(b):
				tmp.append(j[-1])


			for i in tmp:
				comp = components[v1.position]
				hapseq1[comp].append(node_seq_list[i])

			current_node = tmp[-1]
			while current_node in edge_connections and len(edge_connections[current_node]) ==1:
				hapseq1[comp].append(node_seq_list[edge_connections[current_node][0]])
				current_node = edge_connections[current_node][0]

			#TODO: handle ambiguos cases
			b = locus_branch_mapping[v2.position][v2.allele]
			if v2.allele == -2:
				b = locus_branch_mapping[v2.position][1]
			tmp =list()
			if v2.position == 1:
				tmp.append(b[0][0])
			for p,j in enumerate(b):
				tmp.append(j[-1])

			for i in tmp:
				comp = components[v2.position]
				hapseq2[comp].append(node_seq_list[i])
					
			current_node = tmp[-1]
			while current_node in edge_connections and len(edge_connections[current_node]) ==1:
				hapseq2[comp].append(node_seq_list[edge_connections[current_node][0]])
				current_node = edge_connections[current_node][0]


	for k,v in hapseq1.items():
		hap1=''
		hap2=''
		for v in hapseq1[k]:
			hap1=hap1+v
		for v in hapseq2[k]:
			hap2=hap2+v
		print("I am in component" + str(k))
		print(hap1)
		print(hap2)
		
# partition all set of reads by considering the haplotypes from most significant reads. 
def haplotag(pred_superreads, read_set, components, iteration):
	phases = []
	
	for s1,s2 in zip(*pred_superreads):
		VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'position', 'phase1', 'phase2'])
		extract_phase = VariantCallPhase(components[s1.position], s1.position, s1.allele, s2.allele)
		phases.append(extract_phase) #TODO: check it
		
	variantpos_to_allele1 = {
		v.position:int(v.allele) for v, v2 in zip(*pred_superreads) if v.allele!=-2
	}
	variantpos_to_allele2 = {
		v2.position:int(v2.allele) for v, v2 in zip(*pred_superreads) if v2.allele!=-2
	}
	
	variantpos_to_phaseset1 = {
		v.position:components[v.position] for v, v2 in zip(*pred_superreads) if v.allele!=-2
	}
	variantpos_to_phaseset2 = {
		v2.position:components[v2.position] for v, v2 in zip(*pred_superreads) if v2.allele!=-2
	}
	read_to_haplotype = {}
	#read_set = read_reads(readset_reader, chromosome_name, variants, sample, fasta)
	for read in read_set:
		# mapping: phaseset --> phred scaled difference between costs of assigning reads to haplotype 0 or 1
		haplotype_costs = defaultdict(int)
		haplotype_costs[0] = 0
		haplotype_costs[1] = 0
		for v in read:
			if v.position not in variantpos_to_allele1 or v.position not in variantpos_to_allele2:
				continue
			phaseset1 = variantpos_to_allele1[v.position]
			phaseset2 = variantpos_to_allele2[v.position]
	
			if v.allele != phaseset1:
				haplotype_costs[0] += 10
			#else:
				#haplotype_costs[0] -= 10
			if v.allele != phaseset2:
				haplotype_costs[1] += 10
			#else:
				#haplotype_costs[1] -= 10
			
		l =[]
		for k,v in haplotype_costs.items():
			l.append(k)
			l.append(v)
		#l = list(haplotype_costs.items())
		if l[1] < l[3]:
			read_to_haplotype[read.name] = (0, l[1], l[0])
		else:
			read_to_haplotype[read.name] = (1, l[3], l[2])
		#l.sort(key=lambda t:-abs(t[1]))
		#phaseset, quality = l[0]
		#if quality != 0:
			#haplotype = 0 if quality > 0 else 1
			#read_to_haplotype[read.name] = (haplotype, abs(quality), phaseset)
			
	f = open('predicted_all_read_partionting' + str(iteration), 'w')
	accessible_positions = sorted(read_set.get_positions())
	overall_components = find_components(accessible_positions, read_set)
	for read in read_set:
		phaseset = components[read[0].position] + 1
		print(read.name, read_to_haplotype[read.name][2], read_to_haplotype[read.name][0], phaseset, file =f )
	#for k,v in read_to_haplotype.items():
		#print(k, v[2], v[0], file=f)


def compute_read_partitioning_accuracy(predicted_file, true_file):
	true_hap1 =[]
	true_hap2= []
	f = open(true_file, 'r')
	
	for line in open(true_file):
		if line.split("\t")[1] ==1:
			true_hap1.append(ine.split("\t")[0])
		else:
			true_hap2.append(ine.split("\t")[0])
		
	pred_hap1 = defaultdict(list)
	pred_hap2 = defaultdict(list)
	total=0
	blocks =set()
	for line in open(predicted_file):
		tokens= line.split(" ")
		blocks.add(token[1])
		if tokens[2] ==1:
			pred_hap1[tokens[1]].append(tokens[0])
		else:
			pred_hap2[tokens[1]].append(tokens[0])
		total+=1
	
	count = 0
	for k in blocks:
		count += min(set(pred_hap1[k]).intersection(set(true_hap1)), set(pred_hap1[k]).intersection(set(true_hap2)))
		count += min(set(pred_hap2[k]).intersection(set(true_hap1)), set(pred_hap2[k]).intersection(set(true_hap2)))
	percent_partitionining_accuracy =  count/total
	
count_width = 9
class SwitchFlips:
	def __init__(self, switches=0, flips=0):
		self.switches = switches
		self.flips = flips
	def __iadd__(self, other):
		self.switches += other.switches
		self.flips += other.flips
		return self
	def __repr__(self):
		return 'SwitchFlips(switches={}, flips={})'.format(self.switches, self.flips)
	def __str__(self):
		return '{}/{}'.format(self.switches, self.flips)


class PhasingErrors:
	def __init__(self, switches=0, hamming=0, switch_flips=None):
		self.switches = switches
		self.hamming = hamming
		self.switch_flips = SwitchFlips() if switch_flips is None else switch_flips
	def __iadd__(self, other):
		self.switches += other.switches
		self.hamming += other.hamming
		self.switch_flips += other.switch_flips
		return self
	def __repr__(self):
		return 'PhasingErrors(switches={}, hamming={}, switch_flips={})'.format(self.switches, self.hamming, self.switch_flips)


def complement(s):
	t = { '0': '1', '1':'0' }
	return ''.join(t[c] for c in s)


def hamming(s0, s1):
	assert len(s0) == len(s1)
	return sum( c0!=c1 for c0, c1 in zip(s0, s1) )


def switch_encoding(phasing):
	return ''.join( ('0' if phasing[i-1]==phasing[i] else '1') for i in range(1,len(phasing)) )


def compute_switch_flips(phasing0, phasing1):
	assert len(phasing0) == len(phasing1)
	s0 = switch_encoding(phasing0)
	s1 = switch_encoding(phasing1)
	result = SwitchFlips()
	switches_in_a_row = 0
	for i, (p0, p1) in enumerate(zip(s0, s1)):
		if p0 != p1:
			switches_in_a_row += 1
		if (i + 1 == len(s0)) or (p0 == p1):
			result.flips += switches_in_a_row // 2
			result.switches += switches_in_a_row % 2
			switches_in_a_row = 0
	if False:
		print('switch_flips():')
		print('   phasing0={}'.format(phasing0))
		print('   phasing1={}'.format(phasing1))
		print('         s0={}'.format(s0))
		print('         s1={}'.format(s1))
		print('   switches={}, flips={}'.format(result.switches, result.flips))
	return result


def compare_block(phasing_pred1, phasing_pred2, phasing_true1, phasing_true2):
	"""Input are two strings over {0,1}. Output is a PhasingErrors object."""
	return PhasingErrors(
		switches = min(min(hamming(switch_encoding(phasing_pred1), switch_encoding(phasing_true1)), hamming(switch_encoding(phasing_pred1), switch_encoding(phasing_true2))), min(hamming(switch_encoding(phasing_pred2), switch_encoding(phasing_true1)), hamming(switch_encoding(phasing_pred2), switch_encoding(phasing_true2)))),
		hamming = min(min(hamming(phasing_pred1, phasing_true1), hamming(phasing_pred1, phasing_true2)), min(hamming(phasing_pred2, phasing_true1), hamming(phasing_pred2, phasing_true2))), switch_flips=SwitchFlips()
		)
      
def fraction2percentstr(nominator, denominator):
	if denominator == 0:
		return '--'
	else:
		return '{:.2f}%'.format(nominator*100.0/denominator)


def safefraction(nominator, denominator):
	if denominator == 0:
		return float('nan')
	else:
		return nominator/denominator


def create_bed_records(chromosome, phasing0, phasing1, positions, annotation_string):
	"""Determines positions of switch errors between two phasings
	and yields one BED record per switch error (encoded as a tuple).
	The annotation_string is added to each record."""
	assert len(phasing0) == len(phasing1) == len(positions)
	switch_encoding0 = switch_encoding(phasing0)
	switch_encoding1 = switch_encoding(phasing1)
	for i, (sw0, sw1) in enumerate(zip(switch_encoding0, switch_encoding1)):
		if sw0 != sw1:
			yield (chromosome, positions[i]+1, positions[i+1]+1, annotation_string)


def print_errors(errors, phased_pairs, print_hamming=False):
	print('    phased pairs of variants assessed:', str(phased_pairs).rjust(count_width))
	print('                        switch errors:', str(errors.switches).rjust(count_width))
	print('                    switch error rate:', fraction2percentstr(errors.switches, phased_pairs).rjust(count_width))
	print('            switch/flip decomposition:', str(errors.switch_flips).rjust(count_width) )
	print('                     switch/flip rate:', fraction2percentstr(errors.switch_flips.switches+errors.switch_flips.flips, phased_pairs).rjust(count_width))
      
pairwise_comparison_results_fields = [
	'intersection_blocks',
	'covered_variants',
	'all_assessed_pairs',
	'all_switches',
	'all_switch_rate',
	'all_switchflips',
	'all_switchflip_rate',
	'largestblock_assessed_pairs',
	'largestblock_switches',
	'largestblock_switch_rate',
	'largestblock_switchflips',
	'largestblock_switchflip_rate',
	'largestblock_hamming',
	'largestblock_hamming_rate'
]
PairwiseComparisonResults = namedtuple('PairwiseComparisonResults', pairwise_comparison_results_fields)

# TO compare the true and predicted superreads.
def compare(predicted_superreads, true_superreads, components):
	phases_pred1 = defaultdict()
	phases_pred2 = defaultdict()
	phases_true1 = defaultdict()
	phases_true2 = defaultdict()
	
	for s1,s2 in zip(*predicted_superreads):
		#VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'position','phase'])
		#extract_phase = VariantCallPhase(components[s1.position], s1.position, s1.allele)
		if s1.allele !=-2 and s2.allele !=-2:
			phases_pred1[s1.position] = s1.allele
			phases_pred2[s1.position] = s2.allele 

	for s1,s2 in zip(*true_superreads):
		#VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'position','phase'])
		#extract_phase = VariantCallPhase(components[s1.position], s1.position, s1.allele)
		if s1.allele !=-2 and s2.allele !=-2:
			phases_true1[s1.position] = s1.allele
			phases_true2[s1.position] = s2.allele 

	
	blocks = defaultdict(list)
	for k,v in phases_pred1.items():
		if phases_pred1[k]!=-2 or phases_pred2[k]!=-2 or phases_true1[k]!=-2 or phases_true2[k]!=-2:
			blocks[components[k]].append(k)
	intersection_block_variants = sum(len(b) for b in blocks.values() if len(b) > 1)
	
	longest_block = None
	longest_block_errors = None
	phased_pairs = 0
	bed_records = []

	total_errors = PhasingErrors()
	for block in blocks.values():
		if len(block) < 2:
			continue
		phasing_pred1 = ''.join( str(phases_pred1[i]) for i in block )
		phasing_pred2 = ''.join( str(phases_pred2[i]) for i in block )
		phasing_true1 = ''.join( str(phases_true1[i]) for i in block )
		phasing_true2 = ''.join( str(phases_true2[i]) for i in block )
		block_positions = [ i for i in block ]
		errors = compare_block(phasing_pred1, phasing_pred2, phasing_true1, phasing_true2)
		#bed_records.extend(create_bed_records(chromosome, phasing0, phasing1, block_positions, '{}<-->{}'.format(*dataset_names)))
		total_errors += errors
		phased_pairs += len(block) - 1
		if (longest_block is None) or (len(block) > longest_block):
			longest_block = len(block)
			longest_block_errors = errors
	print('              ALL INTERSECTION BLOCKS:', '-'*count_width)
	print_errors(total_errors, phased_pairs)
	print('           LARGEST INTERSECTION BLOCK:', '-'*count_width)
	print_errors(longest_block_errors, longest_block-1)
	print('                     Hamming distance:', str(longest_block_errors.hamming).rjust(count_width))
	print('                 Hamming distance [%]:', fraction2percentstr(longest_block_errors.hamming, longest_block).rjust(count_width))
	print('%phased accuracy', safefraction(total_errors.switches, phased_pairs))
	print('number of blocks', len(blocks.values()))
	return PairwiseComparisonResults(
		intersection_blocks = len(blocks.values()),
		covered_variants = intersection_block_variants,
		all_assessed_pairs = phased_pairs,
		all_switches = total_errors.switches,
		all_switch_rate = safefraction(total_errors.switches, phased_pairs),
		all_switchflips = -1,
		all_switchflip_rate = -1,
		largestblock_assessed_pairs = longest_block-1,
		largestblock_switches = longest_block_errors.switches,
		largestblock_switch_rate = safefraction(longest_block_errors.switches, longest_block - 1),
		largestblock_switchflips = -1,
		largestblock_switchflip_rate = -1,
		largestblock_hamming = -1,
		largestblock_hamming_rate = -1
	), bed_records



def run_phaseg(locus_file, gam_file, vg_file, true_haps):
	"""
	Run WhatsHap.

	gam_file -- path to GAM file
	locus_file -- path to LOCUS file
	"""
	recombrate=1.26
	max_coverage = 20
	all_heterozygous = False
	distrust_genotypes = True
	with ExitStack() as stack:
		node_seq_list, edge_connections = vg_graph_reader(vg_file)
		all_reads, alleles_per_pos, locus_branch_mapping = vg_reader(locus_file, gam_file)

		#print(all_reads)
		selected_indices = readselection(all_reads, max_coverage)
		selected_reads = all_reads.subset(selected_indices)
		#selected_reads = slice_reads(all_reads, max_coverage)
		#print('positions from all reads')
		#print(len(all_reads.get_positions()))
		print("reads after read-selection")
		print(len(selected_reads))
		print("positions covered by atleast one read after read selection")
		print(len(selected_reads.get_positions()))

		accessible_positions = sorted(selected_reads.get_positions())
		print("readset after read_selection")
		#for read in selected_reads:
			#print(read.name)
		pedigree = Pedigree(NumericSampleIds())
		# compute the number of alleles at each position.
		alleles_per_accessible_pos =[]
		genotype_likelihoods = []
		for pos in accessible_positions:
			if pos in alleles_per_pos:
				n_alleles = alleles_per_pos[pos]  
				possible_genotypes = n_alleles +  ncr(n_alleles, 2)
				genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
		# random input of genotypes, since distrust_genotypes is always ON.
		pedigree.add_individual('individual0', [0]* len(accessible_positions), genotype_likelihoods)
		recombination_costs = uniform_recombination_map(recombrate, accessible_positions)
		# Finally, run phasing algorithm
		#print(selected_reads)
		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()

		cost = dp_table.get_optimal_cost()
		print(superreads_list[0])
		#print(cost)
		read_partitions = dp_table.get_optimal_partitioning()
		#print(read_partitions)
		
		## To generate the connected components and corresponding haplotypes.
		print("in components")
		f = open('predicted_read_partionting', 'w')
		overall_components = find_components(accessible_positions, selected_reads)
		for read, haplotype in zip(selected_reads, read_partitions):
			phaseset = overall_components[read[0].position] + 1
			print(read.name, phaseset, haplotype, file=f)
		#phaset is blockid

		n_phased_blocks = len(set(overall_components.values()))
		print('No. of phased blocks: %d', n_phased_blocks)
		largest_component = find_largest_component(overall_components)
		if len(largest_component) > 0:
			print('Largest component contains %d variants',len(largest_component))
		
		## To generate contig sequences
		sample = 0
		superreads, components = dict(), dict()
		superreads[sample] = superreads_list[0]
		components[sample] = overall_components
		
		# evaluation partition all the reads based on one iteration
		print('partition all the reads based on haplotypes from one iteration')
		haplotag(superreads_list[0], all_reads, overall_components, 1)


		#generate_hap_contigs(superreads, components, node_seq_list, locus_branch_mapping, edge_connections)
		
		#For phasing accuracy, read true haps and generate corresponding superreads
		all_reads_true, alleles_per_pos_true, locus_branch_mapping_true = vg_reader(locus_file, true_haps)
		# Finally, run phasing algorithm for true haplotypes
		dp_table_true = PedigreeDPTable(all_reads_true, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list_true, transmission_vector_true = dp_table_true.get_super_reads()
		# to compute the phasing accuracy
		compare(superreads_list[0], superreads_list_true[0], overall_components)
		# To perform iterative whatshap phasing
		remaining_reads =[]
		for read in all_reads:
			remaining_reads.append(read.name)
		prev_superreads = superreads_list[0]
		for read in selected_reads:
			remaining_reads.remove(read.name)
		while len(remaining_reads)>0:
			print('iteration')
			iterative_reaset =  ReadSet()
			for read in all_reads:
				if read.name in remaining_reads:
					iterative_reaset.add(read)

				
			selected_indices = readselection(iterative_reaset, max_coverage)
			selected_reads = iterative_reaset.subset(selected_indices)
			for read in prev_superreads:
				selected_reads.add(read)
				remaining_reads.append(read.name)
			accessible_positions = sorted(selected_reads.get_positions())
			selected_reads.sort()
			pedigree = Pedigree(NumericSampleIds())
			# compute the number of alleles at each position.
			alleles_per_accessible_pos =[]
			genotype_likelihoods = []
			for pos in accessible_positions:
				if pos in alleles_per_pos:
					n_alleles = alleles_per_pos[pos]  
					possible_genotypes = n_alleles +  ncr(n_alleles, 2)
					genotype_likelihoods.append(None if all_heterozygous else PhredGenotypeLikelihoods([0]* possible_genotypes))
			# random input of genotypes, since distrust_genotypes is always ON.
			pedigree.add_individual('individual0', [0]* len(accessible_positions), genotype_likelihoods)
			recombination_costs = uniform_recombination_map(recombrate, accessible_positions)
			# Finally, run phasing algorithm
			#print(selected_reads)
			dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
			superreads_list, transmission_vector = dp_table.get_super_reads()
			for read in selected_reads:
				remaining_reads.remove(read.name)
			prev_superreads = superreads_list[0]
			
		print('I am final')
		accessible_positions = sorted(all_reads.get_positions())
		overall_components = find_components(accessible_positions, all_reads)
		haplotag(superreads_list[0], all_reads, overall_components, all_iter)
		compare(superreads_list[0], superreads_list_true[0], overall_components)
		print(superreads_list[0])


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT', help='read alignments in GAM file ')
	arg('vg_file', metavar='GRAPH', help='node-sequence association')
	arg('true_haps', metavar='TRUE_HAPS', help='compare phasing with true haps in GAM format.')
	#TODO: add k parameter

def main(args):
	run_phaseg(**vars(args))

