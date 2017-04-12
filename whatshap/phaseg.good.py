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


from contextlib import ExitStack
from .vcf import VcfReader, PhasedVcfWriter
from . import __version__
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
	locus_branch_mapping=defaultdict()
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
			print(l.visits)
			# handle forward and backward case of nodes
			current_startsnarl = l.snarl.start.node_id
			current_startsnarl_orientation = l.snarl.start.backward
			current_endsnarl = l.snarl.end.node_id
			current_endsnarl_orientation = l.snarl.end.backward
			path_in_bubble =[]
			print(l.visits[0].snarl.start)
			path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
			path_in_bubble.append(tuple ((l.visits[0].node_id, l.snarl.end.node_id)))
			for i in range(0,len(l.visits)-1):
				path_in_bubble.append(tuple((l.visits[i]), l.visits[i+1]))
				#path_in_bubble.append(tuple ((x.node_id,x.snarl.end) for x in l.visits))
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
	print(locus_branch_mapping)


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
						print(j)
						reverse_mapping[j].append([k,i, len(v)]) # in complex bubbles, a node can map to multiple branches.

	print(reverse_mapping)


	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	readset=ReadSet()
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			# hard-coded source id, mapping quality and other values.
			read=Read(g.name, 0, 0, 0) # create read for each read alignment
			prev_tmp=[]
			prev_locus= -1
			for i in range(0,len(g.path.mapping)-1):
			#for i in g.path.mapping: # go over the mapping in a read
				# TODO: check for forward or reverse strand, we may not need it for DAG.
				edge = tuple((g.path.mapping[i].position.node_id, g.path.mapping[i+1].position.node_id)) # go over nodes in a mapping
				if edge in reverse_mapping: # handle start and sink node.
					qualities = [10]* reverse_mapping[edge][0][2]
					node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge]] # consider (locus, branch)
					tmp = [x for x in node_inf]
					interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
					if len(tmp)==1 and len(prev_tmp)==0: # simple bubble with multiple paths.
						qualities[tmp[0][1]] = 0
						read.add_variant(tmp[0][0], tmp[0][1], qualities) # if any new locus of branch is encountered, enter a variant in read.
						locus = tmp[0][0]
					elif len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
						qualities[interset_tmp[0][1]] = 0
						read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
						locus= interset_tmp[0][0]
					elif len(prev_tmp) > 0 and len(tmp) ==1: # complicated bubbles, but one node can uniquely identify the branch.
						qualities[tmp[0][1]] = 0
						read.add_variant(tmp[0][0], tmp[0][1], qualities)
						locus = tmp[0][0]

					if prev_locus!=locus:
						prev_tmp = []
					else:
						for i in tmp:
							prev_tmp.append(i)

			if len(read) >= 2:
				readset.add(read)
	for read in readset:
			read.sort()
	readset.sort()
	return readset, alleles_per_pos

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

def run_phaseg(locus_file, gam_file):
	"""
	Run WhatsHap.

	gam_file -- path to GAM file
	locus_file -- path to LOCUS file
	"""
	recombrate=1.26
	max_coverage = 15
	all_heterozygous = False
	distrust_genotypes = True
	with ExitStack() as stack:
		all_reads, alleles_per_pos = vg_reader(locus_file, gam_file)
		print(all_reads)
		selected_indices = readselection(all_reads, max_coverage)
		selected_reads = all_reads.subset(selected_indices)
		accessible_positions = sorted(selected_reads.get_positions())
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
		print(selected_reads)
		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()
		cost = dp_table.get_optimal_cost()
		print(superreads_list[0])
		print(cost)
		read_partitions = dp_table.get_optimal_partitioning()
		print(read_partitions)


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT', help='read alignments in GAM file ')
	#TODO: add k parameter

def main(args):
	run_phaseg(**vars(args))

