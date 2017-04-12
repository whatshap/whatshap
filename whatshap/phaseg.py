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
			# handle forward and backward case of nodes
			current_startsnarl = l.snarl.start.node_id
			current_startsnarl_orientation = l.snarl.start.backward
			current_endsnarl = l.snarl.end.node_id
			current_endsnarl_orientation = l.snarl.end.backward
			path_in_bubble =[]
			if len(l.visits) ==1 and l.visits[0].backward == False: # consider only hets 
				path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
				path_in_bubble.append(tuple ((l.visits[0].node_id, l.snarl.end.node_id)))
			#if len(l.visits) >1:
				#for i in range(0,len(l.visits)-1):
					#path_in_bubble.append(tuple((l.visits[i].node_id, l.visits[i+1].node_id)))
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
				#if locus_count > 5:
					#break
	#print(locus_branch_mapping)
	print('The number of hets:')
	print(len(locus_branch_mapping))
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

	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	readset=ReadSet()
	count =0
	duplicated = 0
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			# hard-coded source id, mapping quality and other values.
			val1 = True
			val2 = False

			count1 =0
			count2=0
			#if g.name == "m150910_184604_00127_c100822732550000001823176011031536_s1_p0/29973/0_4989" or g.name =="m150910_184604_00127_c100822732550000001823176011031536_s1_p0/83922/31432_33600" or g.name =="m150912_012316_00127_c100861772550000001823190702121672_s1_p0/131989/32651_34485" or g.name =="m150815_025203_00127_c100823152550000001823177111031543_s1_p0/68101/28755_35578" or g.name =="m150811_224417_00127_c100823122550000001823177111031570_s1_p0/18019/834_5339" or g.name =="m150811_224417_00127_c100823122550000001823177111031570_s1_p0/41465/0_11463" or g.name == "m150911_220012_00127_c100861772550000001823190702121671_s1_p0/139294/42242_42775" or g.name =="m150811_092723_00127_c100844062550000001823187612311514_s1_p0/90852/0_2265" or g.name =="m150814_233250_00127_c100823152550000001823177111031542_s1_p0/58220/0_3937" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/105221/21433_25660" or g.name == "m150912_012316_00127_c100861772550000001823190702121672_s1_p0/120236/0_7569" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/147019/0_8535" or g.name =="m150814_233250_00127_c100823152550000001823177111031542_s1_p0/69772/12374_15347" or g.name == "m150815_061119_00127_c100823152550000001823177111031544_s1_p0/141285/0_2434" or g.name == "m150911_220012_00127_c100861772550000001823190702121671_s1_p0/9290/11072_15672" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/30157/0_4006" or g.name == "m150910_184604_00127_c100822732550000001823176011031536_s1_p0/59834/3718_9567" or g.name =="m150814_233250_00127_c100823152550000001823177111031542_s1_p0/47790/0_17299" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/141068/3846_8373" or g.name == "m150814_201337_00127_c100823152550000001823177111031541_s1_p0/143293/0_12359" or g.name == "m150911_220012_00127_c100861772550000001823190702121671_s1_p0/39147/0_14935" or g.name == "m150910_220412_00127_c100822732550000001823176011031537_s1_p0/52956/0_5569" or g.name == "m150911_220012_00127_c100861772550000001823190702121671_s1_p0/110288/27962_29316" or g.name =="m150814_201337_00127_c100823152550000001823177111031541_s1_p0/96793/4271_8504" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/22285/12653_25111" or g.name == "m150811_092723_00127_c100844062550000001823187612311514_s1_p0/65923/0_10033" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/121604/13638_16688" or g.name == "m150811_092723_00127_c100844062550000001823187612311514_s1_p0/14667/18279_22278" or g.name == "m150911_220012_00127_c100861772550000001823190702121671_s1_p0/3853/0_6406" or g.name == "m150912_012316_00127_c100861772550000001823190702121672_s1_p0/131989/16299_25508" or g.name == "m150811_224417_00127_c100823122550000001823177111031570_s1_p0/157438/0_8938" or g.name == "m150910_220412_00127_c100822732550000001823176011031537_s1_p0/6739/7003_13868":
				#continue
			read=Read(g.name, 0, 0, 0) # create read for each read alignment
			prev_tmp=[]
			prev_locus= -1
			locus = -1
			for i in range(0,len(g.path.mapping)):
				if g.path.mapping[i].position.is_reverse != val1:
					val1 = False
					break
				else:
					count1 = count1 +1
					
			if count1 == len(g.path.mapping):
				count = count+1
				print(g.name)
				
			for i in range(0,len(g.path.mapping)):
				if g.path.mapping[i].position.is_reverse != val2:
					val2 = True
					break
				else:
					count2 = count2 +1
					
			if count2 == len(g.path.mapping):
				count = count+1
				print(g.name)
			if val1 ==val2:
				for i in range(0,len(g.path.mapping)-1):
				#for i in g.path.mapping: # go over the mapping in a read
					# TODO: check for forward or reverse strand, we may not need it for DAG.
					edge1 = tuple((g.path.mapping[i].position.node_id, g.path.mapping[i+1].position.node_id)) # go over nodes in a mapping
					edge2 = tuple((g.path.mapping[i+1].position.node_id, g.path.mapping[i].position.node_id)) # go over nodes in a mapping
					#print(edge)
					if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
						if edge1 in reverse_mapping:
							qualities = [10]* reverse_mapping[edge1][0][2]
							node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge1]] # consider (locus, branch)
						else:
							qualities = [10]* reverse_mapping[edge2][0][2]
							node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge2]]
						tmp = [x for x in node_inf]
						interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
						#if len(tmp)==1 and len(prev_tmp)==0: # simple bubble with multiple paths.
							#print('helloa')
							#qualities[tmp[0][1]] = 0
							#read.add_variant(tmp[0][0], tmp[0][1], qualities) # if any new locus of branch is encountered, enter a variant in read.
							#locus = tmp[0][0]
						if len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
							#print('hellob')
							qualities[interset_tmp[0][1]] = 0
							read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
							locus= interset_tmp[0][0]
						#elif len(prev_tmp) > 0 and len(tmp) ==1: # complicated bubbles, but one node can uniquely identify the branch.
							#print('helloc')
							#qualities[tmp[0][1]] = 0
							#read.add_variant(tmp[0][0], tmp[0][1], qualities)
							#locus = tmp[0][0]
							
						if prev_locus!=locus:
							prev_tmp = []
						else:
							for i in tmp:
								prev_tmp.append(i)
						prev_locus = locus


				if len(read) >= 2:
					readset.add(read)
	print("all forward")
	print(count)
	#print(readset)
	readset1=ReadSet()
	for read in readset:
		if read.sort() ==1:
			duplicated = duplicated +1
			continue
		else:
			readset1.add(read)

	readset1.sort()
	print("duplicated")
	print(duplicated)
	print("reads")
	print(len(readset1))
	return readset1, alleles_per_pos


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
	max_coverage = 20
	all_heterozygous = False
	distrust_genotypes = True
	with ExitStack() as stack:
		all_reads, alleles_per_pos = vg_reader(locus_file, gam_file)
		#print(all_reads)
		selected_indices = readselection(all_reads, max_coverage)
		selected_reads = all_reads.subset(selected_indices)
		print('positions from all reads')
		print(len(all_reads.get_positions()))
		print(len(selected_reads.get_positions()))

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
		#print(selected_reads)
		dp_table = PedigreeDPTable(selected_reads, recombination_costs, pedigree, distrust_genotypes, accessible_positions)
		superreads_list, transmission_vector = dp_table.get_super_reads()
		#positions = selected_reads.get_positions()
		positions = all_reads.get_positions()
		sorted(positions)
		vcf_indices = {position: index for index, position in enumerate(positions)}
		
		SNP_read_map = defaultdict(list)
		
		for index, read in enumerate(selected_reads):
			for variant in read:
				if variant.position in positions:
					snp_index = vcf_indices[variant.position]
					SNP_read_map[snp_index].append(index)
				else:
				      continue

		read_count = list()
		for a in range(0, len(positions)):
			readscommonatpos1_pos2= []
			if a < len(positions)-1 :
				pos1 = positions[a]
				pos2 = positions[a+1]
			x = vcf_indices[pos1]
			for i in SNP_read_map[x]:
				y = vcf_indices[pos2]
				for p in SNP_read_map[y]:
					if i==p:
						readscommonatpos1_pos2.append(i)
						readscommonatpos1_pos2.append(p)
						break
			read_count.append(len(set(readscommonatpos1_pos2)))
			  
		print("read_count")
		print(read_count)
		cost = dp_table.get_optimal_cost()
		#print(superreads_list[0])
		#print(cost)
		read_partitions = dp_table.get_optimal_partitioning()
		#print(read_partitions)
		
		out0 = open('pred_hap0.txt', 'w')
		out1 = open('pred_hap1.txt', 'w')
		
		for i,j in zip(read_partitions, selected_reads):
			if i ==0:
				out0.write(j.name + "\n")
			else:
				out1.write(j.name + "\n")


def add_arguments(parser):
	arg = parser.add_argument
	# Positional arguments
	arg('locus_file', metavar='LOCUS', help='variants in LOCUS file to phase')
	arg('gam_file', metavar='PHASEINPUT', help='read alignments in GAM file ')
	#TODO: add k parameter

def main(args):
	run_phaseg(**vars(args))

