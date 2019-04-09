#from copy import deepcopy
#from math import floor, ceil
import itertools as it
import logging
from collections import defaultdict
from .core import clustering_DP, ReadSet, Read, HaploThreader
#import numpy as np
#from .testhelpers import string_to_readset
#from .phase import find_components
#from statistics import mean
#import operator
import random
from queue import *

logger = logging.getLogger(__name__)

def calc_consensus_blocks(readset, clustering, cluster_blocks, cut_positions, consensus):
	#cluster_consensus = get_cluster_consensus(readset, clustering)

	consensus_blocks = cluster_blocks	
	for i, block in enumerate(consensus_blocks):
		offset = cut_positions[i-1] if i > 0 else 0	
		for j, hap in enumerate(block):
			for pos, clust in enumerate(hap):
				consensus_blocks[i][j][pos] = consensus[offset+pos][clust] if clust != None else -1
				#consensus_blocks[i][j][pos] = consensus[clust][offset+pos] if clust != None else -1
				
	return consensus_blocks

def get_cluster_consensus(readset, clustering):
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1
		
	# Create consensus matrix
	cluster_consensus = [[{} for j in range(num_vars)] for i in range(len(clustering))]
	
	# Count all alleles
	for c in range(len(clustering)):
		for read in clustering[c]:
			alleles = [(index[var.position], var.allele) for var in readset[read]]
			for (pos, allele) in alleles:
				if allele not in cluster_consensus[c][pos]:
					cluster_consensus[c][pos][allele] = 0
				cluster_consensus[c][pos][allele] += 1
				
	# Compute consensus for all positions and clusters
	for c in range(len(clustering)):
		for pos in range(num_vars):
			alleles = cluster_consensus[c][pos]
			max_count = -1
			max_allele = -1
			for key in alleles:
				if alleles[key] > max_count:
					max_count = alleles[key]
					max_allele = key
			cluster_consensus[c][pos] = max_allele		
	return cluster_consensus


def subset_clusters(readset, clustering,ploidy, sample, genotypes, single_block, cpp_threading):
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1

	num_clusters = len(clustering)
	#print('number of variants: ', num_vars)
	#print('number of clusters: ', num_clusters)
	
	#cov_map maps each position to the list of cluster IDS that appear at this position
	cov_positions = defaultdict()
	cov_map = defaultdict(list)
	for c_id in range(num_clusters):
		covered_positions = []
		for read in clustering[c_id]:
			for pos in [index[var.position] for var in readset[read] ]:		
				if pos not in covered_positions:
					covered_positions.append(pos)
		for p in covered_positions:
			cov_map[p].append(c_id)
#	assert(len(cov_map.keys()) == num_vars)
	print("map of clusters at every position computed")
	
	#compute for each position the amount of clusters that 'cover' this position
	for key in cov_map.keys():
		cov_positions[key] = len(cov_map[key])
	cov_positions_sorted = sorted(cov_positions.items(), key=lambda x: x[1], reverse=True)
	
	#restrict the number of clusters at each position to a maximum of 2*ploidy clusters with the largest number of reads	
	for key in cov_map.keys():
		largest_clusters = sorted(cov_map[key], key=lambda x: len(clustering[x]), reverse=True)[:min(len(cov_map[key]), 2*ploidy)]
		cov_map[key] = largest_clusters

	#create dictionary mapping a clusterID to a pair (starting position, ending position)
	positions = {}
	counter = 0
	for c_id in range(num_clusters):
		read = clustering[c_id][0]
		start = index[readset[read][0].position]
		end = index[readset[read][-1].position]
		for read in clustering[c_id]:
			readstart = index[readset[read][0].position]
			readend = index[readset[read][-1].position]
			if (readstart < start):
				start = readstart
			if (readend > end):
				end = readend
		positions[c_id] = (start,end)
	assert(len(positions) == num_clusters)
	print("map of clusters to their read positions computed")

	#for every cluster and every variant position, compute the relative coverage
	coverage = [[0]*num_vars for i in range(num_clusters)]
	for c_id in range(num_clusters):
		for read in clustering[c_id]:
			for pos in [index[var.position] for var in readset[read]]:			
				coverage[c_id][pos] += 1
	print("total coverage computed")
	
	coverage_sum = [sum([coverage[i][j] for i in range(num_clusters)]) for j in range(num_vars)]
	
	for c_id in range(num_clusters):
		coverage[c_id] = [((coverage[c_id][i])/coverage_sum[i]) if coverage_sum[i]!= 0 else 0 for i in range(num_vars)]
	assert(len(coverage) == num_clusters)
	assert(len(coverage[0]) == num_vars)
	print("relative coverage computed")
	
	#compute the consensus sequences for every variant position and every cluster for integrating genotypes in the DP
#	consensus = get_cluster_consensus(readset, clustering)
	consensus = get_cluster_consensus_local(readset, clustering, cov_map, positions)
	print("consensus sequences computed")

	if not cpp_threading:
		geno_map = defaultdict(list)
		counter = 0
		for pos in range(num_vars):
			c_ids = sorted(cov_map[pos])
			c_tuples = sorted(list(it.combinations_with_replacement(c_ids, ploidy)))
			geno_tuples = [tup for tup in c_tuples if (compute_tuple_genotype(consensus,tup, pos) == genotypes[pos])]
			if (len(geno_tuples) == 0):
				#TODO add option to use the next best genotypes if also the next list is empty
				geno_tuples = [tup for tup in c_tuples if (compute_tuple_genotype_soft(consensus,tup, pos, genotypes[pos]) == 1)]
				if (len(geno_tuples) == 0):
					geno_tuples = c_tuples
			geno_map[pos] = geno_tuples
		print("No cluster with fitting genotype: ", counter)
		
		#perform the dynamic programming to fill the scoring matrix (in Cython to speed up computation)
		scoring = clustering_DP(num_vars,clustering,coverage, positions, cov_map, ploidy, genotypes, consensus, geno_map)

		#start the backtracing
		path = []	
		#find the last column that contains a minimum other than INT_MAX (1000000, respectively) 
		start_col = find_backtracing_start(scoring, num_vars)
		print("starting in column: ", start_col)
		last_min_idx = min((val, idx) for (idx, val) in enumerate([i[0] for i in scoring[start_col]]))[1]
		print("DP Score = "+str(min((val, idx) for (idx, val) in enumerate([i[0] for i in scoring[start_col]]))[0]))

		#instead of using all tuples, compute only the tuples relevant for the position start_col
		clusters_tuples = geno_map[start_col]
		conf_clusters_tuples = list(it.chain.from_iterable(sorted(list(set(it.permutations(x)))) for x in clusters_tuples))

		path.append(conf_clusters_tuples[last_min_idx])
		#append stored predecessor
		for i in range(start_col,0,-1):
			pred = scoring[i][last_min_idx][1]
			#compute variants relevant for position pred
			tups = geno_map[i-1]	
			conf_tups = list(it.chain.from_iterable(sorted(list(set(it.permutations(x)))) for x in tups))

			path.append(conf_tups[pred])
			last_min_idx = pred
		#path has been assembled from the last column and needs to be reversed
		path.reverse()
		assert(len(path) == start_col+1)
		
		#determine cut positions: Currently, a cut position is created every time a path changes the used cluster from one position to the next
		if single_block:
			cut_positions = [len(readset.get_positions())-1]
		else:
			cut_positions = []
			for i in range(len(path)-1):
				dissim = 0
				for j in range(0,ploidy):
					#if (i < len(path) -2):
		#				if path[i][j] != path[i+1][j] and path[i][j] != path[i+2][j]:
		#					dissim +=1
		#			else:
		#				if path[i][j] != path[i+1][j]:
		#					dissim += 1
					if path[i][j] != path[i+1][j]:
						dissim += 1
				if (dissim >= 1):
					cut_positions.append(i)
		print("cut positions: ", cut_positions)
	else:
		print("Precomputing cut positions based on low coverage and linkage")
		# cut threshold
		cut_threshold = 8 #TODO: Compute depending on ploidy
		min_block_size = 2
		
		# compute absolute coverage for every position and the number of links between each pair of positions
		pos_coverage = [0]*num_vars
		link_coverage = [dict() for i in range(num_vars)]
		for c_id in range(num_clusters):
			for read in clustering[c_id]:
				for pos in [index[var.position] for var in readset[read] ]:
					pos_coverage[pos] += 1
					for pos2 in [index[var.position] for var in readset[read] ]:
						if pos2 not in link_coverage[pos]:
							link_coverage[pos][pos2] = 0
						link_coverage[pos][pos2] += 1
		
		# form connected sets
		pos_clust = [0]*num_vars
		clust_count = 0
		for pos in range(num_vars):
			current_cluster = 0
			if pos_clust[pos] > 0:
				continue
			else:
				clust_count += 1
				current_cluster = clust_count
				pos_clust[pos] = current_cluster
			pending = Queue()
			pending.put(pos)
			while not pending.empty():
				cur = pending.get()
				for linked in sorted(link_coverage[cur]):
					if link_coverage[cur][linked] >= cut_threshold and pos_clust[linked] == 0:
						pos_clust[linked] = current_cluster
						pending.put(linked)
		
		# make intervals of connected positions
		current_cluster = 0
		current_begin = 0
		intervals = []
		block_starts = set()
		for pos in range(num_vars):
			#print("Pos "+str(pos)+": cov="+str(pos_coverage[pos])+" clust="+str(pos_clust[pos])+" linklast="+str(link_coverage[pos][pos-1] if (pos-1) in link_coverage[pos] else 0))
			if pos_coverage[pos] < cut_threshold:
				if pos - current_begin >= min_block_size:
					intervals.append((current_begin, pos))
				else:
					for i in range(pos-1, current_begin, -1):
						block_starts.add(i)
				current_begin = pos + 1
				current_cluster = 0
				block_starts.add(pos)
			else:
				if current_cluster != pos_clust[pos]:
					if pos - current_begin >= min_block_size:
						intervals.append((current_begin, pos))
					else:
						for i in range(pos-1, current_begin, -1):
							block_starts.add(i)
					current_begin = pos
					current_cluster = pos_clust[pos]
					block_starts.add(pos)
		if num_vars - current_begin >= min_block_size:
			intervals.append((current_begin, num_vars))
		else:
			for i in range(num_vars-1, current_begin, -1):
				block_starts.add(i)
		#print("block starts positions: ", block_starts)
		
		# arrange data
		position_wise_coverage = [] #rearrange the content of "coverage" in a position-wise way, such that the i-th coverage refers to the i-th cluster in cov_map
		cov_map_as_list = []
		compressed_consensus = [] #for every position, give a list of consensi, such that the i-th consensus refers to the i-th cluster in cov_map
		for pos in range(num_vars):
			coverage_list = []
			consensus_list = []
			for i in range(len(cov_map[pos])):
				coverage_list.append(coverage[cov_map[pos][i]][pos])
				consensus_list.append(consensus[pos][cov_map[pos][i]])
			position_wise_coverage.append(coverage_list)
			cov_map_as_list.append(cov_map[pos])
			compressed_consensus.append(consensus_list)
			
		# run threader
		threader = HaploThreader(ploidy, 32.0, 8.0, True, 0)
		path = threader.computePaths(list(sorted(list(block_starts))), cov_map_as_list, position_wise_coverage, compressed_consensus, genotypes)
			
		assert(len(path) == num_vars)
		
		# remove single variant jumps from the path
		#for i in range(1, len(path)-1):
		#	if i in block_starts:
		#		continue
		#	dissim = 0
		#	for j in range(0,ploidy):
		#		if path[i-1][j] != path[i+1][j]:
		#			dissim += 1
		#	if dissim == 0:
		#		for j in range(0,ploidy):
		#			if path[i-1][j] in consensus[i]:
		#				path[i][j] = path[i-1][j]
					
		# if wrong genotype, correct by changing multi-cluster-entries
		#corrections = dict()
		#
		#fconsensus = get_cluster_consensus_local(readset, clustering, cov_map, positions, fractional = True)
		#for i in range(0, len(path)):
		#	phased_genotype = sum([consensus[i][path[i][j]] if path[i][j] in consensus[i] else -999 for j in range(ploidy)])
		#	geno_diff = genotypes[i] - phased_genotype
		#	if abs(geno_diff) == 1:
		#		# count multiplicity of all present clusters
		#		count = dict()
		#		present = set()
		#		max_count = 0
		#		for j in range(ploidy):
		#			if path[i][j] not in count:
		#				count[path[i][j]] = 0
		#				present.add(path[i][j])
		#			count[path[i][j]] += 1
		#			max_count = max(max_count, count[path[i][j]])
		#			
		#		if (max_count < 2):
		#			# only single-hap clusters, nothing we can do here
		#			continue
		#			
		#		# determine highest deviation from fractional consensus in desired direction
		#		max_cid = -1
		#		max_deviation = -float("inf")
		#		for cid in present:
		#			deviation = (fconsensus[i][cid] - consensus[i][cid]) * count[cid] * (geno_diff)
		#			if deviation > max_deviation:
		#				max_deviation = deviation
		#				max_cid = cid
		#				
		#		if max_cid not in corrections:
		#			corrections[max_cid] = []
		#		
		#		corrections[max_cid].append((i, genotypes[i] - phased_genotype))
		#		print("Correction needed: Pos="+str(i)+" cluster="+str(max_cid)+" direction="+str(geno_diff))
					
		
	
		#determine cut positions: Currently, a cut position is created every time a path changes the used cluster from one position to the next
		cut_positions = []
		for i in range(1, len(path)):
			if i in block_starts:
				cut_positions.append(i)
			else:
				dissim = 0
				for j in range(0,ploidy):
					if path[i][j] != path[i-1][j]:
						dissim += 1
				if (not single_block and dissim >= 2):
					cut_positions.append(i)
		print("cut positions: ", cut_positions)
	
	#experimental: Computing longer blocks and discarding cut positions close to each other
#	cuts = []
#	for c in range(len(cut_positions)-2):
#		dissim = 0
#		current = cut_positions[c]
#		next = cut_positions[c+1]
#		prev = cut_positions[c-1]		
#		for j in range(0,ploidy):
#			if path[current][j] != path[next][j]:
#				dissim += 1
#		if (dissim == 1):
#			dist_bwd = current - prev
#			dist_fwd = cut_positions[c+2] - next
#			if dist_bwd < 4 or dist_fwd < 4:
#				cuts.append(cut_positions[c])
#	cut_positions = cuts[:]
#	print("cut positions after: ", cut_positions)

	#divide path of clusters into <ploidy> separate paths
	haps = []
	for i in range(ploidy):
		temp = [tup[i] for tup in path]
		haps.append(temp)
	assert(len(haps) == ploidy)

	#for testing purposes: compute paths by randomly connecting clusters
#	haps = [[] for i in range(ploidy)]
#	for pos in range(num_vars):
#		possible_clusters = cov_map[pos]
#		for i in range(ploidy):
#			random_cluster = random.choice(possible_clusters)
#			haps[i].append(random_cluster)
#	assert(len(haps[0]) == num_vars)

	# Return cluster blocks: List of blocks, each blocks ranges from one cut position to the next and contains <ploidy> sequences
	# of cluster ids that indicate which cluster goes to which haplotype at which position.
	cluster_blocks = []
	old_pos = 0
	for cut_pos in cut_positions:
		cluster_blocks.append([hap[old_pos:cut_pos] for hap in haps])
		old_pos = cut_pos
	if (len(cut_positions) == 0):
		last_cut = 0
	else:	
		last_cut = cut_positions[len(cut_positions)-1]
	cluster_blocks.append([hap[last_cut:] for hap in haps])
	assert(len(cluster_blocks) == len(cut_positions)+1)

	consensus_blocks = calc_consensus_blocks(readset, clustering, cluster_blocks, cut_positions, consensus)

	haplotypes = []
	for i in range(ploidy):
		hap = ""
		alleles_as_strings = []
		for block in consensus_blocks:
			for allele in block[i]:
				if allele == -1:
					alleles_as_strings.append("n")
				else:
					alleles_as_strings.append(str(allele))
		hap = hap.join(alleles_as_strings)
		haplotypes.append(hap)	

	#write new VCF file	
	superreads, components = dict(), dict()
	
	accessible_positions = sorted(readset.get_positions())
#	accessible_positions = local_positions
	cut_positions = []
	overall_components = {}
	last_cut = 0
	for haploblock in consensus_blocks[:len(consensus_blocks)]:
		next_cut = last_cut + len(haploblock[0])
		#print("last_cut = "+str(last_cut)+" next_cut = "+str(next_cut))
		cut_positions.append(accessible_positions[last_cut])
		
		for pos in range(last_cut, next_cut):
			overall_components[accessible_positions[pos]] = accessible_positions[last_cut]
			overall_components[accessible_positions[pos]+1] = accessible_positions[last_cut]
		last_cut = next_cut
	components[sample] = overall_components
	readset = ReadSet()
	for i in range(ploidy):
		read = Read('superread {}'.format(i+1), 0, 0)
		# insert alleles
		for j,allele in enumerate(haplotypes[i]):
			if (allele=="n"):
				continue
			allele = int(allele)
			qual = [10,10]
			qual[allele] = 0
			read.add_variant(accessible_positions[j], allele, qual)
		readset.add(read)

	superreads[sample] = readset

	return(cut_positions, cluster_blocks, components, superreads, coverage, path, haplotypes)

def get_cluster_consensus_local(readset, clustering, cov_map, positions, fractional = False):
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1
		
#	cluster_consensus = []
#	for c in range(len(clustering)):
#		row = defaultdict(list)
#		consensus = defaultdict(int)
#		for read in clustering[c]:
#			for (pos, allele) in [(index[var.position], var.allele) for var in readset[read]]:
#				
#				row[pos][0] += allele
#				row[pos][1] += 1
#		for key in row.keys():
#			max_allele = row[key][0]/row[key][1]
#			if max_allele < 0.5:
#				consensus[key] = 0
#			else:
#				consensus[key] = 1
#		del row
#		cluster_consensus.append(consensus)

#	cluster_consensus = [[{} for j in range(num_vars)] for i in range(len(clustering))]						
#
#	for c in range(len(clustering)):
#		for read in clustering[c]:
#			alleles = [(index[var.position], var.allele) for var in readset[read]]
#			for (pos, allele) in alleles:
#				if allele not in cluster_consensus[c][pos]:
#					cluster_consensus[c][pos][allele] = 0
#				cluster_consensus[c][pos][allele] += 1
#				
#	# Compute consensus for all positions and clusters
##	for c in range(len(clustering)):
##		#for pos in range(positions[c][0], positions[c][1]+1):
##		for read in clustering[c]:
##			for pos in [index[var.position] for var in readset[read]]:
#	posi_map = defaultdict(list)
#	for c in range(len(clustering)):
#		posis = []
#		for read in clustering[c]:
#			positions = [index[var.position] for var in readset[read]]
#			posis.extend(positions)	
#		pos =sorted(list(set(posis)))
#		posi_map[c] = pos
#		
#	for c in range(len(clustering)):
#		for pos in posi_map[c]:
#			alleles = cluster_consensus[c][pos]
#			max_count = -1
#			max_allele = -1
#			for key in alleles:
#				if alleles[key] > max_count:
#					max_count = alleles[key]
#					max_allele = key
#			cluster_consensus[c][pos] = max_allele

	relevant_pos = [[] for i in range(len(clustering))]
	for pos in range(num_vars):
		for c in cov_map[pos]:
			relevant_pos[c].append(pos)

	if fractional:
		clusterwise_consensus = [get_single_cluster_consensus_frac(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	else:
		clusterwise_consensus = [get_single_cluster_consensus(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	whole_consensus = []
	for pos in range(num_vars):
		#if (pos%1000 == 0):
		#	print("computing consensus for position: ", pos)
		newdict = defaultdict()
		for c in cov_map[pos]:
			newdict[c] = clusterwise_consensus[c][pos]
		whole_consensus.append(newdict)

	#cluster_consensus = []
	#for pos in range(num_vars):
	#	if (pos%1000 == 0):
	#		print("computing consensus for position: ", pos)
	#	newdict = defaultdict()
	#	for c in cov_map[pos]:
	#		allele0, allele1 = 0,0
	#		for read in clustering[c]:
	#			alleles = [(index[var.position], var.allele) for var in readset[read]]
	#			#todo: get read[position].allele directly without iterating through every var in read?
	#			for (posi, allele) in alleles:	
	#				if (posi == pos):				
	#					if allele == 1:
	#						allele1 += 1
	#					else:
	#						allele0 += 1			
	#		
	#		max_allele = 0
	#		if (allele1 > allele0):	
	#			max_allele = 1	
	#		else:
	#			max_allele = 0
	#		newdict[c] = max_allele
	#	cluster_consensus.append(newdict)
		
	#for pos in range(num_vars):
	#	for c in cov_map[pos]:
	#		assert(cluster_consensus[pos][c] == whole_consensus[pos][c])
	
	return whole_consensus

def get_single_cluster_consensus(readset, cluster, index, relevant_pos):
	# Count zeroes and one for every position
	num_zero = {}
	num_one = {}
	for read in cluster:
		for var in readset[read]:
			pos = index[var.position]
			if pos not in num_zero:
				num_zero[pos] = 0
				num_one[pos] = 0
			if var.allele == 0:
				num_zero[pos] += 1
			else:
				num_one[pos] += 1
				
	# Determine majority allele for every position
	cluster_consensus = {}
	for pos in relevant_pos:
		if pos in num_zero:
			if num_one[pos] > num_zero[pos]:
				cluster_consensus[pos] = 1
			else:
				cluster_consensus[pos] = 0
		else:
			cluster_consensus[pos] = 0
	return cluster_consensus

def get_single_cluster_consensus_frac(readset, cluster, index, relevant_pos):
	# Count zeroes and one for every position
	num_zero = {}
	num_one = {}
	for read in cluster:
		for var in readset[read]:
			pos = index[var.position]
			if pos not in num_zero:
				num_zero[pos] = 0
				num_one[pos] = 0
			if var.allele == 0:
				num_zero[pos] += 1
			else:
				num_one[pos] += 1
				
	# Determine majority allele for every position
	cluster_consensus = {}
	for pos in relevant_pos:
		if pos in num_zero:
			cluster_consensus[pos] = float(num_one[pos]) / float(num_one[pos] + num_zero[pos])
		else:
			cluster_consensus[pos] = 0.5
	return cluster_consensus

def get_cluster_consensusfrac_local(readset, clustering, cov_map, positions):
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1

	relevant_pos = [[] for i in range(len(clustering))]
	for pos in range(num_vars):
		for c in cov_map[pos]:
			relevant_pos[c].append(pos)

	clusterwise_consensus = [get_single_cluster_consensus(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	whole_consensus = []
	for pos in range(num_vars):
		newdict = defaultdict()
		for c in cov_map[pos]:
			newdict[c] = clusterwise_consensus[c][pos]
		whole_consensus.append(newdict)
	
	return whole_consensus

def compute_tuple_genotype(consensus,tup, var):
	genotype = 0
	for i in tup:
		#allele = consensus[i][var]
		allele = consensus[var][i]
		genotype += allele
	return(genotype)

def compute_tuple_genotype_soft(consensus,tup, var, geno):
	genotype = 0
	for i in tup:
		allele = consensus[var][i]
		#allele = consensus[i][var]
		genotype += allele
	res = max((geno-genotype),(genotype-geno))
	return(res)

def find_backtracing_start(scoring, num_vars):
	minimum = 1000000
	last_col = num_vars-1
	res = False
	for i in scoring[last_col]:
		if i[0] < minimum:
			res=True
	if res:
		return(last_col)
	else:
		return(find_backtracing_start(scoring,last_col-1))

