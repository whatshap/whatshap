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

def subset_clusters(readset, clustering,ploidy, sample, genotypes, single_block, cpp_threading, dynamic_switchcost):
	logger.debug("Computing cluster coverages and consensi ..")
	
	# Map genome positions to [0,l)
	index, rev_index = get_position_map(readset)
	num_vars = len(rev_index)
	num_clusters = len(clustering)
	
	# for each position, compute the relevant list of clusters (sorted by descendant coverage)
	cov_map = get_pos_to_clusters_map(readset, clustering, index, ploidy)

	#create dictionary mapping a clusterID to a pair (starting position, ending position)
	positions = get_cluster_start_end_positions(readset, clustering, index)

	#for every cluster and every variant position, compute the relative coverage	
	coverage = get_coverage(readset, clustering, index)
	
	#compute the consensus sequences for every variant position and every cluster for integrating genotypes in the DP
	consensus = get_cluster_consensus_local(readset, clustering, cov_map, positions)

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
		#print("No cluster with fitting genotype: ", counter)
		
		#perform the dynamic programming to fill the scoring matrix (in Cython to speed up computation)
		scoring = clustering_DP(num_vars,clustering,coverage, positions, cov_map, ploidy, genotypes, consensus, geno_map)

		#start the backtracing
		path = []	
		#find the last column that contains a minimum other than INT_MAX (1000000, respectively) 
		start_col = find_backtracing_start(scoring, num_vars)
		#print("starting in column: ", start_col)
		last_min_idx = min((val, idx) for (idx, val) in enumerate([i[0] for i in scoring[start_col]]))[1]
		#print("DP Score = "+str(min((val, idx) for (idx, val) in enumerate([i[0] for i in scoring[start_col]]))[0]))

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
		logger.debug("Cut positions: ", cut_positions)
	else:
		#block_starts = compute_linkage_based_block_starts(readset, index, ploidy)
		block_starts = [0]
		
		# arrange data
		compressed_coverage = [] #rearrange the content of "coverage" in a position-wise way, such that the i-th coverage refers to the i-th cluster in cov_map
		cov_map_as_list = []
		compressed_consensus = [] #for every position, give a list of consensi, such that the i-th consensus refers to the i-th cluster in cov_map
		for pos in range(num_vars):
			coverage_list = []
			consensus_list = []
			for i in range(len(cov_map[pos])):
				coverage_list.append(coverage[pos][cov_map[pos][i]])
				consensus_list.append(consensus[pos][cov_map[pos][i]])
			compressed_coverage.append(coverage_list)
			cov_map_as_list.append(cov_map[pos])
			compressed_consensus.append(consensus_list)
			
		# compute cluster dissimilarity on consensus
		coverage_abs = get_coverage_absolute(readset, clustering, index)
		c_to_c_dissim = [[[0 for j in range(len(cov_map[p+1]))] for i in range(len(cov_map[p]))] for p in range(num_vars-1)]
		
		cluster_zeroes = [dict() for c_id in range(num_clusters)]
		cluster_ones = [dict() for c_id in range(num_clusters)]
		for pos in range(num_vars):
			for c_id in consensus[pos]:
				cluster_zeroes[c_id][pos] = coverage_abs[pos][c_id] * (1-consensus[pos][c_id])
				cluster_ones[c_id][pos] = coverage_abs[pos][c_id] * consensus[pos][c_id]
				
		avg_dissim = 0
		count = 0
		for var in range(num_vars-1):
			for i, c1 in enumerate(cov_map[var]):
				for j, c2 in enumerate(cov_map[var+1]):
					if dynamic_switchcost:
						same = 0
						diff = 0
						for pos in range(max(0, var-15), min(num_vars-1, var+15)):
							if pos in cluster_zeroes[c1] and pos in cluster_zeroes[c2]:
								same += cluster_zeroes[c1][pos] * cluster_zeroes[c2][pos] + cluster_ones[c1][pos] * cluster_ones[c2][pos]
								diff += cluster_zeroes[c1][pos] * cluster_ones[c2][pos] + cluster_ones[c1][pos] * cluster_zeroes[c2][pos]
						c_to_c_dissim[var][i][j] = diff / (same+diff) if diff > 0 else 0
						if c1 != c2:
							avg_dissim += c_to_c_dissim[var][i][j]
							count += 1
					else:
						c_to_c_dissim[var][i][j] = 1 if c1 != c2 else 0
						
		avg_dissim = avg_dissim / count if count > 0 else 1.0
			
		# run threader
		threader = HaploThreader(ploidy, 32.0/avg_dissim, 8.0, True, 0)
		path = threader.computePathsBlockwise(block_starts, cov_map_as_list, compressed_coverage, compressed_consensus, genotypes, c_to_c_dissim)
			
		assert(len(path) == num_vars)
	
		#determine cut positions: Currently, a cut position is created every time a path changes the used cluster from one position to the next
		cut_positions = [0]
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
		logger.debug("Cut positions: ", cut_positions)
	
	# compute haplotypes
	haplotypes = []
	for i in range(ploidy):
		hap = ""
		alleles_as_strings = []
		for pos in range(len(path)):
			c_id = path[pos][i]
			allele = consensus[pos][c_id] if c_id in consensus[pos] else -1
			if allele == -1:
				alleles_as_strings.append("n")
			else:
				alleles_as_strings.append(str(allele))
		haplotypes.append(hap.join(alleles_as_strings))
		
	return (cut_positions, path, haplotypes)

def get_position_map(readset):
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0

	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1
		
	return index, rev_index

def get_pos_to_clusters_map(readset, clustering, pos_index, ploidy):
	#cov_map maps each position to the list of cluster IDS that appear at this position
	num_clusters = len(clustering)
	cov_map = defaultdict(list)
	for c_id in range(num_clusters):
		covered_positions = []
		for read in clustering[c_id]:
			for pos in [pos_index[var.position] for var in readset[read] ]:		
				if pos not in covered_positions:
					covered_positions.append(pos)
		for p in covered_positions:
			cov_map[p].append(c_id)

	#restrict the number of clusters at each position to a maximum of 2*ploidy clusters with the largest number of reads	
	for key in cov_map.keys():
		largest_clusters = sorted(cov_map[key], key=lambda x: len(clustering[x]), reverse=True)[:min(len(cov_map[key]), 2*ploidy)]
		cov_map[key] = largest_clusters
	return cov_map

def get_coverage(readset, clustering, pos_index):
	num_vars = len(pos_index)
	num_clusters = len(clustering)
	coverage = [dict() for pos in range(num_vars)]
	coverage_sum = [0 for pos in range(num_vars)]
	for c_id in range(num_clusters):
		for read in clustering[c_id]:
			for pos in [pos_index[var.position] for var in readset[read]]:
				if c_id not in coverage[pos]:
					coverage[pos][c_id] = 0
				coverage[pos][c_id] += 1
				coverage_sum[pos] += 1

	for pos in range(num_vars):
		for c_id in coverage[pos]:
			coverage[pos][c_id] = coverage[pos][c_id]/coverage_sum[pos]
	
	return coverage

def get_coverage_absolute(readset, clustering, pos_index):
	num_vars = len(pos_index)
	num_clusters = len(clustering)
	coverage = [dict() for pos in range(num_vars)]
	for c_id in range(num_clusters):
		for read in clustering[c_id]:
			for pos in [pos_index[var.position] for var in readset[read]]:
				if c_id not in coverage[pos]:
					coverage[pos][c_id] = 0
				coverage[pos][c_id] += 1
	
	return coverage

def get_cluster_start_end_positions(readset, clustering, pos_index):
	num_clusters = len(clustering)
	positions = {}
	counter = 0
	for c_id in range(num_clusters):
		read = clustering[c_id][0]
		start = pos_index[readset[read][0].position]
		end = pos_index[readset[read][-1].position]
		for read in clustering[c_id]:
			readstart = pos_index[readset[read][0].position]
			readend = pos_index[readset[read][-1].position]
			if (readstart < start):
				start = readstart
			if (readend > end):
				end = readend
		positions[c_id] = (start,end)
	assert(len(positions) == num_clusters)
	return positions

def get_cluster_consensus_local(readset, clustering, cov_map, positions, fractional = False):
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
	if fractional:
		clusterwise_consensus = [get_single_cluster_consensus_frac(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	else:
		clusterwise_consensus = [get_single_cluster_consensus(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	whole_consensus = []
	for pos in range(num_vars):
		newdict = defaultdict()
		for c in cov_map[pos]:
			newdict[c] = clusterwise_consensus[c][pos]
		whole_consensus.append(newdict)	
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

def compute_linkage_based_block_starts(readset, pos_index, ploidy):
	num_vars = len(pos_index)

	# cut threshold
	min_block_size = 1
	if ploidy == 2:
		cut_threshold = 1
	else:
		cut_threshold = ploidy*ploidy
		for i in range(ploidy-1, ploidy*ploidy):
			# chance to cover at most ploidy-2 haplotypes with i reads
			cut_threshold = i
			if ploidy * pow((ploidy-2)/ploidy, i) < 0.02:
				cut_threshold = i
				break
	logger.debug("Cut position threshold: coverage >= {}".format(cut_threshold))

	# compute absolute coverage for every position and the number of links between each pair of positions
	pos_coverage = [0]*num_vars
	link_coverage = [dict() for i in range(num_vars)]
	for read in readset:
		for pos in [pos_index[var.position] for var in read]:
			pos_coverage[pos] += 1
			for pos2 in [pos_index[var.position] for var in read]:
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
	
	block_list = list(block_starts)
	block_list.sort()
	return block_list