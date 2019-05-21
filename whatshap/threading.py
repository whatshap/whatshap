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

def subset_clusters(readset, clustering, ploidy, sample, genotypes, block_cut_sensitivity, cython_threading):
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
	consensus = get_local_cluster_consensus(readset, clustering, cov_map, positions)

	if cython_threading:
		geno_map = defaultdict(list)
		counter = 0
		for pos in range(num_vars):
			c_ids = sorted(cov_map[pos])
			c_tuples = sorted(list(it.combinations_with_replacement(c_ids, ploidy)))
			geno_tuples = [tup for tup in c_tuples if compute_tuple_genotype_dist(consensus,tup, pos, genotypes[pos]) == 0]
			if (len(geno_tuples) == 0):
				#TODO add option to use the next best genotypes if also the next list is empty
				geno_tuples = [tup for tup in c_tuples if compute_tuple_genotype_dist(consensus,tup, pos, genotypes[pos]) == 1]
				if (len(geno_tuples) == 0):
					geno_tuples = c_tuples
			geno_map[pos] = geno_tuples
		#print("No cluster with fitting genotype: ", counter)
		
		#perform the dynamic programming to fill the scoring matrix (in Cython to speed up computation)
		scoring = clustering_DP(num_vars,clustering,coverage, positions, cov_map, ploidy, consensus, geno_map)

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

	else:
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
				
		count = 0
		for var in range(num_vars-1):
			for i, c1 in enumerate(cov_map[var]):
				for j, c2 in enumerate(cov_map[var+1]):
					same = 0
					diff = 0
					for pos in range(max(0, var-10), min(num_vars-1, var+10)):
						if pos in cluster_zeroes[c1] and pos in cluster_zeroes[c2]:
							same += cluster_zeroes[c1][pos] * cluster_zeroes[c2][pos] + cluster_ones[c1][pos] * cluster_ones[c2][pos]
							diff += cluster_zeroes[c1][pos] * cluster_ones[c2][pos] + cluster_ones[c1][pos] * cluster_zeroes[c2][pos]
					c_to_c_dissim[var][i][j] = diff / (same+diff) if diff > 0 else 0
					if c1 != c2:
						count += 1
			
		# run threader
		threader = HaploThreader(ploidy, 48.0, 8.0, True, 0)
		path = threader.computePathsBlockwise([0], cov_map_as_list, compressed_coverage, compressed_consensus, genotypes, c_to_c_dissim)
		assert(len(path) == num_vars)
		
	# end threading
	
	# determine cut positions
	cut_positions = [0]
	if block_cut_sensitivity >= 3:
		if block_cut_sensitivity >= 6:
			dissim_threshold = 1
			rise_fall_dissim = ploidy+1
		elif block_cut_sensitivity == 5:
			dissim_threshold = 2
			rise_fall_dissim = ploidy+1
		elif block_cut_sensitivity == 4:
			dissim_threshold = ploidy+1
			rise_fall_dissim = ploidy+1
		else:
			dissim_threshold = 2
			rise_fall_dissim = 0
		
		copynrs = []
		for i in range(0, len(path)):
			copynr = dict()
			for j in range(0, ploidy):
				if path[i][j] not in copynr:
					copynr[path[i][j]] = 0
				copynr[path[i][j]] += 1
			copynrs.append(copynr)
			
		cpn_rising = [False for c_id in range(num_clusters)]
		cpn_falling = [False for c_id in range(num_clusters)]
		
		for i in range(1, len(path)):
			dissim = 0
			for j in range(0, ploidy):
				old_c = path[i-1][j]
				new_c = path[i][j]
				if old_c != new_c:
					rise_fall = False
					# check if previous cluster went down from copy number >= 2 to a smaller one >= 1
					if old_c in copynrs[i] and copynrs[i][old_c] >= 1:
						cpn_falling[old_c] = True
						if cpn_rising[old_c]:
							rise_fall = True
					# check if new cluster went up from copy number >= 1 to a greater one >= 2
					if new_c in copynrs[i-1] and copynrs[i-1][new_c] >= 1:
						cpn_rising[new_c] = True
					# check if one cluster has been rising and then falling in the current block
					if rise_fall:
						dissim += rise_fall_dissim
						
					# count general switches
					dissim += 1
			if dissim >= dissim_threshold:
				cut_positions.append(i)
				cpn_rising = [False for c_id in range(num_clusters)]
				cpn_falling = [False for c_id in range(num_clusters)]
				
	logger.debug("Cut positions: {}".format(cut_positions))
	
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

def get_local_cluster_consensus(readset, clustering, cov_map, positions):
	return [{c_id : pos_cons[c_id][0] for c_id in pos_cons} for pos_cons in get_local_cluster_consensus_withfrac(readset, clustering, cov_map, positions)]

def get_local_cluster_consensus_withfrac(readset, clustering, cov_map, positions):
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
	
	clusterwise_consensus = [get_single_cluster_consensus_frac(readset, clustering[i], index, relevant_pos[i]) for i in range(len(clustering))]
	whole_consensus = []
	for pos in range(num_vars):
		newdict = defaultdict()
		for c in cov_map[pos]:
			newdict[c] = clusterwise_consensus[c][pos]
		whole_consensus.append(newdict)	
	return whole_consensus

def get_single_cluster_consensus_frac(readset, cluster, index, relevant_pos):
	# Count zeroes and one for every position
	poswise_allelecount = dict()
	num_zero = {}
	num_one = {}
	for read in cluster:
		for var in readset[read]:
			pos = index[var.position]
			if pos not in poswise_allelecount:
				poswise_allelecount[pos] = dict()
			if var.allele not in poswise_allelecount[pos]:
				poswise_allelecount[pos][var.allele] = 0
			poswise_allelecount[pos][var.allele] += 1
			
	# Determine majority allele
	cluster_consensus = {}
	for pos in relevant_pos:
		if pos in poswise_allelecount:
			max_allele = 0
			max_count = 0
			sum_count = 0
			for allele in poswise_allelecount[pos]:
				cur_count = poswise_allelecount[pos][allele]
				sum_count += cur_count
				if cur_count > max_count:
					max_allele = allele
					max_count = cur_count
			cluster_consensus[pos] = (max_allele, max_count / sum_count)
		else:
			cluster_consensus[pos] = (0, 1.0)
	
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

def compute_tuple_genotype_dist(consensus,tup, var, geno):
	diff = 0
	tup_geno = dict()
	for i in tup:
		allele = consensus[var][i]
		if allele not in tup_geno:
			tup_geno[allele] = 0
		tup_geno[allele] += 1
	for allele in geno:
		if allele not in tup_geno:
			diff += geno[allele]
		else:
			diff += abs(geno[allele] - tup_geno[allele])
	return diff

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

def compute_linkage_based_block_starts(readset, pos_index, ploidy, single_linkage = False):
	num_vars = len(pos_index)

	# cut threshold
	min_block_size = 1
	if ploidy == 2 or single_linkage:
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
	
	# start by looking at neighbouring
	link_to_next = [0 for i in range(num_vars)]
	starts = []
	ends = []
	for read in readset:
		pos_list = [pos_index[var.position] for var in read]
		starts.append(pos_list[0])
		ends.append(pos_list[-1])
		for i in range(len(pos_list)-1):
			if pos_list[i] + 1 == pos_list[i+1]:
				link_to_next[pos_list[i]] += 1
				
	pos_clust = [0 for i in range(num_vars)]
	for i in range(1, num_vars):
		if link_to_next[i-1] >= cut_threshold:
			pos_clust[i] = pos_clust[i-1]
		else:
			pos_clust[i] = pos_clust[i-1] + 1
	num_clust = pos_clust[-1]+1
			
	# find linkage between clusters
	link_coverage = [dict() for i in range(num_clust)]
	for i, (start, end) in enumerate(zip(starts, ends)):
		covered_pos_clusts = set([pos_clust[pos_index[var.position]] for var in readset[i]])
		for p1 in covered_pos_clusts:
			for p2 in covered_pos_clusts:
				if p2 not in link_coverage[p1]:
					link_coverage[p1][p2] = 0
				link_coverage[p1][p2] += 1
					
	# merge clusters
	merged_clust = [-1 for i in range(num_clust)]
	new_num_clust = 0
	for i in range(num_clust):
		if merged_clust[i] >= 0:
			continue
		q = Queue()
		q.put(i)
		while not q.empty():
			cur = q.get()
			merged_clust[cur] = new_num_clust
			for linked in link_coverage[cur]:
				if merged_clust[linked] < 0 and link_coverage[cur][linked] >= cut_threshold:
					q.put(linked)
		new_num_clust += 1
		
	# determine cut positions
	cuts = [0]
	for i in range(1, num_vars):
		if merged_clust[pos_clust[i]] != merged_clust[pos_clust[i-1]]:
			cuts.append(i)

	return cuts