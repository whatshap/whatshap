from copy import deepcopy
from math import floor, ceil
import itertools as it
import logging
from collections import defaultdict
from .core import clustering_DP
import numpy as np

logger = logging.getLogger(__name__)

def clusters_to_haps(readset, clustering, ploidy, coverage_padding = 12, copynumber_max_artifact_len = 1.0, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = False):
	logger.info("   Processing "+str(len(clustering))+" read clusters ...")
	cluster_blocks = []
	# Compute cluster blocks
	if len(clustering) >= ploidy:
		logger.info("   Processing "+str(len(clustering))+" read clusters ...")
		coverage, copynumbers, cluster_blocks, cut_positions = clusters_to_blocks(readset, clustering, ploidy, coverage_padding, copynumber_max_artifact_len, copynumber_cut_contraction_dist, single_hap_cuts)
	else:
		logger.info("   Processing "+str(len(clustering))+" read clusters ... could not create "+str(ploidy)+" haplotypes from this!")
		return []
		
	# Compute consensus blocks
	logger.info("   Received "+str(len(cluster_blocks))+" cluster blocks. Generating consensus sequences ...")

	consensus_blocks = calc_consensus_blocks(readset, clustering, cluster_blocks, cut_positions)
	
	gaps = sum([sum([sum([1 for i in hap if hap[i] == -1]) for hap in block]) for block in consensus_blocks])
	logger.info("   Phased "+str(ploidy)+" haplotypes over "+str(len(readset.get_positions()))+" variant positions, using "+str(len(consensus_blocks))+" blocks with "+str(gaps)+" undefined sites.")
	
	return consensus_blocks

def calc_consensus_blocks(readset, clustering, cluster_blocks, cut_positions):
	cluster_consensus = get_cluster_consensus(readset, clustering)
	consensus_blocks = deepcopy(cluster_blocks)
	for i, block in enumerate(consensus_blocks):
		offset = cut_positions[i-1] if i > 0 else 0
		for j, hap in enumerate(block):
			for pos, clust in enumerate(hap):
				consensus_blocks[i][j][pos] = cluster_consensus[clust][offset+pos] if clust != None else -1
				
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


def subset_clusters(readset, clustering,ploidy):

	# Sort a deep copy of clustering
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))
	readlen = avg_readlength(readset)
	
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1
	print("clustering: ", clustering)
	clusterlengths = [len(c) for c in clustering]
	print("cluster lengths: ", sorted(clusterlengths))
	print("length clustering before: ", len(clustering))
	
	#Subset the set of clusters to clusters with the highest amount of reads
	subset_clustering = [c for c in clustering if len(c) > 20]
	subset_cids = [c_id for c_id in range(len(clustering)) if len(clustering[c_id]) > 20]
	cluster_dict = {}
	print("subset_clustering: ", subset_clustering)
	print("length of subset_clustering: ", len(subset_clustering))
	for i in range(len(clustering)):
		for j in range(len(subset_clustering)):
			if (clustering[i]==subset_clustering[j]):
				cluster_dict[j] = i
	print("cluster dict: ", cluster_dict)
	

	#for every cluster in clustering, compute its sequence of starting and ending positions
	#create dictionary mapping the clusterID to a list of pairs (starting position, ending position)
	print("number of variants: ", num_vars)
	positions = defaultdict(list)
	for c_id in range(0, len(subset_clustering)):
		read_id = 0
		for read in subset_clustering[c_id]:
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			positions[c_id].append((start,end))
	
	#for every cluster and every variant position, compute the relative coverage
	num_vars = len(index)
	coverage = [[0]*num_vars for i in range(len(subset_clustering))]
	for c_id in range(0, len(subset_clustering)):
		read_id = 0
		for read in subset_clustering[c_id]:
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			for pos in range(start, end+1):
				coverage[c_id][pos] += 1

	
	coverage_sum = [sum([coverage[i][j] for i in range(len(subset_clustering))]) for j in range(num_vars)]

	for c_id in range(0, len(subset_clustering)):
	#	coverage[c_id] = [sum(coverage[c_id][i:i+1]) / sum(coverage_sum[0:num_vars]) for i in range(num_vars)]
		coverage[c_id] = [((coverage[c_id][i])/coverage_sum[i]) if coverage_sum[i]!= 0 else 0 for i in range(num_vars)]	

	cov_2 = [coverage[c_id][6] for c_id in range(len(subset_clustering))]
	print("cov 2: ", cov_2)
	solu = []	
	for i in range(num_vars):
		test = False
		for j in range(len(subset_clustering)):
			if (coverage[j][i]==1):
				test = True
		if test:
			solu.append(i)
	print("coverage equals 1: ", solu)
	print("subset clustering: ", subset_clustering)
	

	#map the clusters to a list of IDs
	c_ids = [c_id for c_id in range(0,len(subset_clustering))]
	print("c_ids: ",c_ids)
	print("subset c_ids: ", subset_cids)
	cluster_tuples = [(i,j,k,l) for i in c_ids for j in c_ids for k in c_ids for l in c_ids]
	print("cluster_tuples: ", cluster_tuples[:10], "length: ", len(cluster_tuples))
	#map the cluster tuples to indices from [0, len(tuples)] 	
	c_tuples = {}
	for i in range(0, len(cluster_tuples)):
		c_tuples[i] = cluster_tuples[i]
	
	#initialize matrix: make list for every possible 4-tuple of clusters in clustering to keep the matrix twodimensional
	scoring = [[0]*num_vars for i in range(len(cluster_tuples))]	

	#rows: tuple of cluster IDs, columns: var positions
	#initialize first column: only coverage costs
	for c_tuple in range(len(cluster_tuples)):
		scoring[c_tuple][0] = (cov_costs(cluster_tuples[c_tuple], 0, coverage), -1)
	first_col = [scoring[i][0][0] for i in range(len(cluster_tuples))]
	#map each variant to a list of cluster tuple indices that do not cover the variant
	uncovered = defaultdict(set)
	scoring = clustering_DP(num_vars,cluster_tuples,coverage, positions)
	for i in range(num_vars):
		test = []		
		for j in range(len(cluster_tuples)):
			if (scoring[i][j][0] < 1000000):
				test.append(scoring[i][j])
		if (len(test) == 0):
			print("no clusters: ", i)
			break
#	for var in range(1,num_vars):
#		print("var: ", var)
#		for c_tuple in range(len(cluster_tuples)):
#			if (cov_costs(cluster_tuples[c_tuple], var, coverage) == 1000):
#				scoring[c_tuple][var] = (1000,c_tuple)
#				uncovered[var].add(c_tuple)
#			#minimum of (the column before plus switch costs) plus coverage costs
#			else:
#				#pred = [(scoring[i][var-1][0]+switch_costs(cluster_tuples[c_tuple], cluster_tuples[i],positions,var-1)) for i in range(len(cluster_tuples))]
#				pred = []				
#				for i in range(len(cluster_tuples)):
#					if (i in uncovered[var]):
#						pred.append(1000)
#					else:
#						pred.append(scoring[i][var-1][0]+switch_costs(cluster_tuples[c_tuple], cluster_tuples[i],positions,var-1))
#									
#				minimum, minimum_index = min((val, idx) for (idx,val) in enumerate(pred))
#				scoring[c_tuple][var] = ((cov_costs(cluster_tuples[c_tuple], var, coverage) + minimum), minimum_index)

	#start backtracing from the minimum in the last column
	path = []
	start_col = find_backtracing_start(scoring, num_vars, cluster_tuples)
	print("start col: ", start_col)
	#last_col = [scoring[num_vars-5][i][0] for i in range(len(cluster_tuples))]
	last_col = [scoring[start_col][i][0] for i in range(len(cluster_tuples))]
	first_col = [scoring[0][i][0] for i in range(len(cluster_tuples))]

	second_col = [scoring[1][i][0] for i in range(len(cluster_tuples))]

	last_min_idx = int(min((val, idx) for (idx, val) in enumerate(last_col))[1])
	path.append(cluster_tuples[last_min_idx])
	#append stored predecessor
	for i in range(num_vars-5,0,-1):
		pred = scoring[i][int(last_min_idx)][1]
		path.append(cluster_tuples[int(pred)])
		last_min_idx = pred
	path.reverse()

	counter = 0
	for i in path:
		if (i[0]==i[1]==i[2]==i[3]):
			counter +=1
	print("homzg: ", counter)
	#determine cut positions: When at least two haplotypes change the cluster from position i to i+1, a new cut is made
	#if only one haplotype changes the cluster, both clusters are assumed to belong together (?)
	cut_positions = []
	for i in range(len(path)-1):
		dissim = 0
		for j in range(0,4):
			if path[i][j] != path[i+1][j]:
				dissim +=1
		if (dissim > 1):
			cut_positions.append(i)
	print("cut positions: ", cut_positions)
	#divide path of clusters into 4 paths of haplotypes
	hap1 = [cluster_dict[tup[0]] for tup in path]
	print("hap1: ", hap1, "length: ", len(hap1))
	hap2 = [cluster_dict[tup[1]] for tup in path]
#	print("hap2: ", hap2, "length: ", len(hap2))
	hap3 = [cluster_dict[tup[2]] for tup in path]
#	print("hap3: ", hap3, "length: ", len(hap3))
	hap4 = [cluster_dict[tup[3]] for tup in path]
#	print("hap4: ", hap4, "length: ", len(hap4))
	print("num vars: ", num_vars)
#	assert(len(hap1)==len(hap2)==len(hap3)==len(hap4)==num_vars)

	haps=[hap1,hap2,hap3,hap4]
	# Return cluster blocks: List of blocks, each blocks ranges from one cut position to the next and contains <ploidy> sequences
	# of cluster ids, which indicat e, which cluster goes to which haplotype at which position.
	cluster_blocks = []
	old_pos = 0
	for cut_pos in cut_positions:
		cluster_blocks.append([hap[old_pos:cut_pos] for hap in haps])
		old_pos = cut_pos
	print("cluster_blocks: ", cluster_blocks[0])
	return(coverage,cut_positions, cluster_blocks)

def find_backtracing_start(scoring, num_vars, cluster_tuples):
	[scoring[num_vars-5][i][0] for i in range(len(cluster_tuples))]
	minimum = 1000000
	last_col = num_vars-1
	res = False
	for i in range(len(cluster_tuples)):
		if scoring[last_col][i][0] < minimum:
			res = True
	if res:
		return(last_col)
	else:
		return(find_backtracing_start(scoring,last_col-1,cluster_tuples))
			

def cov_costs(c_tuple, var, coverage):
	costs = 0
	exp_cn = 0
	#compute copy numbers for every cluster in c_tuple
	for i in range(0,4):
		cov = coverage[c_tuple[i]][var]
		#if cluster does not cover the position var:
		if (cov == 0):
			#return (float('inf'))
			return (1000)
		#else compare the expected copy number to the real one
		else:
			if (cov > 0 and cov < 0.25):
				exp_cn = 0
			if (cov >= 0.25 and cov < 0.5):
				exp_cn = 1
			if (cov >= 0.5 and cov < 0.75):
				exp_cn = 2
			if (cov >= 0.75 and cov < 1):
				exp_cn = 3
		cn = c_tuple.count(c_tuple[i])
		if (exp_cn != cn):
			costs+= 1
	return(costs)

def switch_costs(c_tuple1, c_tuple2, positions, var):
	costs = 0
	#switch costs depend on the position: if var is the end of c_tuple1 or var+1 is the beginning of c_tuple2, switching is free
	for i in range(0,4):
		starts_list = [j[0] for j in positions[c_tuple2[i]]]
		ends_list = [j[1] for j in positions[c_tuple1[i]]]
		if (var not in ends_list and var+1 not in starts_list):
			if (c_tuple1[i] != c_tuple2[i]):
				costs += 1
	return(costs)
				
			
def clusters_to_blocks(readset, clustering, ploidy, coverage_padding = 12, copynumber_max_artifact_len = 1.0, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = False):
	# Sort a deep copy of clustering
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))
	readlen = avg_readlength(readset)
	
	# Map genome positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1

	# Get relative coverage for each cluster at each variant position
	coverage = calc_coverage(readset, clusters, coverage_padding, index)
	logger.info("      ... computed relative coverage for all clusters")

	# Assign haploid copy numbers to each cluster at each variant position
	copynumbers = calc_haploid_copynumbers(coverage, num_vars, ploidy)
	logger.info("      ... assigned local copy numbers for all clusters")
	
	# Optimize copynumbers
	postprocess_copynumbers(copynumbers, rev_index, num_vars, ploidy, readlen, copynumber_max_artifact_len, copynumber_cut_contraction_dist)
	
	# Compute cluster blocks
	cut_positions, cluster_blocks = calc_cluster_blocks(readset, copynumbers, num_vars, ploidy, single_hap_cuts)
	logger.info("   Cut positions:")
	print(cut_positions)
	
	return coverage, copynumbers, cluster_blocks, cut_positions

def avg_readlength(readset):
	# Estiamtes the average read length in base pairs
	if len(readset) > 0:
		return sum([read[-1].position - read[0].position for read in readset]) / len(readset)
	else:
		return 0
	
def calc_coverage(readset, clustering, padding, index):
	# Determines for every variant position the relative coverage of each cluster. 
	# "Relative" means, that it is the fraction of the total coverage over all clusters at a certain position
	num_vars = len(index)
	coverage = [[0]*num_vars for i in range(len(clustering))]
	for c_id in range(0, len(clustering)):
		read_id = 0
		for read in clustering[c_id]:
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			for pos in range(start, end+1):
				coverage[c_id][pos] += 1
	
	coverage_sum = [sum([coverage[i][j] for i in range(len(clustering))]) for j in range(num_vars)]
	for c_id in range(0, len(clustering)):
		coverage[c_id] = [sum(coverage[c_id][max(0, i-padding):min(i+padding+1, num_vars)]) / sum(coverage_sum[max(0, i-padding):min(i+padding+1, num_vars)]) for i in range(num_vars)]
	
	# cov[i][j] = relative coverage of cluster i at variant position j
	return coverage

def calc_haploid_copynumbers(coverage, num_vars, ploidy):
	result = deepcopy(coverage)
	
	for pos in range(num_vars):
		# Sort relative coverages in descending order, keep original index at first tuple position
		cn = sorted([(i, coverage[i][pos]) for i in range(len(coverage))], key = lambda x: x[1], reverse=True)

		# Only look at k biggest values (rest will have ploidy 0 anyways): Only neighbouring integers can be optimal, otherwise error > 1/k
		possibilities = [[floor(cn[i][1]*ploidy), ceil(cn[i][1]*ploidy)] for i in range(ploidy)]
		min_cost = len(coverage)
		min_comb = [1]*ploidy
		for comb in it.product(*possibilities):
			if sum(comb) <= ploidy:
				cur_cost = sum([abs(comb[i]/ploidy - cn[i][1]) for i in range(ploidy)])
				if cur_cost < min_cost:
					min_cost = cur_cost
					min_comb = comb
		for i in range(ploidy):
			cn[i] = (cn[i][0], min_comb[i])
		for i in range(ploidy, len(cn)):
			cn[i] = (cn[i][0], 0)

		# Sort computed copy numbers by index
		cn.sort(key = lambda x: x[0])
		new_cov = [cn[i][1] for i in range(len(cn))]
		
		# Write copy numbers into result
		for c_id in range(len(cn)):
			result[c_id][pos] = cn[c_id][1]
	
	return result

def postprocess_copynumbers(copynumbers, rev_index, num_vars, ploidy, readlen, artifact_len, contraction_dist):
	# Construct intervals of distinct assignments, start with first interval
	start = 0
	interval = []
	intervals = []
	for i in range(len(copynumbers)):
		for j in range(copynumbers[i][0]):
			interval.append(i)
	interval.sort() # Sorted multiset of clusters present at the start
	
	for pos in range(1, num_vars):
		# Create multiset of present clusters for current position
		current = []
		for i in range(len(copynumbers)):
			for j in range(copynumbers[i][pos]):
				current.append(i)
		current.sort()
		if current != interval:
			# if assignment changes, append old interval to list and open new one
			intervals.append((start, pos, deepcopy(interval)))
			start = pos
			interval = current
	# Append last opened interval to list
	intervals.append((start, num_vars, interval))
	logger.info("      ... divided variant location space into "+str(len(intervals))+" intervals")
	
	# Phase 1: Remove intermediate deviation
	if artifact_len > 0.0:
		i = 0
		while i < len(intervals):
			to_merge = -1
			for j in range(i+1, len(intervals)):
				if intervals[i][2] == intervals[j][2]:
					len1 = intervals[i][1] - intervals[i][0]
					len2 = intervals[j][1] - intervals[j][0]
					gap = intervals[j][0] - intervals[i][1]
					gaplen = rev_index[intervals[j][0]] - rev_index[intervals[i][1]]
					if len1 > 2*gap and len2 > 2*gap and gaplen < readlen * artifact_len:
						to_merge = j
						break
			if to_merge > i:
				#print("merging interval "+str(i)+": "+str(intervals[i])+" with interval "+str(to_merge)+": "+str(intervals[to_merge]))
				intervals = intervals[:i] + [(intervals[i][0], intervals[to_merge][1], intervals[i][2])] + intervals[to_merge+1:]
				for k in range(intervals[i][0]+1, intervals[i][1]):
					for l in range(len(copynumbers)):
						copynumbers[l][k] = copynumbers[l][intervals[i][0]]
			else:
				i += 1
		logger.info("      ... reduced intervals to "+str(len(intervals))+" by removing intermediate artifacts")
	
	# Questionable: Close k-1-gaps
	if False:
		i = 0
		while i < len(intervals) - 1:
			set1 = set(intervals[i][2])
			set2 = set(intervals[i+1][2])
			#print(str(i)+" int1="+str(intervals[i][2])+" int2="+str(intervals[i+1][2]))
			if set1 >= set2 and len(intervals[i][2]) == len(intervals[i+1][2]) + 1:
				intervals = intervals[:i] + [(intervals[i][0], intervals[i+1][1], intervals[i][2])] + intervals[i+2:]
				for k in range(intervals[i][0]+1, intervals[i][1]):
					for l in range(len(copynumbers)):
						copynumbers[l][k] = copynumbers[l][intervals[i][0]]
			else:
				i += 1
		logger.info("      ... reduced intervals to "+str(len(intervals))+" by extending intervals to uncovered sites")
			
	# Phase 2: Remove intervals with less than 1*read_length
	if contraction_dist > 0.0:
		i = 0
		while i < len(intervals):
			to_merge = -1
			for j in range(i+2, len(intervals)):
				if rev_index[intervals[j][0]] - rev_index[intervals[i][1]] <= readlen * contraction_dist:
					to_merge = j
				else:
					break
			if to_merge > i:
				middle = (intervals[to_merge][0] + intervals[i][1]) // 2
				#print("merging interval "+str(i)+": "+str(intervals[i])+" with interval "+str(to_merge)+": "+str(intervals[to_merge]))
				intervals = intervals[:i] + [(intervals[i][0], middle, intervals[i][2]), (middle, intervals[to_merge][1], intervals[to_merge][2])] + intervals[to_merge+1:]
				for k in range(intervals[i][0]+1, intervals[i][1]):
					for l in range(len(copynumbers)):
						copynumbers[l][k] = copynumbers[l][intervals[i][0]]
				for k in range(intervals[i+1][0], intervals[i+1][1]-1):
					for l in range(len(copynumbers)):
						copynumbers[l][k] = copynumbers[l][intervals[i+1][1]-1]
			else:
				i += 1
		logger.info("      ... reduced intervals to "+str(len(intervals))+" by shifting proximate interval bounds")
	
def calc_cluster_blocks(readset, copynumbers, num_vars, ploidy, single_hap_cuts = False):
	# Get all changes for each position
	num_clust = len(copynumbers)
	
	# For each position: Get clusters which increase/decrease their copynumber right here. Position 0 contains all initial clusters as increasing.
	increasing = []
	decreasing = []
	increasing.append([c for c in range(num_clust) if copynumbers[c][0] > 0])
	decreasing.append([])
	for pos in range(1, num_vars):
		increasing.append([c for c in range(num_clust) if copynumbers[c][pos-1] < copynumbers[c][pos]])
		decreasing.append([c for c in range(num_clust) if copynumbers[c][pos-1] > copynumbers[c][pos]])
		
	# Layout the clusters into haplotypes and indicate the actual set of cut positions
	cut_positions = []
	haps = [[None]*num_vars for i in range(ploidy)]
	
	# Usually, if a cluster increased in copynumber, i.e. from 1 to 2, then there is a potential block border. The same case is a decrease in copynumber.
	# What really forces a new block is, if the copynumber developes like this: x -> y -> z with 0 < x < y > z > 0, i.e. an increase, followed by a decrease. 
	# In this case it is unclear, how the haplotypes from copynumber x match the ones from copynumber z. So we need to cut, either for the shift from x to y 
	# or from y to z or in between. The opposite case (0 < x > y < z > 0, y > 0) has the same problem. Consecutive increases (x < y < z) or decreases (x > y > z)
	# force no problems. We therefore track for each cluster, if at the current position it is allowed to increase or decrease its copynumber.
	increase_disallowed = set()
	decrease_disallowed = set()
	
	# Start with clusters, that are present at position 0
	h = 0
	for c in increasing[0]:
		for i in range(copynumbers[c][0]):
			haps[h][0] = c
			h += 1
	
	# Iterate over all positions
	for pos in range(1, num_vars):
		open_haps = []
		assigned = [0]*num_clust
		for h in range(ploidy):
			c = haps[h][pos-1]
			# If cluster does not decrease or if it decreases, but the new copynumber is not reached yet, just continue with current cluster on h
			if c != None and (c not in decreasing[pos] or assigned[c] < copynumbers[c][pos]):
				haps[h][pos] = c
				assigned[c] += 1
			else:
				# Else, report haplotype slot as open
				open_haps.append(h)
				
		# Assign unmatched cluster copynumbers to remaining, open slots
		i = 0
		for c in range(num_clust):
			for x in range(copynumbers[c][pos] - assigned[c]):
				haps[open_haps[i]][pos] = c
				i += 1
				
		# Determine, whether this position is a new block
		blocking_clusters = 0
		msg = ""
		for h in range(ploidy):
			# Check for the following: If more than one haplotype needs a new block, we have to add this position as a cut position. If it is only one,
			# we can just assume that the old and the new cluster belong together, implied by the other clusters not needing a change.
			# Exception: If a cluster increases or decreases its copynumber without permission, we always cut
			c = haps[h][pos]
			if haps[h][pos-1] == c or c == None or haps[h][pos-1] == None:
				# No cluster change -> nothing to do
				continue
			# Increase number of block forces by one, in all cases.
			blocking_clusters += 1
			if copynumbers[c][pos-1] == 0:
				# new cluster
				msg = msg + "Cluster "+str(c)+" is new and replacing an old cluster. "
				pass
			elif c in increasing[pos] and c in increase_disallowed:
				# copynumber increases, but c already had a decrease before
				msg = msg + "Cluster "+str(c)+" is increasing, but has decresed since the last cut. "
				blocking_clusters += ploidy
			elif c in decreasing[pos] and c in decrease_disallowed:
				# copynumber increases, but c already had a decrease before
				msg = msg + "Cluster "+str(c)+" is decreasing, but has incresed since the last cut. "
				blocking_clusters += ploidy
			else:
				msg = msg + "Cluster "+str(c)+" is in/decreasing ("+str(copynumbers[c][pos-1])+","+str(copynumbers[c][pos])+") "
		
		if (blocking_clusters <= 1 and not single_hap_cuts) or blocking_clusters == 0:
			#if blocking_clusters == 1:
			#	print("Potential cut at pos="+str(pos)+": " + msg + "But single discontinued cluster can be resolved.")
			# No block, but disallow increase/decrease for clusters with changed copynumber
			for c in range(num_clust):
				if c in increasing[pos] and copynumbers[c][pos-1] > 0:
					decrease_disallowed.add(c)
				elif c in decreasing[pos] and copynumbers[c][pos] > 0:
					increase_disallowed.add(c)
		else:
			# Add cut positions and remove prohibitions regarding copynumber in/decreases
			#print(msg + "Making a cut.")
			cut_positions.append(pos)
			decrease_disallowed.clear()
			increase_disallowed.clear()
	
	# Return cluster blocks: List of blocks, each blocks ranges from one cut position to the next and contains <ploidy> sequences
	# of cluster ids, which indicat e, which cluster goes to which haplotype at which position.
	cluster_blocks = []
	old_pos = 0
	for cut_pos in cut_positions:
		cluster_blocks.append([hap[old_pos:cut_pos] for hap in haps])
		old_pos = cut_pos
		
	return cut_positions, cluster_blocks
