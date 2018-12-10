import numpy as np
import itertools as it
from math import ceil, floor
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import shutil
import os
from .core import Read, ReadSet, CoreAlgorithm, LightCompleteGraph
from .readscoring import score, locality_sensitive_score, parse_haplotype
from .kclustifier import clusters_to_haps, clusters_to_blocks, avg_readlength, calc_consensus_blocks
from .vectorerror import vector_error, vector_error_blockwise
from copy import deepcopy

def draw_heatmaps(readset, clustering, heatmap_folder):

	# Sort a deep copy of clustering
	clusters = sorted(deepcopy(clustering), key = lambda x: -len(x))

	# Construct real and predicted haplotype per read
	true_hap = [parse_haplotype(read.name) for read in readset]

	# Map variant positions to [0,l)
	index = {}
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1
		
	# Plot heatmaps
	if os.path.exists(heatmap_folder):
		shutil.rmtree(heatmap_folder)
	os.makedirs(heatmap_folder)

	for c_id in range(0, len(clusters)):
		fig = plt.figure(figsize=(12, 6), dpi=100)
		legend_handles = {}
		SCALING_FACTOR = 1

		read_id = 0
		for read in clusters[c_id]:   
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			read_id += 1
			
			color_code = 'C'+str(true_hap[read]) if true_hap[read] != '-1' else 'black'
			if color_code not in legend_handles:
				legend_handles[color_code] = mpatches.Patch(color=color_code, label=true_hap[read])
			plt.hlines(y=read_id*SCALING_FACTOR, xmin=start, xmax=end, color=color_code)
		
		plt.legend(handles=legend_handles.values())
		axes = plt.gca()
		axes.set_xlim([0, num_vars])
		heatmap_path = heatmap_folder+"/cluster_"+str(c_id)+".png"
		fig.savefig(heatmap_path)			
		fig.clear()
		
def draw_superheatmap(readset, clustering, var_table, path, genome_space = False):
	# Sort a deep copy of clustering
	#clusters = sorted(deepcopy(clustering), key = lambda x: -len(x))
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))

	# Construct real and predicted haplotype per read
	true_hap = [parse_haplotype(read.name) for read in readset]

	# Map variant positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	min_pos = float("inf")
	max_pos = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1
		min_pos = min(min_pos, position)
		max_pos = max(max_pos, position)
		
	min_pos = min(readset.get_positions()) if genome_space else 0
	max_pos = max(readset.get_positions()) if genome_space else num_vars

	# Plot heatmaps
	fig = plt.figure(figsize=(num_vars/100, len(readset)/200), dpi=200)
	fig = plt.figure(figsize=(num_vars/100, len(readset)/200), dpi=200)
	legend_handles = {}
	y_offset = 0
	y_margin = 5
	
	# Plot haplotype dissimilarity
	if var_table != None:
		plot_haplotype_dissimilarity(legend_handles, y_offset, y_margin, index, rev_index, readset, var_table, genome_space)
	
	y_offset = 0
	
	# Plot heatmaps
	for c_id in range(0, len(clusters)):
		if len(clusters[c_id]) < 5:
			continue
		read_id = 0
		for read in clusters[c_id]:   
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			read_id += 1
			
			color_code = 'C'+str(true_hap[read]) if true_hap[read] != '-1' else 'black'
			if color_code not in legend_handles:
				legend_handles[color_code] = mpatches.Patch(color=color_code, label=true_hap[read])
			if genome_space:
				plt.hlines(y = read_id + y_offset, xmin = rev_index[start], xmax = rev_index[end], color = color_code)
			else:
				plt.hlines(y = read_id + y_offset, xmin = start, xmax = end, color = color_code)
			
		y_offset += len(clusters[c_id]) + y_margin
		plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color=(0.5, 0.5, 0.5, 0.5))
		y_offset += y_margin
		
	plt.legend(handles=legend_handles.values(), loc='lower center', ncol=len(legend_handles))
	axes = plt.gca()
	axes.set_xlim([min_pos, max_pos])
	fig.savefig(path)
	fig.clear()
	
def draw_cluster_coverage(readset, clustering, path):
	ploidy = 4
	
	# Sort a deep copy of clustering
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))

	# Map variant positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1

	# Plot heatmaps
	fig = plt.figure(figsize=(num_vars/25, len(readset)/200), dpi=200)
	
	# Get coverage from external function
	coverage, copynumbers, cluster_blocks, cut_positions = clusters_to_blocks(readset, clustering, ploidy, coverage_padding = 7, copynumber_max_artifact_len = 1.0, copynumber_cut_contraction_dist = 0.5)
	
	for c_id in range(0, len(clusters)):
		if len(clusters[c_id]) >= 10:
			plt.plot(list(range(num_vars)), copynumbers[c_id], lw=1)
	
	axes = plt.gca()
	axes.set_xlim([0, num_vars])
	fig.savefig(path)
	fig.clear()
	
def draw_cluster_blocks(readset, clustering, cluster_blocks, cut_positions, var_table, path, genome_space):
	# Sort a deep copy of clustering
	#clusters = sorted(deepcopy(clustering), key = lambda x: -len(x))
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))

	# Construct real and predicted haplotype per read
	true_hap = [parse_haplotype(read.name) for read in readset]

	# Map variant positions to [0,l)
	index = {}
	rev_index = []
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		rev_index.append(position)
		num_vars += 1

	# Plot heatmaps
	fig = plt.figure(figsize=(num_vars/100, len(readset)/200), dpi=200)
	legend_handles = {}
	y_offset = 0
	y_margin = 5
	
	# Plot haplotype dissimilarity
	if var_table != None:
		plot_haplotype_dissimilarity(legend_handles, y_offset, y_margin, index, rev_index, readset, var_table, genome_space)
	
	y_offset = 0
	
	# Plot cluster blocks
	for i, block in enumerate(cluster_blocks):
		y_offset = 0
		# x-Positions for current block in plot
		end = cut_positions[i]
		start = cut_positions[i-1] if i > 0 else 0
		for j, hap in enumerate(block):
			# Split each haplotype by cluster: [(cluster1, start1, end1), (cluster2, start2, end2), ... ]
			sections = []
			section_start = 0
			assert end - start == len(hap)
			for k in range(len(hap)-1):
				if hap[k] != hap[k+1]:
					if hap[k] != None:
						sections.append((hap[k], start + section_start, start + k+1))
					section_start = k+1
			if hap[-1] != None:
				sections.append((hap[-1], start + section_start, end))
			
			# Collect all reads for each section, but cut off reads at section borders. Discard reads, which do not overlap the section at all.
			# Reads are saved as: [(truth1, start1, end1), (truth2, start2, end2), ... ]
			reads = []
			for section in sections:
				#print("block="+str(i)+" hap="+str(j)+" section="+str(section))
				for read_id in clusters[section[0]]:
					read_start = max(section[1], index[readset[read_id][0].position])
					read_end = min(section[2]-1, index[readset[read_id][-1].position])
					#print("    read=("+str(read_start)+","+str(read_end)+")")
					if (read_end >= read_start):
						reads.append((true_hap[read_id], read_start, read_end))
			
			# Plot reads:
			#print("len(reads)="+str(len(reads)))
			for k, read in enumerate(reads):
				color_code = 'C'+str(read[0]) if read[0] != '-1' else 'black'
				if color_code not in legend_handles:
					legend_handles[color_code] = mpatches.Patch(color=color_code, label=true_hap[read])
				plt.hlines(y=k+y_offset, xmin=read[1], xmax=read[2], color=color_code)
				
			# Haplotype seperator
			y_offset += len(reads) + y_margin
			plt.hlines(y = y_offset, xmin = start, xmax = end, color=(0.5, 0.5, 0.5, 0.5))
			y_offset += y_margin
			
		# Plot block seperators
		plt.vlines(x = start, ymin = 0, ymax = y_offset - y_margin, color=(0.5, 0.5, 0.5, 0.5))
		plt.vlines(x = end, ymin = 0, ymax = y_offset - y_margin, color=(0.5, 0.5, 0.5, 0.5))
		
	plt.legend(handles=legend_handles.values(), loc='lower center', ncol=len(legend_handles))
	axes = plt.gca()
	axes.set_xlim([0, num_vars])
	fig.savefig(path)
	fig.clear()
	
def plot_haplotype_dissimilarity(legend_handles, y_offset, y_margin, index, rev_index, readset, var_table, genome_space = False):
	
	num_vars = len(readset.get_positions())
	min_pos = min(readset.get_positions()) if genome_space else 0
	max_pos = max(readset.get_positions()) if genome_space else num_vars
	readlen = avg_readlength(readset)
	
	# Plot heatmaps
	y_offset = 0
	y_margin = 5
	
	# Plot haplotype dissimilarity
	#tmp_table = deepcopy(var_table)
	#tmp_table.subset_rows_by_position(readset.get_positions())
	#phase_rows = [variant.phase for variant in tmp_table.phases[0]]
	#phase_vectors = [[row[i] for row in phase_rows] for i in range(len(phase_rows[0]))]
	phase_vectors = get_phase(readset, var_table)
	chunk = 24
	padding = readlen/2 # dissimilarity of position i is averaged over interval of 2*padding base pairs (not positions)

	if genome_space:
		# get variant density
		dens_pos = list(range(min_pos+padding, max_pos-padding, 2*padding))
		dens_list = [len(list(filter(lambda x: pos-padding <= x <= pos+padding, rev_index))) for pos in dens_pos]
		max_dens = max(dens_list)
		print("max_dens="+str(max_dens))

		y_offset -= 104 + y_margin
		plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color='black', lw=1)
		plt.hlines(y=y_offset+104, xmin=min_pos, xmax=max_pos, color='black', lw=1)
		plt.plot(dens_pos, [100*x/max_dens + y_offset for x in dens_list], lw=1, color='blue')

	# determines for each position, over which interval of positions the average must be taken
	intervals = []
	for i in range(num_vars):
		left = right = i
		pos = rev_index[i]
		while left - 1 >= 0 and rev_index[left-1] >= pos - padding:
			left -= 1
		while right + 1 < num_vars and rev_index[right+1] <= pos + padding:
			right += 1
		intervals.append([left, right])

	# One plot for each pair of haplotypes
	for i, j in it.combinations(range(len(phase_vectors)), 2):
		y_offset -= 104 + y_margin
		colors = ['C'+str(i), 'C'+str(j)]
		if colors[0] not in legend_handles:
			legend_handles[colors[0]] = mpatches.Patch(color=colors[0], label=i)
		if colors[1] not in legend_handles:
			legend_handles[colors[1]] = mpatches.Patch(color=colors[1], label=j)
		dist = [y_offset + 2 + 100*v for v in haplodist(phase_vectors[i], phase_vectors[j], intervals)]
		plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color='black', lw=1)
		plt.hlines(y=y_offset+104, xmin=min_pos, xmax=max_pos, color='black', lw=1)
		for k in range(ceil(num_vars/chunk)):
			start = k*chunk
			end = min(num_vars, (k+1)*chunk+1)
			if genome_space:
				plt.plot(rev_index[start:end], dist[start:end], lw=1, color=colors[k % 2])
			else:
				plt.plot(list(range(start, end)), dist[start:end], lw=1, color=colors[k % 2])
				
def get_phase(readset, var_table):
	tmp_table = deepcopy(var_table)
	tmp_table.subset_rows_by_position(readset.get_positions())
	phase_rows = [variant.phase for variant in tmp_table.phases[0]]
	return [[row[i] for row in phase_rows] for i in range(len(phase_rows[0]))]
	
def haplodist(h1, h2, intervals):
	if len(h1) != len(h2):
		return -1
	n = len(h1)
	return [relative_hamming_dist(h1[intervals[i][0]:min(n, intervals[i][1]+1)], h2[intervals[i][0]:min(n, intervals[i][1]+1)]) for i in range(0, n)]

def relative_hamming_dist(seq1, seq2):
	if len(seq1) != len(seq2):
		return -1
	else:
		return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]]) / len(seq1)
	
def construct_ground_truth(readset):
	true_hap = [parse_haplotype(read.name) for read in readset]
	real_clusters = [[read_id for read_id in range(len(readset)) if true_hap[read_id] == hap_id] for hap_id in range(4)]
	num_vars = len(readset.get_positions())
	trivial_cluster_blocks = [[[hap_id]*num_vars for hap_id in range(4)]]
	return calc_consensus_blocks(readset, real_clusters, trivial_cluster_blocks, [])[0]

def construct_ground_truth_blockwise(readset, cut_positions):
	true_hap = [parse_haplotype(read.name) for read in readset]
	real_clusters = [[read_id for read_id in range(len(readset)) if true_hap[read_id] == hap_id] for hap_id in range(4)]
	num_vars = len(readset.get_positions())
	trivial_cluster_blocks = []
	start = 0
	for pos in cut_positions:
		trivial_cluster_blocks.append([[hap_id]*(pos-start) for hap_id in range(4)])
		start = pos
	return calc_consensus_blocks(readset, real_clusters, trivial_cluster_blocks, cut_positions)
	
def cluster_and_draw(output, readset, ploidy, errorrate, min_overlap, var_table=None):
	print("Clustering reads for homogeneity plots.")
	print("Computing similarities ...")
	similarities = score(readset, ploidy, errorrate, min_overlap)
	#similarities = locality_sensitive_score(readset, ploidy, min_overlap)

	print("Constructing graph ...")
	# create read graph object
	graph = LightCompleteGraph(len(readset),True)

	# insert edges into read graph
	n_reads = len(readset)
	for id1 in range(n_reads):
		for id2 in range(id1+1, n_reads):
			if similarities.get(id1, id2) != 0:
				graph.setWeight(id1, id2, similarities.get(id1, id2))

	print("Solving cluster editing ...")
	# run cluster editing
	clusterediting = CoreAlgorithm(graph)	
	readpartitioning = clusterediting.run()
	
	print("Merging clusters ...")
	#coverage, copynumbers, cluster_blocks, cut_positions = clusters_to_blocks(readset, readpartitioning, ploidy, coverage_padding = 7, copynumber_max_artifact_len = 0.5, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = True)
	#consensus_blocks = clusters_to_haps(readset, readpartitioning, ploidy, coverage_padding = 7, copynumber_max_artifact_len = 0.5, copynumber_cut_contraction_dist = 0.5, single_hap_cuts = True)

	#truth = get_phase(readset, var_table)
	#truth_cons = construct_ground_truth(readset)
	
	#truth_blocks = construct_ground_truth_blockwise(readset, cut_positions)
	#for hap in truth:
	#	print("".join(list(map(lambda x: str(x) if x >= 0 else ".", hap))))
	#for hap in truth_cons:
	#	print("".join(list(map(lambda x: str(x) if x >= 0 else ".", hap))))
	#vec_error = vector_error_blockwise(consensus_blocks, truth, 1, 100000)
	#print("Vector error = "+str(vec_error))
	
	#vec_error = vector_error_blockwise(consensus_blocks, truth_cons, 1, 100000)
	#print("Vector error = "+str(vec_error))
	
	print("Generating plots ...")
	#draw_heatmaps(readset, readpartitioning, output+".heatmaps/")
	#draw_superheatmap(readset, readpartitioning, var_table, output+".superheatmap.png")
	draw_superheatmap(readset, readpartitioning, var_table, output+".superheatmapg.png", genome_space = False)
	#draw_cluster_coverage(readset, readpartitioning, output+".coverage.png")
	#draw_cluster_blocks(readset, readpartitioning, cluster_blocks, cut_positions, var_table, output+".haploblocks.png", genome_space = False)
	print("... finished")