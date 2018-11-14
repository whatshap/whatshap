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
from .readscoring import score, parse_haplotype
from .kclustifier import k_clustify
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
		
def draw_superheatmap(readset, clustering, var_table, path):
	# Sort a deep copy of clustering
	#clusters = sorted(deepcopy(clustering), key = lambda x: -len(x))
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))

	# Construct real and predicted haplotype per read
	true_hap = [parse_haplotype(read.name) for read in readset]

	# Map variant positions to [0,l)
	index = {}
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1

	# Plot heatmaps
	fig = plt.figure(figsize=(num_vars/100, len(readset)/200), dpi=200)
	legend_handles = {}
	y_offset = 0
	y_margin = 5
	
	# Plot haplotype dissimilarity
	if var_table != None:
		tmp_table = deepcopy(var_table)
		tmp_table.subset_rows_by_position(readset.get_positions())
		phase_rows = [variant.phase for variant in tmp_table.phases[0]]
		phase_vectors = [[row[i] for row in phase_rows] for i in range(len(phase_rows[0]))]
		chunk = 24
		for i, j in it.combinations(range(len(phase_vectors)), 2):
			y_offset -= 104 + y_margin
			colors = ['C'+str(i), 'C'+str(j)]
			if colors[0] not in legend_handles:
				legend_handles[colors[0]] = mpatches.Patch(color=colors[0], label=i)
			if colors[1] not in legend_handles:
				legend_handles[colors[1]] = mpatches.Patch(color=colors[1], label=j)
			dist = [y_offset + 2 + 100*v for v in haplodist(phase_vectors[i], phase_vectors[j], 15)]
			#print(str(i)+":"+str(j)+" = "+str(haplodist(phase_vectors[i], phase_vectors[j], 15)[1300:1450]))
			plt.hlines(y=y_offset, xmin=0, xmax=num_vars, color='black', lw=1)
			plt.hlines(y=y_offset+104, xmin=0, xmax=num_vars, color='black', lw=1)
			for k in range(ceil(num_vars/chunk)):
				start = k*chunk
				end = min(num_vars, (k+1)*chunk+1)
				plt.plot(list(range(start, end)), dist[start:end], lw=1, color=colors[k % 2])
	
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
			plt.hlines(y=read_id+y_offset, xmin=start, xmax=end, color=color_code)
			
		y_offset += len(clusters[c_id]) + y_margin
		plt.hlines(y=y_offset, xmin=0, xmax=num_vars, color=(0.5, 0.5, 0.5, 0.5))
		y_offset += y_margin
		
	plt.legend(handles=legend_handles.values(), loc='lower center', ncol=len(legend_handles))
	axes = plt.gca()
	axes.set_xlim([0, num_vars])
	fig.savefig(path)
	fig.clear()
	
def haplodist(h1, h2, padding):
	if len(h1) != len(h2):
		return -1
	n = len(h1)
	return [relative_hamming_dist(h1[max(0,i-padding):min(n,i+padding+1)], h2[max(0,i-padding):min(n,i+padding+1)]) for i in range(0, n)]

def relative_hamming_dist(seq1, seq2):
	if len(seq1) != len(seq2):
		return -1
	else:
		return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]]) / len(seq1)
	
def draw_cluster_coverage(readset, clustering, path):
	# Sort a deep copy of clustering
	clusters = sorted(deepcopy(clustering), key = lambda x: min([readset[i][0].position for i in x]))

	# Map variant positions to [0,l)
	index = {}
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1

	# Plot heatmaps
	fig = plt.figure(figsize=(num_vars/25, len(readset)/200), dpi=200)
	
	# Coverage plots
	coverage = [[0]*num_vars for i in range(len(clusters))]
	for c_id in range(0, len(clusters)):
		read_id = 0
		for read in clusters[c_id]:
			start = index[readset[read][0].position]
			end = index[readset[read][-1].position]
			for pos in range(start, end+1):
				coverage[c_id][pos] += 1
	
	coverage_sum = [sum([coverage[i][j] for i in range(len(clusters))]) for j in range(num_vars)]
	#coverage_sum = [80 for j in range(num_vars)]
	for c_id in range(0, len(clusters)):
		padding = 40
		coverage[c_id] = [sum(coverage[c_id][max(0, i-padding):min(i+padding+1, num_vars)]) / sum(coverage_sum[max(0, i-padding):min(i+padding+1, num_vars)]) for i in range(num_vars)]
		#coverage[c_id] = [coverage[c_id][i] / coverage_sum[i] for i in range(num_vars)]

	for pos in range(num_vars):
		assign = assign_ploidy(4, [coverage[c_id][pos] for c_id in range(len(clusters))])
		for c_id in range(len(clusters)):
			coverage[c_id][pos] = assign[c_id]+c_id/(10*len(clusters))
	
	for c_id in range(0, len(clusters)):
		if len(clusters[c_id]) >= 10:
			plt.plot(list(range(num_vars)), coverage[c_id], lw=1)
	
	axes = plt.gca()
	axes.set_xlim([0, num_vars])
	fig.savefig(path)
	fig.clear()
	
def assign_ploidy(ploidy, rel_cov):
	cov = sorted([(i, rel_cov[i]) for i in range(len(rel_cov))], key = lambda x: x[1], reverse=True)
	
	# Only look at k biggest values (rest will have ploidy 0 anyways): Only neighbouring integers can be optimal, otherwise error > 1/k
	possibilities = [[floor(cov[i][1]*ploidy), ceil(cov[i][1]*ploidy)] for i in range(ploidy)]
	#print(possibilities)
	min_cost = len(rel_cov)
	min_comb = [1]*ploidy
	#for comb in it.product(*possibilities):
	#	print(comb)
	for comb in it.product(*possibilities):
		if sum(comb) == ploidy:
			cur_cost = sum([abs(comb[i]/ploidy - cov[i][1]) for i in range(ploidy)])
			if cur_cost < min_cost:
				min_cost = cur_cost
				min_comb = comb
	#print(str(cov[:ploidy]) + " -> " + str(min_comb))
	for i in range(ploidy):
		cov[i] = (cov[i][0], min_comb[i])
	for i in range(ploidy, len(cov)):
		cov[i] = (cov[i][0], 0)
		
	cov.sort(key = lambda x: x[0])
	return [cov[i][1] for i in range(len(cov))]				
	
def cluster_and_draw(output, readset, ploidy, errorrate, min_overlap, var_table=None):
	print("Clustering reads for homogeneity plots.")
	print("Computing similarities ...")
	similarities = score(readset, ploidy, errorrate, min_overlap)

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
	#kclusters = k_clustify(similarities, readpartitioning, ploidy)

	print("Generating plots ...")
	#draw_heatmaps(readset, readpartitioning, output+".heatmaps/")
	#draw_superheatmap(readset, readpartitioning, var_table, output+".superheatmap.png")
	draw_cluster_coverage(readset, readpartitioning, output+".superheatmap.png")
	print("... finished")
