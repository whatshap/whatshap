import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
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

def cluster_and_draw(output, readset, ploidy, errorrate, min_overlap):
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
			graph.setWeight(id1, id2, similarities[id1][id2 - id1 - 1])

	print("Solving cluster editing ...")
	# run cluster editing
	clusterediting = CoreAlgorithm(graph)	
	readpartitioning = clusterediting.run()
	
	print("Merging clusters ...")
	#kclusters = k_clustify(similarities, readpartitioning, ploidy)

	print("Generating plots ...")
	draw_heatmaps(readset, readpartitioning, output+".heatmaps/")
	print("... finished")
