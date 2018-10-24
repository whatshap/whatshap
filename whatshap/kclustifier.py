from copy import deepcopy

def k_clustify(sim, clusters, k):
	if (len(clusters) == k):
		return clusters
	elif (len(clusters) > k):
		return k_clustify_merge(sim, clusters, k)
	else:
		return k_clustify_split(sim, clusters, k)
	
def k_clustify_merge(sim, clusters, k):
	# Preprocessing
	m = len(clusters)
	clust_dist = [[calc_inter_cluster_weight(sim, clusters[c1], clusters[c2]) for c2 in range(c1+1, m)] for c1 in range(m-1)]
	clust_dist_abs = [[calc_inter_cluster_weight_abs(sim, clusters[c1], clusters[c2]) for c2 in range(c1+1, m)] for c1 in range(m-1)]
	sorted_by_size = sorted(list(range(m)), key = lambda x: -len(clusters[x]))

	# Lists to manage clusters
	active = [sorted_by_size[0]]
	unprocessed = [c for c in range(m) if c not in active]
	cluster_list = [[c] for c in range(m)]
	print("active: "+str(active))
	print("unprocessed: "+str(unprocessed))
	print("cluster_list: "+str(cluster_list))

	# Find initial k clusters
	while len(active) < k:
		c = find_most_overlapping_cluster(clust_dist_abs, unprocessed, active)
		active.append(c)
		active.sort()
		unprocessed.remove(c)

	print("active: "+str(active))
	print("unprocessed: "+str(unprocessed))
	print("cluster_list: "+str(cluster_list))

	# Until <certain condition> add the currently most overlappint cluster and remove on out of k+1
	while(len(unprocessed) > 0):
		# add new cluster to active clusters (so we have k+1 now)
		c = find_most_overlapping_cluster(clust_dist_abs, unprocessed, active)
		if (len(clusters[c]) < 0.1*len(clusters[sorted_by_size[k-1]])):
			cluster_list[c] = []
			unprocessed.remove(c)
			continue
		assert (c in unprocessed)
		active.append(c)
		assert(len(active) == k+1)
		active.sort()
		unprocessed.remove(c)

		# find cluster to merge. i1 < i2. i1 and i2 are NO cluster ids, but indices for "active" !
		i1, i2 = find_cheapest_merge(clust_dist, active)
		c1 = active[i1]
		c2 = active[i2]
		active.remove(c2)
		cluster_list[c1].extend(cluster_list[c2])
		cluster_list[c2] = []
		
		# add distances between clusters, according to the merge
		assert(len(unprocessed) == len(set(unprocessed)))
		assert(len(active) == len(set(active)))
		for x in set(unprocessed).union(set(active)).difference(set([c1,c2])):
			sim[min(x,c1)][max(x,c1)-min(x,c1)-1] += sim[min(x,c2)][max(x,c2)-min(x,c1)-1]
			sim[min(x,c2)][max(x,c2)-min(x,c2)-1] = 0
		
		#print("active: "+str(active))
		#print("unprocessed: "+str(unprocessed))
		#print("cluster_list: "+str(cluster_list))

	# Construct result
	merged_clusters = [merged_cluster for merged_cluster in cluster_list if len(merged_cluster) > 0]
	print("merged_clusters: "+str(merged_clusters))
	return [[ item for i in merged_cluster for item in clusters[i] ] for merged_cluster in merged_clusters]

def k_clustify_split(sim, clusters, k):
	#TODO: Implement!
	return deepcopy(clusters)

def find_most_overlapping_cluster(clust_overlap, unprocessed, active_list):
	# Searches for the cluster, which has the highest overlap to all clusters in the list
	# Different overlaps to clusters in the list are aggegrated by the minimum operator
	best_c = -1
	best_overlap = -float("inf")
	for c in unprocessed:
		overlap = min([clust_overlap[min(x,c)][max(x,c)-min(x,c)-1] for x in active_list])
		if (overlap > best_overlap):
			best_c = c
			best_overlap = overlap
	return best_c

def find_cheapest_merge(clust_dist, active_list):
	c1, c2 = -1, -1
	best_score = -float("inf")
	for a1 in range(len(active_list)-1):
		for a2 in range(a1+1, len(active_list)):
			assert (a1 < len(active_list))
			assert (a2 < len(active_list))
			assert (active_list[a1] < len(clust_dist))
			assert (active_list[a2] - active_list[a1] - 1 < len(clust_dist[ active_list[a1] ]))
			#print(str(a1)+"\t"+str(a2)+"\t"+str(len(active_list)))
			#print(str(active_list[a1])+"\t"+str(active_list[a2])+"\t"+str(len(clust_dist))+"\t"+str(len(clust_dist[active_list[a1]])))
			score = clust_dist[ active_list[a1] ][ active_list[a2] - active_list[a1] - 1]
			if score > best_score:
				best_score = score
				c1, c2 = a1, a2
	return c1, c2

def calc_inter_cluster_weight(sim, cluster1, cluster2):
	sum = 0
	for read1 in cluster1:
		for read2 in cluster2:
			if (read1 < read2):
				sum += sim[read1][read2-read1-1]
			else:
				sum += sim[read2][read1-read2-1]
	return sum

def calc_inter_cluster_weight_abs(sim, cluster1, cluster2):
	sum = 0
	for read1 in cluster1:
		for read2 in cluster2:
			if (read1 < read2):
				sum += abs(sim[read1][read2-read1-1])
			else:
				sum += abs(sim[read2][read1-read2-1])
	return sum

def read_similarities(path):
	read_names = []
	sim = []
	with open(path, "r") as f:
		lines = f.read().splitlines()
		num_reads = int(lines[0])
		for i in range(1, num_reads+1):
			read_names.append(lines[i])
		for i in range(num_reads+1, 2*num_reads):
			tokens = lines[i].split()
			#l = []
			#for token in tokens:
			#	print(token)
			#	l.append(float(token))
			#sim.append(l)
			sim.append([float(s) for s in tokens])
	return read_names, sim

def read_clustering(path, read_names):
	clusters = []
	with open(path, "r") as f:
		lines = f.read().splitlines()
		num_reads = int(lines[0])
		assert(num_reads == len(read_names))
		num_clusters = int(lines[1])
		clusters = [[] for i in range(num_clusters)]
		for i in range(num_reads):
			tokens = lines[i+3].split()
			clusters[int(tokens[1])].append(read_names.index(tokens[0]))
	return clusters

def write_clusters(path, read_names, clusters):
	to_cluster = [-1]*len(read_names)
	for i in range(len(clusters)):
		for read in clusters[i]:
			to_cluster[read] = i
	entries = []
	for i in range(len(read_names)):
		if to_cluster[i] >= 0:
			entries.append((read_names[i], to_cluster[i]))
	#entries = [(read_names[i], to_cluster[i]) for i in range(len(read_names))]
	entries.sort(key = lambda x: x[1])
	
	with open(path, "w") as f:
		f.write(str(len(read_names)))
		f.write("\n")
		f.write(str(len(clusters)))
		f.write("\n")
		f.write("name\tcluster")
		f.write("\n")
		for i in range(len(entries)):
			f.write(entries[i][0]+"\t"+str(entries[i][1]))
			f.write("\n")

def test_main():
	read_names, sim = read_similarities("/home/sven/workspace/wabi-k-clusterediting/k-cluster-editing-paper/snakemake/auxi/sim.2.0.005-p.rg")
	clusters = read_clustering("/home/sven/workspace/wabi-k-clusterediting/k-cluster-editing-paper/snakemake/auxi/sim.2.0.005-p.clust.txt", read_names)
	k_clusters = k_clustify(sim, clusters, 8)
	write_clusters("/home/sven/workspace/wabi-k-clusterediting/k-cluster-editing-paper/snakemake/auxi/sim.2.0.005-p.k-clust.txt", read_names, k_clusters)

if __name__ == "__main__":
    test_main()
