import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import savefig
from .readscoring import calc_overlap_and_diffs, parse_haplotype, score, locality_sensitive_score

def draw_plots_dissimilarity(readset, path, min_overlap = 5, steps = 100):
	num_reads = len(readset)
	overlap, diffs = calc_overlap_and_diffs(readset)
	haps = [parse_haplotype(readset[i].name) for i in range(num_reads)]
	dissims_same = []
	dissims_diff = []

	for i, j in itertools.combinations(range(num_reads), 2):
		if (overlap.get(i, j) >= min_overlap):
			d = diffs.get(i, j) / overlap.get(i, j)
			if (haps[i] == haps[j]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(path, dissims_same, dissims_diff, steps, [0.0, 1.0], "Dissimilarity", "Read-pair comparison")
	
def draw_plots_scoring(readset, path, ploidy, error_rate, min_overlap = 5, steps=120, dim=[-60, 60]):
	num_reads = len(readset)
	overlap, diffs = calc_overlap_and_diffs(readset)
	similarities = score(readset, ploidy, error_rate, min_overlap)
	#similarities = locality_sensitive_score(readset, ploidy, min_overlap)
	haps = [parse_haplotype(readset[i].name) for i in range(num_reads)]
	dissims_same = []
	dissims_diff = []

	for i, j in itertools.combinations(range(num_reads), 2):
		if (overlap.get(i, j) >= min_overlap):
			d = similarities.get(i, j)
			if (haps[i] == haps[j]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(path, dissims_same, dissims_diff, steps, dim, "Similarity score", "Read-pair comparison")
	
def draw_column_dissimilarity(readset, path, steps = 100):
	num_reads = len(readset)
	alleles = [[0]*4 for i in readset.get_positions()]
	index = {}
	num_vars = 0
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1
		
	for read in readset:
		for variant in read:
			pos = index[variant.position]
			allele = variant.allele
			if allele > len(alleles[pos]):
				for i in range(len(allele[pos]), allele):
					allele[pos].append(0)
			alleles[index[variant.position]][variant.allele] += 1
			
	sim1 = [max(alleles[i]) / sum(alleles[i]) for i in range(len(alleles))]
	sim2 = [min([alleles[i][j] for j in range(len(alleles[i])) if alleles[i][j] > 0]) / sum(alleles[i]) for i in range(len(alleles))]
	createHistogram(path, sim1, sim2, steps, [0.0, 1.0], "Frequency of most frequent allele", "Column-wise comparison", name1='most freqeunt', name2='least frequent')

#Counts the fraction of ones in each column of the matrix
def createHistogram(path, same, diff, steps, dim, x_label, title, name1='same', name2='diff'):
	hist = {}
	left_bound = dim[0]
	right_bound = dim[1]
	bins = [left_bound + i*(right_bound-left_bound)/steps for i in range(steps+1)]
	plt.hist(same, bins, alpha=0.5, label=name1)
	if len(diff) > 0:
		plt.hist(diff, bins, alpha=0.5, label=name2)
	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel("Frequency")
	plt.legend(loc='upper center')
	savefig(path, bbox_inches='tight')
	plt.close()
