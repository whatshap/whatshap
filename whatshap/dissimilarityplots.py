import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import savefig
from .readscoring import calc_overlap_and_diffs, parse_haplotype, score

def draw_plots_dissimilarity(readset, path, min_overlap = 5, steps = 100):
	num_reads = len(readset)
	overlap, diffs = calc_overlap_and_diffs(readset)
	haps = [parse_haplotype(readset[i].name) for i in range(num_reads)]
	dissims_same = []
	dissims_diff = []

	for i, j in itertools.combinations(range(num_reads), 2):
		if (overlap[i][j-i-1] >= min_overlap):
			d = diffs[i][j-i-1] / overlap[i][j-i-1]
			d = similarities[i][j-i-1]
			if (haps[i] == haps[j]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(path, dissims_same, dissims_diff, steps, [0.0, 1.0], "Dissimilarity")
	
def draw_plots_scoring(readset, path, ploidy, error_rate, min_overlap = 5, steps=120, dim=[-60, 60]):
	num_reads = len(readset)
	overlap, diffs = calc_overlap_and_diffs(readset)
	similarities = score(readset, ploidy, error_rate, min_overlap)
	haps = [parse_haplotype(readset[i].name) for i in range(num_reads)]
	dissims_same = []
	dissims_diff = []

	for i, j in itertools.combinations(range(num_reads), 2):
		if (overlap[i][j-i-1] >= min_overlap):
			d = similarities[i][j-i-1]
			if (haps[i] == haps[j]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(path, dissims_same, dissims_diff, steps, dim, "Similarity score")

#Counts the fraction of ones in each column of the matrix
def createHistogram(path, same, diff, steps, dim, x_label):
	hist = {}
	left_bound = dim[0]
	right_bound = dim[1]
	bins = [left_bound + i*(right_bound-left_bound)/steps for i in range(steps+1)]
	plt.hist(same, bins, alpha=0.5, label='same')
	plt.hist(diff, bins, alpha=0.5, label='diff')
	plt.title("Read-pair comparison")
	plt.xlabel(x_label)
	plt.ylabel("Frequency")
	plt.legend(loc='upper center')
	savefig(path, bbox_inches='tight')
	plt.close()
