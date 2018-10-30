import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import savefig
from .readscoring import calc_overlap_and_diffs, parse_haplotype

def draw_plots(readset, path, min_overlap = 5):
	num_reads = len(readset)
	overlap, diffs = calc_overlap_and_diffs(readset)
	haps = [parse_haplotype(readset[i].name) for i in range(num_reads)]
	dissims_same = []
	dissims_diff = []

	for i, j in itertools.combinations(range(num_reads), 2):
		if (overlap[i][j-i-1] >= min_overlap):
			d = diffs[i][j-i-1] / overlap[i][j-i-1]
			if (haps[i] == haps[j]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(dissims_same, dissims_diff, 100, path)

#Counts the fraction of ones in each column of the matrix
def createHistogram(same, diff, steps, path):
	hist = {}
	right_bound = 1.0
	bins = [i*right_bound/steps for i in range(steps+1)]
	plt.hist(same, bins, alpha=0.5, label='same')
	plt.hist(diff, bins, alpha=0.5, label='diff')
	plt.title("Read-pair dissimilarities")
	plt.xlabel("Dissimilarity")
	plt.ylabel("Frequency")
	plt.legend(loc='upper center')
	savefig(path, bbox_inches='tight')
	plt.close()
