import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import savefig

def draw_plots(readset, path):
	num_reads = len(readset)
	num_vars = 0
	index = {}
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1

	m, starts, ends, haps = calc_matrix_and_bounds(readset, index, num_vars)
	dissims_same = []
	dissims_diff = []
	for row1, row2 in itertools.combinations(range(len(m)), 2):
		mi = max(starts[row1], starts[row2])
		ma = min(ends[row1], ends[row2])
		d = dissim(m[row1], m[row2], mi, ma)
		if (d >= 0.0):
			if (haps[row1] == haps[row2]):
				dissims_same.append(d)
			else:
				dissims_diff.append(d)
	createHistogram(dissims_same, dissims_diff, 100, path)

def calc_matrix_and_bounds(readset, index, num_vars):
	num_reads = len(readset)
	m = []
	starts = [num_vars+1]*num_reads
	ends = [-1]*num_reads
	haps = []
	for i in range(num_reads):
		line = [-1]*num_vars
		for variant in readset[i]:
			if (variant.position in index):
				starts[i] = min(starts[i], index[variant.position])
				ends[i] = max(ends[i], index[variant.position])
				line[index[variant.position]] = variant.allele
		m.append(line)
		haps.append(parse_readname(readset[i].name))
	return m, starts, ends, haps

def parse_readname(name):
	tokens = name.split("_")
	if (tokens[-2] == "HG00514" and tokens[-1] == "HAP1"):
		return 0
	elif (tokens[-2] == "HG00514" and tokens[-1] == "HAP2"):
		return 1
	elif (tokens[-2] == "NA19240" and tokens[-1] == "HAP1"):
		return 2
	elif (tokens[-2] == "NA19240" and tokens[-1] == "HAP2"):
		return 3
	print("Error: Could not parse haplotype from read "+name)
	return -1
	
# Compute the relative dissimilarity between two rows
def dissim(row1, row2, start, end):
	pos = 0
	diff = 0
	for i in range(start, end+1, 1):
		if (row1[i] > 0 and row2[i] > 0):
			pos += 1
			if (row1[i] != row2[i]):
				diff += 1
	if (pos < 5):
		return -1
	else:
		return diff/pos

#Counts the fraction of ones in each column of the matrix
def createHistogram(same, diff, steps, path):
	hist = {}
	right_bound = 1.0
	bin_size = int(steps + 1)
	bins = [i*right_bound/steps for i in range(steps+1)]
	plt.hist(same, bins, alpha=0.5, label='same')
	plt.hist(diff, bins, alpha=0.5, label='diff')
	plt.title("Read-pair dissimilarities")
	plt.xlabel("Dissimilarity")
	plt.ylabel("Frequency")
	plt.legend(loc='upper center')
	savefig(path, bbox_inches='tight')
