import operator
from scipy.stats import binom_test
from collections import defaultdict

from .core import Read, ReadSet, Variant
from whatshap.polyphaseplots import parse_haplotype, get_phase
from .threading import get_position_map

def correct_rare_alleles(readset, ploidy, max_dist = 10, alpha = 0.01, ground_truth_haplotypes = None):
	# alpha = level of significant for rare alleles. lower value = detect less, but more confident alleles

	num_vars = len(readset.get_positions())
	
	# for every pair of relevant positions, count the absolute frequency of the different allele-pairs
	#semantics: allelefreq[first position][second position - first position - 1][allele tuple] = absolute frequency
	allelefreq = [[defaultdict(int) for i in range(max_dist)] for j in range(num_vars - 1)]
	
	# for every pair of relevant positions, count the absolute frequency of all alleles
	totalfreq = [[0 for i in range(max_dist)] for j in range(num_vars - 1)]
	
	# mapping between genome and allele matrix positions
	index, rev_index = get_position_map(readset)
	
	# count absolute allele-pair frequencies over all reads
	for read in readset:
		for i in range(len(read)):
			for j in range(i+1, min(i+max_dist+1, len(read))):
				first_pos = index[read[i].position]
				second_pos = index[read[j].position] - first_pos - 1
				if second_pos < max_dist:
					allele_pair = (int(read[i].allele), int(read[j].allele))
					allelefreq[first_pos][second_pos][allele_pair] += 1
					totalfreq[first_pos][second_pos] += 1

	# for every pair of relevant positions: perform statistical to find significant allele-pairs
	rare = [[defaultdict(lambda: False) for i in range(max_dist)] for j in range(num_vars - 1)]
	for pos in range(num_vars-1):
		for offset in range(max_dist):
			for allele_pair in allelefreq[pos][offset]:
				k = allelefreq[pos][offset][allele_pair]
				n = totalfreq[pos][offset]
				# test is hypothesis test of binary distribution. default hypothesis is that every allele_pair is at least as
				# frequent as we would except it to be, if at least one haplotype has this allele_pair (i.e. 1/ploidy)
				p_val = binom_test(k, n=n, p=1/ploidy, alternative="less")
				if p_val < alpha:
					rare[pos][offset][allele_pair] = True

	# Check for available ground truth
	compare = False
	if ground_truth_haplotypes != None and len(ground_truth_haplotypes) == ploidy and len(ground_truth_haplotypes[0]) == num_vars:
		compare = True

	# iterate over reads and correct minimal set of alleles in every read
	corrections = 0
	false_corrections = 0
	alleles = 0
	false_alleles = 0
	for read in readset:
		'''
		pairs = store all rare position pairs
		x = and for every position inside the read, how often it occurs in these pairs
		m = largest entry of x
		'''
		pairs = []
		x = [0 for i in range(len(read))]
		m = 0
		
		# find rare position pairs
		for i in range(len(read)):
			for j in range(i+1, min(i+max_dist, len(read))):
				first_pos = index[read[i].position]
				second_pos = index[read[j].position] - first_pos - 1
				if second_pos < max_dist:
					allele_pair = (int(read[i].allele), int(read[j].allele))
					if rare[first_pos][second_pos][allele_pair]:
						x[i] += 1
						x[j] += 1
						m = max(m, x[i], x[j])
						pairs.append((i,j))
						
		# solve vertex cover heuristically: nodes = positions, edges = rare position pairs
		while m >= 2:
			for i in range(len(read)):
				found = False
				if x[i] == m and not found:
					# determine best allele for replacement
					solved_pairs = [(j,k) for (j,k) in pairs if j == i or k == i]
					allele_count = defaultdict(int)
					for (j,k) in solved_pairs:
						pos = index[read[j].position]
						offset = index[read[k].position] - pos - 1
						for allele_pair in allelefreq[pos][offset]:
							if j == i:
								allele_count[allele_pair[1]] += allelefreq[pos][offset][allele_pair]
							else:
								allele_count[allele_pair[0]] += allelefreq[pos][offset][allele_pair]
					replacement = int(max(allele_count.items(), key=operator.itemgetter(1))[0])
					
					read[i] = Variant(position=read[i].position, allele = 1 - read[i].allele, quality = read[i].quality)
					#read[i] = Variant(position=read[i].position, allele = replacement, quality = read[i].quality)
					
					# Just for debugging and development!!!
					if compare:
						true_hap = parse_haplotype(read.name)
						var_pos = index[read[i].position]
						new_allele = read[i].allele
						correct_allele = int(ground_truth_haplotypes[true_hap][var_pos])
						if correct_allele != new_allele:
							false_corrections += 1

					corrections += 1
					pairs = [(j,k) for (j,k) in pairs if j != i and k != i]
					x = [0 for i in range(len(read))]
					m = 0
					for (j,k) in pairs:
						x[j] += 1
						x[k] += 1
						m = max(m, x[j], x[j])
					found = True
					break
		
		# Just for debugging and development!!!
		if compare:
			for i in range(len(read)):
				true_hap = parse_haplotype(read.name)
				var_pos = index[read[i].position]
				old_allele = read[i].allele
				correct_allele = int(ground_truth_haplotypes[true_hap][var_pos])
				if correct_allele != old_allele:
					false_alleles += 1
				alleles += 1

	# Just for debugging and development!!!
	if compare:
		print("Corrected "+str(corrections)+" alleles.")
		if corrections > 0:
			print("... of which {} were wrong ({:.2f}%)".format(false_corrections, 100*false_corrections/corrections))
		print("Non-Corrected "+str(alleles-corrections)+" alleles.")
		if alleles > 0:
			print("... of which {} were wrong ({:.2f}%)".format(false_alleles-false_corrections, 100*(false_alleles-false_corrections)/(alleles-corrections)))

		print("Global allele error rate: {:.2f}%".format(100*(false_alleles)/(alleles)))
		print("Correction recall:        {:.2f}%".format(100*(corrections-false_corrections)/(corrections-false_corrections+false_alleles)))
		print("Correction precision:     {:.2f}%".format(100*(corrections-false_corrections)/(corrections)))