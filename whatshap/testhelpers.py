"""
Utility functions only used by unit tests
"""
import textwrap
from collections import defaultdict
from whatshap.core import Read, ReadSet, Variant
import math

def string_to_readset(s, w = None, sample_ids = None, source_id=0, scale_quality = None, n_alleles = 2):
	s = textwrap.dedent(s).strip()
	if w is not None:
		w = textwrap.dedent(w).strip().split('\n')
	rs = ReadSet()
	for index, line in enumerate(s.split('\n')):
		if len(line) == 0:
			continue
		if sample_ids is None:
			read = Read('Read {}'.format(index+1), 50, source_id)
		else:
			read = Read('Read {}'.format(index+1), 50, source_id, sample_ids[index])
		for pos, c in enumerate(line):
			if c == ' ':
				continue
			q = 1
			if w is not None:
				q = int(w[index][pos])
			if scale_quality is not None:
				q *= scale_quality
			quality = [q] * n_alleles
			quality[int(c)] = 0
			read.add_variant(position=(pos+1) * 10, allele=int(c), quality=quality)
		assert len(read) > 1, 'Reads covering less than two variants are not allowed'
		rs.add(read)
	print(rs)
	return rs


def string_to_readset_pedigree(s, w = None, scaling_quality = None):
	s = textwrap.dedent(s).strip()
	read_sources = []
	s2 = ''
	for line in s.split('\n'):
		if len(line) == 0:
			continue
		individual = ord(line[0]) - ord('A')
		assert 0 <= individual < 26
		read_sources.append(individual)
		s2 += line[1:] + '\n'
	rs = string_to_readset(s=s2, w=w, sample_ids=read_sources, scale_quality=scaling_quality)
	print('read_sources:', read_sources)
	return rs


def matrix_to_readset(lines) :

	rs = ReadSet()
	index_tracker = 0
	for line in lines :

		s = line.split()
		assert len(s) % 2 == 1, 'Not in matrix format.'

		index = int(s[0])
		index_tracker += 1
		assert index == index_tracker, 'Not in matrix format.'

		read = Read('Read {}'.format(index), 50)
		for i in range(int(len(s) / 2)) :

			offset = int(s[2*i+1])
			for pos, c in enumerate(s[2*i+2]) :
				read.add_variant(position=(offset+pos) * 10, allele=int(c), quality = [1])

		rs.add(read)

	print(rs)
	return rs


def flip_cost(variant, target_value):
	"""Returns cost of flipping the given read variant to target_value."""
#	if variant.allele == target_value:
#		return 0
#	else:
#		return min([x for x in variant.quality if x != 0])

	return variant.quality[target_value]

def is_ambiguous(assignments, ploidy):
	sets = [set() for i in range(ploidy)]
	for assignment in assignments:
		for s, allele in zip(sets, assignment):
			s.add(allele)
	return [len(s) > 1 for s in sets]

def assignment_to_list(assignment, ploidy, n_alleles):
	return tuple( int((assignment/pow(n_alleles,i)) % n_alleles) for i in range(ploidy) )

def column_cost(variants, possible_assignments, ploidy):
	""" Compute cost for one position and return the minimum cost assignment. """
	costs = []
	# go through all allele assignments
	for assignment in possible_assignments:
		# consider cost for each partition
		cost = 0
		for partition in range(ploidy):
			# get allele assigned to current partition:
			allele = assignment[partition];
			cost += sum(flip_cost(v, allele) for v in variants[partition])
		costs.append(cost)
	l = [(cost,i) for i,cost in enumerate(costs)]
	l.sort()
	min_cost = l[0][0]
	best_assignment = list(possible_assignments[l[0][1]])
	counts = defaultdict(int)
	for cost, index in l:
		counts[cost] += 1
	ties = counts[min_cost]
	ambiguous = is_ambiguous([possible_assignments[i] for cost,i in l[:ties]], ploidy)
	for i in range(ploidy):
		if ambiguous[i]:
			best_assignment[i] = -2
	return min_cost, best_assignment

# TODO extend to multiallelic case
def allowed_assignments_for_genotype(genotype, ploidy, n_alleles):
	assignment_count = 1 << ploidy
	result = []
	for i in range(0,assignment_count):
		assignment_list = assignment_to_list(i, ploidy, n_alleles)
		if sum(assignment_list) == genotype:
			result.append(assignment_list)
	return result
	

def brute_force_phase(read_set, ploidy, n_alleles = 2, allowed_genotypes = None):
	"""Solves MEC by enumerating all possible bipartitions."""
	def print(*args): pass

	assert len(read_set) < 10, "Too many reads for brute force"
	positions = read_set.get_positions()
	assignment_count = pow(n_alleles, ploidy)

	# bit i in "partition" encodes to which set read i belongs
	best_partition = None
	best_cost = None
	best_haplotypes = None
	solution_count = 0
	for partition in range(ploidy**len(read_set)):
		print('Looking at partition {{:0>{}b}}'.format(len(read_set)).format(partition))
		# compute cost induced by that partition
		cost = 0
		haplotypes = []
		for index,p in enumerate(positions):
			# find variants covering this position
			variants = [ [] for x in range(ploidy) ]
			for n, read in enumerate(read_set):
				i = (partition // (ploidy**n)) % ploidy
				for variant in read:
					if variant.position == p:
						variants[i].append(variant)
			if allowed_genotypes is None:
				possible_assignments = [ assignment_to_list(i, ploidy, n_alleles) for i in range(0,assignment_count)]
				c, assignment = column_cost(variants, possible_assignments, ploidy)
			else:
				possible_assignments = allowed_assignments_for_genotype(allowed_genotypes[index], ploidy, n_alleles)
				c, assignment = column_cost(variants, possible_assignments, ploidy)
			print('    position: {}, variants: {} --> cost = {}'.format(p, str(variants), c))
			cost += c
			haplotypes.append(assignment)
		print('  --> cost for this partitioning:', cost)
		if (best_cost is None) or (cost < best_cost):
			best_partition = partition
			best_cost = cost
			best_haplotypes = haplotypes
			solution_count = 1
		elif cost == best_cost:
			solution_count += 1

	# Each partition has its inverse with the same cost
	number_of_equal_solutions = math.factorial(ploidy)
	assert solution_count % number_of_equal_solutions == 0
	haplotypes = []
	print(solution_count, number_of_equal_solutions)
	for i in range(ploidy):
		h = ''.join([str(hap[i]) for hap in best_haplotypes])
		haplotypes.append(h)
	return best_cost, [(best_partition//(ploidy**x)) % ploidy for x in range(len(read_set))], solution_count//number_of_equal_solutions, haplotypes
