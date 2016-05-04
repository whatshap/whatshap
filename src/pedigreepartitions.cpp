#include <cassert>

#include "pedigreepartitions.h"

using namespace std;

PedigreePartitions::PedigreePartitions(const Pedigree& pedigree, unsigned int transmission_vector) : pedigree(pedigree), transmission_vector(transmission_vector), haplotype_to_partition_map(pedigree.size(), std::make_pair(-1,-1)) {
	partition_count = 2 * (pedigree.size() - pedigree.triple_count());
	
	// for each individual the index of the triple in which this individual is a child (-1 for "none")
	vector<int> triple_indices(pedigree.size(), -1);
	for (size_t i=0; i<pedigree.triple_count(); ++i) {
		triple_indices[pedigree.get_triples()[i][2]] = i;
	}
	// assign partition numbers to "roots", i.e. to those individuals that don't 
	// have parents in our pedigree
	int p = 0;
	for (size_t i=0; i<pedigree.size(); ++i) {
		if (triple_indices[i] == -1) {
			haplotype_to_partition_map[i] = std::make_pair(p, p+1);
			p += 2;
		}
	}
	for (size_t i=0; i<pedigree.size(); ++i) {
		compute_haplotype_to_partition_rec(i, triple_indices);
	}
}

void PedigreePartitions::compute_haplotype_to_partition_rec(size_t i,  const vector<int>& triple_indices) {
	if (haplotype_to_partition_map[i].first != -1) return;
	int triple_index = triple_indices[i];
	assert(triple_index >=0);
	int parent0 = pedigree.get_triples()[triple_index][0];
	int parent1 = pedigree.get_triples()[triple_index][1];
	compute_haplotype_to_partition_rec(parent0, triple_indices);
	compute_haplotype_to_partition_rec(parent1, triple_indices);
	haplotype_to_partition_map[i] = std::make_pair(
		((transmission_vector >> (2*triple_index)) & 1) ? haplotype_to_partition_map[parent0].first : haplotype_to_partition_map[parent0].second,
		((transmission_vector >> (2*triple_index+1)) & 1) ? haplotype_to_partition_map[parent1].first : haplotype_to_partition_map[parent1].second
	);
}

unsigned int PedigreePartitions::count() const {
	return partition_count;
}

size_t PedigreePartitions::haplotype_to_partition(size_t individual_index, size_t haplotype) const {
	if (haplotype == 0) {
		return haplotype_to_partition_map[individual_index].first;
	} else {
		assert(haplotype == 1);
		return haplotype_to_partition_map[individual_index].second;
	}
}
