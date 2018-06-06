#include <cassert>

#include "pedigreepartitions.h"

using namespace std;

PedigreePartitions::PedigreePartitions(const Pedigree& pedigree, unsigned int transmission_vector) 
	: pedigree(pedigree),
	transmission_vector(transmission_vector),
	ploidy(pedigree.get_ploidy()),
	haplotype_to_partition_map(pedigree.size(), vector<int>(ploidy, -1))
	{

	// combining partitions that are IDB is currently only allowed for diploid cases.
	// For higher ploidy, treat every individual separately
	// TODO is that correct for polyploid case?
	if (ploidy > 2) assert(pedigree.triple_count() == 0);
	partition_count = ploidy * (pedigree.size() - pedigree.triple_count());

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
			for (unsigned int h = 0; h < ploidy; h++){
				haplotype_to_partition_map[i][h] = p;
				p += 1;
			}
		}
	}

	if (ploidy == 2){
		for (size_t i=0; i<pedigree.size(); ++i) {
			compute_haplotype_to_partition_rec(i, triple_indices);
		}
	}
}

// only for diploid case. Assume, in case of polyploidy, there are no trio relationships given.
void PedigreePartitions::compute_haplotype_to_partition_rec(size_t i,  const vector<int>& triple_indices) {
	if (haplotype_to_partition_map[i][0] != -1) return;
	int triple_index = triple_indices[i];
	assert(triple_index >=0);
	int parent0 = pedigree.get_triples()[triple_index][0];
	int parent1 = pedigree.get_triples()[triple_index][1];
	compute_haplotype_to_partition_rec(parent0, triple_indices);
	compute_haplotype_to_partition_rec(parent1, triple_indices);
	haplotype_to_partition_map[i] = {
	  haplotype_to_partition_map[parent0][!(bool)((transmission_vector >> (2*triple_index)) & 1)],
	  haplotype_to_partition_map[parent1][!(bool)((transmission_vector >> (2*triple_index+1)) & 1)]
	};
}


unsigned int PedigreePartitions::count() const {
	return partition_count;
}


size_t PedigreePartitions::haplotype_to_partition(size_t individual_index, size_t haplotype) const {
	assert(haplotype < pedigree.get_ploidy());
	return haplotype_to_partition_map[individual_index][haplotype];
}

std::ostream& operator<<(std::ostream& out, const PedigreePartitions& pp) {
	for (size_t i=0; i<pp.pedigree.size(); ++i) {
		out << "sample" << i << ":";
		for (size_t h=0; h<pp.ploidy; ++h) {
			out << "  hap" << h << "-->" << pp.haplotype_to_partition(i, h);
		}
		out << endl;
	}
	return out;
}
