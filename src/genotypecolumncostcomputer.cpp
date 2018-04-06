#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "vector2d.h"
#include "genotypecolumncostcomputer.h"

using namespace std;

GenotypeColumnCostComputer::GenotypeColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<unsigned int>& read_marks, const Pedigree* pedigree, const PedigreePartitions& pedigree_partitions)
    :column(column),
     column_index(column_index),
     read_marks(read_marks),
     partitioning(0),
     pedigree(pedigree),
     cost_partition(pedigree_partitions.count(),{1.0L,1.0L}),
     pedigree_partitions(pedigree_partitions)

{}

namespace {
  array<long double, 256> precompute_phred_probabilities() {
       array<long double, 256> result;
    result[0] = 0.9999;
    for(auto i = 1; i < 256; ++i){
      result[i] = pow(10, -i/10.0L);
    }
    return result;
  }

  array<long double, 256> phred_probability_small = precompute_phred_probabilities();
  unordered_map<unsigned int, long double> phred_probability;

  long double get_phred_probability(unsigned int phred_score) {
    if(phred_score < 256) {
      return phred_probability_small[phred_score];
    }
    auto it = phred_probability.find(phred_score);
    if(it != phred_probability.end()) {
      return it->second;
    } else {
      auto res = pow(10, -(int)phred_score/10.0L);
      phred_probability.emplace(phred_score, res);
      return res;
    }
  }
}

void GenotypeColumnCostComputer::set_partitioning(unsigned int p) {
    cost_partition.assign(pedigree_partitions.count(), {1.0L,1.0L});
    partitioning = p;
    for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
        auto & entry = **it;
        if(entry.get_allele_type() == Entry::BLANK) {
            continue;
        }
        bool  entry_in_partition1 = (p & ((unsigned int) 1)) == 0;
        unsigned int    ind_id = read_marks[entry.get_read_id()];
        bool is_ref_allele = entry.get_allele_type() == Entry::REF_ALLELE;

        auto proba = get_phred_probability(entry.get_phred_score());
        cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,entry_in_partition1)][!is_ref_allele] *= (1.0L-proba);
        cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,entry_in_partition1)][is_ref_allele] *= proba;
        p = p >> 1;
    }
}

void GenotypeColumnCostComputer::update_partitioning(int bit_to_flip) {
    const Entry& entry = *column[bit_to_flip];
    if(entry.get_allele_type() == Entry::BLANK) {
      return;
    }

    // update the partitioning by flipping the given bit
    partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
    // check if the entry is in partition 1
    bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
    unsigned int ind_id = read_marks[entry.get_read_id()];

    // update the costs
    bool is_ref_allele = entry.get_allele_type() == Entry::REF_ALLELE;

    auto proba = get_phred_probability(entry.get_phred_score());
    cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,entry_in_partition1)][!is_ref_allele] *= (1.0L-proba);
    cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,entry_in_partition1)][is_ref_allele] *= proba;
    cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,!entry_in_partition1)][!is_ref_allele] /= (1.0L-proba);
    cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,!entry_in_partition1)][is_ref_allele] /= proba;
}

long double GenotypeColumnCostComputer::get_cost(unsigned int allele_assignment) {
    long double cost = 1.0L;
    // for the given allele assignment multiply the costs of the partitions
    for(size_t p = 0; p < pedigree_partitions.count(); ++p){
        // get the allele corresponding to the partition
        unsigned int allele = (allele_assignment >> p) & 1;
        cost *= cost_partition[p][allele];
    }
    return cost;
}
