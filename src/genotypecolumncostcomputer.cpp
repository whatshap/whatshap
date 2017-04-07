#include <cassert>
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

void GenotypeColumnCostComputer::set_partitioning(unsigned int p) {
    cost_partition.assign(pedigree_partitions.count(), {1.0L,1.0L});
    partitioning = p;
    for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
        auto & entry = **it;
        bool  entry_in_partition1 = (p & ((unsigned int) 1)) == 0;
        unsigned int    ind_id = read_marks[entry.get_read_id()];
        switch(entry.get_allele_type()) {
        case Entry::REF_ALLELE:
            if(entry_in_partition1){
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
            } else {
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
            }
            break;
        case Entry::ALT_ALLELE:
            if(entry_in_partition1){
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            } else {
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
                cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            }
            break;
        case Entry::BLANK:
            break;
        default:
            assert(false);
        }
        p = p >> 1;
    }
}

void GenotypeColumnCostComputer::update_partitioning(int bit_to_flip) {
    const Entry& entry = *column[bit_to_flip];
    // update the partitioning by flipping the given bit
    partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
    // check if the entry is in partition 1
    bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
    unsigned int ind_id = read_marks[entry.get_read_id()];

    // update the costs
    switch(entry.get_allele_type()){
    case Entry::REF_ALLELE:
        if(entry_in_partition1) { 
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] /= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] /= pow(10,-(long double)(entry.get_phred_score())/10.0L);
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
        } else {
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] /= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] /= pow(10,-(long double)(entry.get_phred_score())/10.0L);
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
        }
        break;
    case Entry::ALT_ALLELE:
        if(entry_in_partition1) {
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] /= pow(10,-(long double)(entry.get_phred_score()/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] /= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] *= pow(10,-(long double)(entry.get_phred_score()/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
        } else {
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][0] /= pow(10,-(long double)(entry.get_phred_score()/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,1)][1] /= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][0] *= pow(10,-(long double)(entry.get_phred_score())/10.0L);
            cost_partition[pedigree_partitions.haplotype_to_partition(ind_id,0)][1] *= (1.0L-pow(10,-(long double)(entry.get_phred_score())/10.0L));
        }
        break;
    case Entry::BLANK:
        break;
    default:
        assert(false);
    }
}

long double GenotypeColumnCostComputer::get_cost(unsigned int allele_assignment) {
    long double cost = 1.0L;
    // for the given allele assignment multiply the costs of the partitions
    for(size_t p = 0; p < pedigree_partitions.count(); p++){
        // get the allele corresponding to the partition
        unsigned int allele = (allele_assignment >> p) & 1;
        cost *= cost_partition[p][allele];
    }
    return cost;
}
