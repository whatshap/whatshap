#include "transitionprobabilitycomputer.h"
#include <math.h>
#include<iostream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(size_t column_index, unsigned int recombcost, const Pedigree* pedigree, const std::vector<PedigreePartitions*>& pedigree_partitions)
    :transmission_configurations(std::pow(4, pedigree->triple_count())),
     allele_assignments(1<<pedigree_partitions[0]->count()),
     transitions_transmissions(transmission_configurations,transmission_configurations,0.0L),
     pedigree(pedigree),
     pedigree_partitions(pedigree_partitions),
     transitions_allele_assignments(transmission_configurations,allele_assignments)
{
    size_t trio_count = pedigree->triple_count();

    // precompute bernoulli distribution
    long double recomb_prob = pow(10,-(long double)(recombcost)/10.0L);
    std::vector<long double> bernoulli;
    bernoulli.reserve(2*trio_count);
    for(unsigned int i=0; i <= 2*trio_count; ++i){
      bernoulli.emplace_back(pow(recomb_prob,i)*pow(1-recomb_prob,2*trio_count-i));
    }

    for(size_t i = 0; i < transmission_configurations; ++i){
        // each row must sum up to 1 and consider also all genotype combinations
        long double normalization_sum = 0.0L;
        for(size_t j = 0; j < transmission_configurations; ++j){
            size_t x = i ^ j;
            // count how many bits are set
            x = popcount(x);
            long double prob = bernoulli[x];
            transitions_transmissions.set(i,j, prob);
            normalization_sum += prob;
        }
        // normalize row
        for(size_t j = 0; j < transmission_configurations; ++j){
            transitions_transmissions.at(i,j) /= normalization_sum;
        }
    }

    // compute transition probabilities corresponding to allele assignments
    for(size_t i = 0; i < transmission_configurations; ++i){
        // maps genotype vectors to the number of possible allele assignments
        std::map< std::vector<unsigned int>, size_t > genotypes_to_haplotype_counts;
        // maps a haplotype to the corresponding genotype
        std::vector< std::vector<unsigned int> > haplotypes_to_genotypes(allele_assignments);
        for(unsigned int a = 0; a < allele_assignments; ++a){
            long double prob = 1.0L;
            vector<unsigned int> genotype_vector;
            for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
                unsigned int partition0 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,0);
                unsigned int partition1 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,1);
                unsigned int allele0 = (a >> partition0) & 1;
                unsigned int allele1 = (a >> partition1) & 1;

                int genotype = allele0 + allele1;
                const PhredGenotypeLikelihoods* gls = pedigree->get_genotype_likelihoods(individuals_index, column_index);
                assert(gls != nullptr);
                prob *= gls->get(genotype);
                genotype_vector.push_back(genotype);
            }

            // keep the results
            genotypes_to_haplotype_counts[genotype_vector] += 1;
            transitions_allele_assignments.set(i,a,prob);
            haplotypes_to_genotypes[a] = genotype_vector;
        }

        // divide each probability by the number of times the genotype vector occurs
        long double normalization_sum = 0.0L;
        for(unsigned int a = 0; a < allele_assignments; ++a){
            transitions_allele_assignments.at(i,a) /= genotypes_to_haplotype_counts[haplotypes_to_genotypes[a]];
            normalization_sum += transitions_allele_assignments.at(i,a);
        }

        // normalize the probabilities
        for(unsigned int a = 0; a < allele_assignments; ++a){
            transitions_allele_assignments.at(i,a) /= normalization_sum;
        }
    }
}

long double TransitionProbabilityComputer::get_prob_transmission(unsigned int t1, unsigned int t2)
{
    assert(t1 < transmission_configurations);
    assert(t2 < transmission_configurations);
    return transitions_transmissions.at(t1,t2);
}

long double TransitionProbabilityComputer::get_prob_allele_assignment(unsigned int t, unsigned int a){
    return transitions_allele_assignments.at(t,a);
}

size_t TransitionProbabilityComputer::popcount(size_t& x) {
    unsigned int count = 0;
    for (;x; x >>= 1) {
        count += x & 1;
    }
    return count;
}
