#include "transitionprobabilitycomputer.h"
#include <math.h>
#include<iostream>
#include <cassert>

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(unsigned int recombcost, unsigned int trio_count, unsigned int allele_assignments, unsigned int transmission_configurations)
    :transition_prob_matrix(transmission_configurations,transmission_configurations,0.0L),
      transmission_configurations(transmission_configurations)

{

    // precompute bernoulli distribution
    long double recomb_prob = pow(10,-(long double)(recombcost)/10.0L);
    std::vector<long double> bernoulli;
    bernoulli.reserve(2*trio_count);
    for(unsigned int i=0; i <= 2*trio_count; i++){
      bernoulli.emplace_back(pow(recomb_prob,i)*pow(1-recomb_prob,2*trio_count-i));
    }

    for(size_t i = 0; i < transmission_configurations; i++){
        // each row must sum up to 1 and consider also all genotype combinations
        long double normalization_sum = 0.0L;
        for(size_t j = 0; j < transmission_configurations; j++){
            size_t x = i ^ j;
            // count how many bits are set
            x = popcount(x);
            long double prob = bernoulli[x];
            transition_prob_matrix.set(i,j, prob);
            normalization_sum += prob*allele_assignments;
        }
        // normalize row
        for(size_t j = 0; j < transmission_configurations; j++){
            transition_prob_matrix.at(i,j) /= normalization_sum;
        }
    }
}

long double TransitionProbabilityComputer::get(unsigned int t1, unsigned int t2)
{
    assert(t1 < transmission_configurations);
    assert(t2 < transmission_configurations);
    return transition_prob_matrix.at(t1,t2);
}

size_t TransitionProbabilityComputer::popcount(size_t& x) {
    unsigned int count = 0;
    for (;x; x >>= 1) {
        count += x & 1;
    }
    return count;
}
