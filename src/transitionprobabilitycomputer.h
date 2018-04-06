#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

#include <map>
#include "vector2d.h"
#include "pedigree.h"
#include "pedigreepartitions.h"

class TransitionProbabilityComputer {
private:
    unsigned int transmission_configurations;
    unsigned int allele_assignments;
    Vector2D<long double> transitions_transmissions;
    size_t popcount(size_t& x);

    const Pedigree* pedigree;
    const std::vector<PedigreePartitions*>& pedigree_partitions;
    // transitions to allele assignments
    Vector2D<long double> transitions_allele_assignments;

public:
    TransitionProbabilityComputer(size_t column_index, unsigned int recombcost, const Pedigree* pedigree, const std::vector<PedigreePartitions*>& pedigree_partitions);
    // get the transision probability for change of transmission vector t1 to t2
    long double get_prob_transmission(unsigned int t1, unsigned int t2);
    long double get_prob_allele_assignment(unsigned int t, unsigned int a);
};

#endif // TRANSITIONPROBABILITYCOMPUTER_H
