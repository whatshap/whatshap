#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

#include "vector2d.h"

class TransitionProbabilityComputer {
private:
    Vector2D<long double> transition_prob_matrix;
    unsigned int transmission_configurations;
    size_t popcount(size_t& x);
public:
    TransitionProbabilityComputer(unsigned int recombcost, unsigned int trio_count, unsigned int allele_assignments);
    // get the transision probability for change of transmission vector t1 to t2
    long double get(unsigned int t1, unsigned int t2);
};

#endif // TRANSITIONPROBABILITYCOMPUTER_H
