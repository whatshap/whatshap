#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

#include <map>
#include "vector2d.h"
#include "column.h"

class TransitionProbabilityComputer {
private:
    std::vector<Vector2D<long double> > transition_probability_matrices;

    /* Maps the values in the group of states (representing one allele pair) to the r_indices values.
    So reordering_map[allele_index] contains the r_index values of the states which have the alleles given by allele_index.
    The r_index values are ordered which means that the first r_index value is the r_index of the first alpha value calculated in the matrix multiplication.
    */ 
    std::vector<std::vector<unsigned int> > reordering_map;

public:
    TransitionProbabilityComputer(const float& recombcost, Column& column, const std::vector<int>& current_allele_reference, const std::vector<int>& next_allele_reference);

    Vector2D<long double> * get_transition_matrix(unsigned int a1, unsigned int a2);

    std::vector<unsigned int> * get_reordering_map(unsigned int a1, unsigned int a2);
};

#endif // TRANSITIONPROBABILITYCOMPUTER_H
