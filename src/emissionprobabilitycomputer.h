#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include "vector2d.h"
#include "column.h"
#include "readset.h"

class EmissionProbabilityComputer {

private:
    std::vector<long double> emission_values;
    unsigned int b_index;
    unsigned int r_index;

public:
    
    EmissionProbabilityComputer(Column* column, std::vector<const Entry*>& entries, const unsigned int n_allele);

    long double get_emission_value(unsigned int b_index, unsigned int allele_1, unsigned int allele_2, const unsigned int n_allele);
};

#endif