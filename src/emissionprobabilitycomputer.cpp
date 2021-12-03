#include "emissionprobabilitycomputer.h"
#include <cmath>
#include <iostream>
#include <cassert>

#include "columnindexingiterator.h"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(Column* column, vector<const Entry*>& entries, const unsigned int n_allele) {
    emission_values.resize(pow(n_allele,2)*pow(2,column->get_read_ids()->size()));
    vector<unsigned int>* read_ids = column->get_read_ids();
    long double value;
    for (int allele_1 = 0; allele_1 < n_allele; allele_1++) {
        for (int allele_2 = 0; allele_2 < n_allele; allele_2++) {
            int a_index = allele_1*n_allele + allele_2;
            value = 1.0;
            unique_ptr<ColumnIndexingIterator> iterator = unique_ptr<ColumnIndexingIterator>(new ColumnIndexingIterator(column));
            while (iterator->has_next()) {
                int bit_changed = -1;
                iterator->advance(&bit_changed);
                vector<int> binaryVector = iterator->get_binary_vector();
                int b_index = iterator->get_b_index();
                int index = b_index*pow(n_allele, 2) + (n_allele*allele_1 + allele_2);
                if (bit_changed == -1) {
                    for (int read_index = 0; read_index < column->get_read_ids()->size(); read_index++) {
                        if (entries.at(read_index)->get_allele_type() == -1) {
                            continue;
                        }
                        value = value*(entries.at(read_index)->get_emission_score()[allele_1]);
                    }
                    emission_values[index] = value;
                }
                int new_bit = binaryVector[bit_changed];
                if (entries.at(bit_changed)->get_allele_type() == -1) {
                        emission_values[index] = value;
                        continue;
                    }
                if (new_bit == 1) {
                    value = value * (entries.at(bit_changed)->get_emission_score()[allele_2]/entries.at(bit_changed)->get_emission_score()[allele_1]);
                }
                else {
                    value = value * (entries.at(bit_changed)->get_emission_score()[allele_1]/entries.at(bit_changed)->get_emission_score()[allele_2]);
                }
                emission_values[index] = value;
            }
        }
    }
}

long double EmissionProbabilityComputer::get_emission_value(unsigned int b_index, unsigned int allele_1, unsigned int allele_2, const unsigned int n_allele) {
    return emission_values.at(b_index*pow(n_allele, 2) + (n_allele*allele_1 + allele_2));
}