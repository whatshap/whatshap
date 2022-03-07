#include "transitionprobabilitycomputer.h"
#include <cmath>
#include <iostream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(const float& recombcost, Column& column, const vector<int>& current_allele_reference, const vector<int>& next_allele_reference) {
    int s = next_allele_reference.size();
    float r = recombcost;
    double pr = (1.0 - exp(-(r/s))) / s;
    double qr = exp(-(r/s)) + pr;
    vector<unsigned int>* read_ids = column.get_read_ids();
    vector<unsigned int>* next_read_ids = column.get_next_read_ids();
    vector<unsigned int>* active_nonterminating_read_ids = column.get_active_nonterminating_read_ids();
    vector<unsigned int>* active_terminating_read_ids = column.get_active_terminating_read_ids();
    // The +1 is to account for cases where all the samples have different alleles and none of them are ref allele
    vector<vector<unsigned int> > allele_to_reference(next_allele_reference.size()+1);
    int n_alleles = 0;
    for (int i = 0; i < next_allele_reference.size(); i++) {
        if (n_alleles < next_allele_reference[i]) {
            n_alleles = next_allele_reference[i];
        }
        if (next_allele_reference[i] == -1) continue;
        allele_to_reference[next_allele_reference[i]].push_back(i);
    }
    n_alleles++;
    allele_to_reference.resize(n_alleles);
    transition_probability_matrices.resize(n_alleles*n_alleles);
    reordering_map.resize(n_alleles*n_alleles);
    vector<unsigned int> paths_i;
    vector<unsigned int> paths_j;
    vector<unsigned int> previous_references;
    unsigned int matrix_index;
    unsigned int f_index;
    long double value;
    vector<unsigned int>::iterator reference_i;
    vector<unsigned int>::iterator reference_j;
    // Iterating over the alleles present at the position for the first haplotype
    for (int i = 0; i < n_alleles; i++) {
        // Extracting the paths which have the allele for the first haplotype
        paths_i = allele_to_reference[i];
        // Iterating over the alleles present at the position for the second haplotype
        for (int j = 0; j < n_alleles; j++) {
            // Extracting the paths which have the allele for the second haplotype
            paths_j = allele_to_reference[j];
            matrix_index = i*n_alleles + j;
            transition_probability_matrices[matrix_index].remake(pow(s, 2), paths_i.size()*paths_j.size(), 0.0);
            f_index = 0;
            // Iterating over the nodes which contain the allele pair
            for (reference_i = paths_i.begin(); reference_i < paths_i.end(); reference_i++) {
                for (reference_j = paths_j.begin(); reference_j < paths_j.end(); reference_j++) {
                    // Iterating over the indices in the previous column
                    for (unsigned int p_index = 0; p_index < pow(s, 2); p_index++) {
                        previous_references = column.index_to_reference_allele(p_index, 0);
                        if ((*reference_i == previous_references[0]) && (*reference_j == previous_references[1])) value = qr*qr;
                        if ((*reference_i != previous_references[0]) && (*reference_j == previous_references[1])) value = pr*qr;
                        if ((*reference_i == previous_references[0]) && (*reference_j != previous_references[1])) value = pr*qr;
                        if ((*reference_i != previous_references[0]) && (*reference_j != previous_references[1])) value = pr*pr;
                        if ((current_allele_reference.at(previous_references[0]) == -1) || (current_allele_reference.at(previous_references[1]) == -1)) value = 0.0;
                        transition_probability_matrices[matrix_index].set(p_index, f_index, value);
                    }
                    reordering_map[matrix_index].push_back(column.reference_allele_to_index(*reference_i, *reference_j));
                    f_index++;
                }
            }
            
        }
    }
}

Vector2D<long double> * TransitionProbabilityComputer::get_transition_matrix(unsigned int a1, unsigned int a2) {
    return &(this->transition_probability_matrices[a1*(int)(sqrt(transition_probability_matrices.size()))+a2]);
}

vector<unsigned int> * TransitionProbabilityComputer::get_reordering_map(unsigned int a1, unsigned int a2) {
    return &(this->reordering_map[a1*(int)(sqrt(reordering_map.size()))+a2]);
}
