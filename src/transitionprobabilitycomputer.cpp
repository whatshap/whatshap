#include "transitionprobabilitycomputer.h"
#include <cmath>
#include <iostream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(const unsigned int& recombcost, Column& column, const vector<unsigned int>& allele_reference) {
    pr = (1 - exp(-recombcost/allele_reference.size()))/allele_reference.size();
    qr = exp(-recombcost/allele_reference.size()) + pr;
    vector<unsigned int>* read_ids = column.get_read_ids();
    vector<unsigned int>* next_read_ids = column.get_next_read_ids();
    vector<unsigned int>* active_nonterminating_read_ids = column.get_active_nonterminating_read_ids();
    vector<unsigned int>* active_terminating_read_ids = column.get_active_terminating_read_ids();
    vector<vector<unsigned int> > allele_to_reference(allele_reference.size());
    unsigned int n_alleles = 0;
    for (int i = 0; i < allele_reference.size(); i++) {
        if (n_alleles < allele_reference[i]) {
            n_alleles = allele_reference[i];
        }
        allele_to_reference[allele_reference[i]].push_back(i);
    }
    n_alleles++;
    allele_to_reference.resize(n_alleles);
    transition_probability_matrices.resize(n_alleles*n_alleles);
    reordering_map.resize(n_alleles*n_alleles);
    vector<unsigned int> paths_i;
    vector<unsigned int> paths_j;
    vector<unsigned int> previous_references;
    unsigned int matrix_index;
    unsigned int count;
    double value;
    vector<unsigned int>::iterator reference_i;
    vector<unsigned int>::iterator reference_j;
    
    // Takes n_haplotype_path^4 to make the transition matrix. Have to work on decreasing this as it affects scalability.
    for (int i = 0; i < n_alleles; i++) {
        paths_i = allele_to_reference[i];
        for (int j = 0; j < n_alleles; j++) {
            paths_j = allele_to_reference[j];
            matrix_index = i*n_alleles + j;
            transition_probability_matrices[matrix_index].remake(pow(allele_reference.size(),2), paths_i.size()*paths_j.size(), 0.0);
            for (unsigned int p_index = 0; p_index < pow(allele_reference.size(),2); p_index++) {
                count = 0;
                previous_references = column.index_to_reference_allele(p_index, 0);
                for (reference_i = paths_i.begin(); reference_i < paths_i.end(); reference_i++) {
                    for (reference_j = paths_j.begin(); reference_j < paths_j.end(); reference_j++) {
                        if ((*reference_i == previous_references[0]) && (*reference_j == previous_references[1])) value = qr*qr;
                        if ((*reference_i != previous_references[0]) && (*reference_j == previous_references[1])) value = pr*qr;
                        if ((*reference_i == previous_references[0]) && (*reference_j != previous_references[1])) value = pr*qr;
                        if ((*reference_i != previous_references[0]) && (*reference_j != previous_references[1])) value = pr*pr;
                        reordering_map[matrix_index].push_back(column.reference_allele_to_index(*reference_i, *reference_j));
                        transition_probability_matrices[matrix_index].set(p_index, count, value);
                        count++;
                    }
                }
            }
        }
    }
}

Vector2D<double> * TransitionProbabilityComputer::get_transition_matrix(unsigned int a1, unsigned int a2) {
    return &(this->transition_probability_matrices[a1*(int)(sqrt(transition_probability_matrices.size()))+a2]);
}

vector<unsigned int> * TransitionProbabilityComputer::get_reordering_map(unsigned int a1, unsigned int a2) {
    return &(this->reordering_map[a1*(int)(sqrt(reordering_map.size()))+a2]);
}
