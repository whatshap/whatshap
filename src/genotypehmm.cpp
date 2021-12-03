#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "genotypecolumncostcomputer.h"
#include "genotypehmm.h"
#include "columnindexingiterator.h"
#include "transitionprobabilitycomputer.h"
#include "emissionprobabilitycomputer.h"
#include "binomial.h"
#include "matrixmultiplication.h"

using namespace std;

GenotypeHMM::GenotypeHMM(ReadSet* read_set, const vector<float>& recombcost, const Pedigree* pedigree, const unsigned int& n_references, const vector<unsigned int>* positions, const vector<unsigned int>* n_allele_positions,  const vector<vector<unsigned int> >* allele_references)
    :read_set(read_set),
     recombcost(recombcost),
     pedigree(pedigree),
     input_column_iterator(*read_set, positions),
     backward_input_column_iterator(*read_set, positions),
     transition_probability_table(input_column_iterator.get_column_count() - 1,nullptr),
     emission_probability_table(input_column_iterator.get_column_count(), nullptr),
     scaling_parameters(input_column_iterator.get_column_count(),-1.0L),
     variant_positions(positions),
     variant_n_allele_positions(n_allele_positions),
     n_references(n_references),
     allele_references(allele_references)
{
   
   genotype_likelihood_table = Vector2D<genotype_likelihood_t>(pedigree->size(),input_column_iterator.get_column_count());
   assert (pedigree->size() == 1);
   for (size_t i = 0; i < pedigree->size(); i ++) {
       for (size_t j = 0; j < input_column_iterator.get_column_count(); j++) {
           genotype_likelihood_table.set(i, j, genotype_likelihood_t(binomial_coefficient(n_allele_positions->at(j)+1, n_allele_positions->at(j)-1)));
       }
   }
   read_set->reassignReadIds();

   assert(input_column_iterator.get_column_count() == backward_input_column_iterator.get_column_count());

   // translate all individual ids to individual indices
   for(size_t i = 0; i<read_set->size(); ++i)
   {
       read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
   }

   //compute forward and backward probabilities
   compute_index();
   compute_backward_prob();
   compute_forward_prob();

}

GenotypeHMM::~GenotypeHMM()
{
    init(forward_pass_column_table,0);
    init(backward_pass_column_table, 0);
    init(hmm_columns,0);
    init(pedigree_partitions,0);
    init(transition_probability_table,0);
}

void GenotypeHMM::clear_forward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(forward_pass_column_table, 1);
}

void GenotypeHMM::clear_backward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(backward_pass_column_table, column_count);
}

unique_ptr<vector<unsigned int> > GenotypeHMM::extract_read_ids(const vector<const Entry *>& entries) {
    unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
    for (size_t i=0; i<entries.size(); ++i) {
        read_ids->push_back(entries[i]->get_read_id());
    }
    return read_ids;
}

void GenotypeHMM::compute_index(){
    size_t column_count = input_column_iterator.get_column_count();
    if(column_count == 0) return;
    init(hmm_columns, column_count);
    // do one forward pass to get the indexers (that are needed in forward and backward pass)
    input_column_iterator.jump_to_column(0);
    unique_ptr<vector<const Entry*> > current_input_column;
    unique_ptr<vector<const Entry*> > next_input_column;
    unique_ptr<vector<unsigned int> > current_read_ids;
    unique_ptr<vector<unsigned int> > next_read_ids;
    Column* current_column = nullptr;
    next_input_column = input_column_iterator.get_next();
    next_read_ids = extract_read_ids(*current_input_column);
    
    for(size_t column_index=0; column_index < input_column_iterator.get_column_count(); ++column_index){
        
        current_input_column = std::move(next_input_column);
        current_read_ids = std::move(next_read_ids);
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
            current_column = new Column(column_index, &n_references, *current_read_ids, *next_read_ids);
            hmm_columns[column_index] = current_column;
            transition_probability_table[column_index] = new TransitionProbabilityComputer(recombcost[column_index], *current_column, allele_references->at(column_index+1));
        } 
        else {
            assert (column_index == input_column_iterator.get_column_count() - 1);
            current_column = new Column(column_index, &n_references, *current_read_ids, vector<unsigned int>{}); 
            hmm_columns[column_index] = current_column;
        }
    }
}

void GenotypeHMM::compute_backward_prob()
{
    clear_backward_table();
    unsigned int column_count = backward_input_column_iterator.get_column_count();

    // if no reads are in the read set, nothing to do
    if(backward_input_column_iterator.get_column_count() == 0){
        return;
    }

    // do backward pass, start at rightmost column
    backward_input_column_iterator.jump_to_column(column_count-1);

    // get the next column (which is left of current one)
    unique_ptr<vector<const Entry*> > current_input_column;
    unique_ptr<vector<unsigned int> > current_read_ids;
    unique_ptr<vector<const Entry*> > next_input_column = backward_input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);

    // backward pass: create sparse table
    size_t k = (size_t)sqrt(column_count);
    for(int column_index = column_count-1; column_index >= 0; --column_index){
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        current_read_ids = std::move(next_read_ids);
        // peek ahead and get the next column
        if (backward_input_column_iterator.has_next()){
            next_input_column = backward_input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
        }

        // compute the backward probabilities
        compute_backward_column(column_index, std::move(current_input_column));

        // check whether to delete the previous column
        if ((k>1) && (column_index < column_count-1) && (((column_index+1)%k) != 0)) {
            delete backward_pass_column_table[column_index+1];
            backward_pass_column_table[column_index+1] = nullptr;
        }
    }
}

void GenotypeHMM::compute_forward_prob()
{
    clear_forward_table();

    // if no reads are in read set, nothing to compute
    if (input_column_iterator.get_column_count() == 0) {
        return;
    }

    // start at leftmost column (= 0th column)
    input_column_iterator.jump_to_column(0);
    // store current and next column
    unique_ptr<vector<const Entry *> > current_input_column;
    unique_ptr<vector<const Entry *> > next_input_column;
    // get the next column ahead of time
    next_input_column = input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);

    // forward pass: create a sparse table, storing values at every sqrt(#columns)-th position
    for (size_t column_index=0; column_index < input_column_iterator.get_column_count(); ++column_index) {
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
        // peek ahead and get the next column
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
        }

        // compute forward probabilities for the current column
        compute_forward_column(column_index,std::move(current_input_column));
    }
}

void GenotypeHMM::compute_backward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column) {
        
    assert(column_index < backward_input_column_iterator.get_column_count());

    // check if column already exists
    if(column_index > 0){
        if (backward_pass_column_table[column_index-1] != nullptr) return;
    }

    Column* current_indexer = hmm_columns[column_index];
    TransitionProbabilityComputer* current_transition_table = transition_probability_table[column_index-1];
    assert(current_indexer != nullptr);
    assert (current_transition_table != nullptr);

    // if current input column was not provided, create it
    if(current_input_column.get() == nullptr) {
        backward_input_column_iterator.jump_to_column(column_index);
        current_input_column = backward_input_column_iterator.get_next();
    }

    vector<long double>* previous_projection_column = nullptr;

    // check if there is a projection column
    if(column_index < backward_input_column_iterator.get_column_count()-1){
        previous_projection_column = backward_pass_column_table[column_index];
    }

    // initialize the new projection column (= current index-1)
    vector<long double>* current_projection_column = nullptr;
    if(column_index > 0){
        current_projection_column = new vector<long double>(hmm_columns[column_index-1]->get_column_size(), 0.0L);
    }
    else {
        return;
    }

    int n_alleles = variant_n_allele_positions->at(column_index);
    Vector2D<long double> emission_probability_computer = Vector2D<long double>(n_alleles, n_alleles);

    // for scaled version of forward backward alg, keep track of the sum of backward
    long double scaling_sum = 0.0L;

    // iterate over all bipartitions of the column on the right. So we are calculating the values in column current_index - 1
    
    unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
    while (iterator->has_next()){
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        // Update the emission probability based on the bipartition defined by the iterator
        update_emission_probability(&emission_probability_computer, bit_changed, *iterator, *current_input_column);

        // Determine the indices that are compatible with with the bipartition at position column_index
        long double backward_prob = 1.0L;
        int b_index = iterator->get_b_index();
        vector<unsigned int> compatible_bipartitions = hmm_columns[column_index-1]->get_backward_compatible_bipartitions(b_index);

        // UPDATING THE BETA VALUES

        // Iterating over the allele pairs
        for (int allele_1 = 0; allele_1 < n_alleles; allele_1++) {
            for (int allele_2 = 0; allele_2 < n_alleles; allele_2++) {
                // Extract the beta values having (allele_1, allele_2) pair in proper ordering (based on the transition probability matrix)
                vector<long double> beta_values;
                Vector2D<double> * transition_matrix = current_transition_table->get_transition_matrix(allele_1, allele_2);
                vector<unsigned int> * reordering_map = current_transition_table->get_reordering_map(allele_1, allele_2);
                for (int i = 0; i < reordering_map->size(); i++) {
                    int r_index = reordering_map->at(i);
                    int index = current_indexer->get_index(b_index, r_index);
                    if (column_index + 1 < backward_input_column_iterator.get_column_count()) {
                        beta_values.push_back(previous_projection_column->at(index));
                        scaling_sum += previous_projection_column->at(index);
                    }
                    else {
                        beta_values.push_back(1.0L);
                        scaling_sum += 1.0L;
                    }
                }
                vector<long double> result = MatMul(beta_values, *transition_matrix, true);
                // Update the beta values in the current_projection_matrix
                for (int i = 0; i < compatible_bipartitions.size(); i++) {
                    int index = compatible_bipartitions[i] * pow(n_references,2);
                    for (unsigned int j = 0; j < pow(n_references, 2); j++) {
                        current_projection_column->at(index+j) = current_projection_column->at(index+j) + result[j]*emission_probability_computer.at(allele_1, allele_2);   
                    }
                }
            }
        }
    }
    
    // go through (old) projection column to scale the values -> when we lookup betas later, they will sum up to 1
    if(previous_projection_column != 0){
        std::transform((*previous_projection_column).begin(), (*previous_projection_column).end(), (*previous_projection_column).begin(), std::bind2nd(std::divides<long double>(), scaling_sum));
    }
    if(current_projection_column != 0){
        std::transform((*current_projection_column).begin(), (*current_projection_column).end(), (*current_projection_column).begin(), std::bind2nd(std::divides<long double>(), scaling_sum));
        backward_pass_column_table[column_index-1] = current_projection_column;
    }
    scaling_parameters[column_index] = scaling_sum;
}

// given the current matrix column, compute the forward probability table
void GenotypeHMM::compute_forward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{
    assert(column_index < input_column_iterator.get_column_count());

    Column* current_indexer = hmm_columns[column_index];
    assert(current_indexer != nullptr);

    // if the current input column was not provided, then create it
    if(current_input_column.get() == nullptr) {
        input_column_iterator.jump_to_column(column_index);
        current_input_column = input_column_iterator.get_next();
    }

    // obtain previous projection column (which is assumed to have already been computed)
    vector<long double>* previous_projection_column = nullptr;
    if (column_index > 0) {
        // forward_pass_column_table is a vector (containing vectors) of size 1 (since we only need to store one column at a time). Jana plans on changing it to a vector containing long double values.
        previous_projection_column = forward_pass_column_table[0];
        assert(previous_projection_column != nullptr);
    }

    // obtain the backward projection table, from where to get the backward probabilities
    // ASSUMED THAT THIS WORKS
    size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
    vector<long double>* backward_probabilities = nullptr;
    if(column_index + 1 < input_column_iterator.get_column_count()){
        backward_probabilities = backward_pass_column_table[column_index];
        // if column is not stored, recompute it
        if(backward_probabilities == nullptr){
            // compute index of next column that has been stored
            size_t next = std::min((unsigned int) ( ((column_index + k) / k) * k ), input_column_iterator.get_column_count()-1);
            for(size_t i = next; i > column_index; --i){
                compute_backward_column(i);
            }
            // last column just computed still needs to be scaled
            std::transform((*backward_pass_column_table[column_index]).begin(), (*backward_pass_column_table[column_index]).end(), (*backward_pass_column_table[column_index]).begin(), std::bind2nd(std::divides<long double>(), scaling_parameters[column_index]));
        }
        backward_probabilities = backward_pass_column_table[column_index];
        assert(backward_probabilities != nullptr);
    }

    // initialize the new projection column (2D: has entry for every bipartition and transmission value)
    vector<long double>* current_projection_column = nullptr;
    // CHECK!!!!!
    if(column_index + 1 < input_column_iterator.get_column_count()){
        current_projection_column = new vector<long double>(current_indexer->get_column_size(), 0.0L);
    }

    int n_alleles = variant_n_allele_positions->at(column_index);
    Vector2D<long double> emission_probability_computer = Vector2D<long double>(n_alleles, n_alleles);
    TransitionProbabilityComputer* current_transition_table;
    if (column_index > 0) {
        current_transition_table = transition_probability_table[column_index - 1];
    }

    // sum of alpha*beta, used to normalize the likelihoods
    long double normalization = 0.0L;

    // iterate over all bipartitions
    unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
    while (iterator->has_next()) {
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        // Update the emission probability based on the bipartition defined by the iterator
        update_emission_probability(&emission_probability_computer, bit_changed, *iterator, *current_input_column);

        // Determine the indices that are compatible with with the bipartition at position column_index
        long double backward_prob = 1.0L;
        int b_index = iterator->get_b_index();
        // Calculating the current_projection_column
        if (column_index == 0) {
            long double tr_prb = 1/current_indexer->get_column_size();
            for (unsigned int i = 0; i < allele_references->at(column_index).size(); i++) {
                int allele_1 = allele_references->at(column_index).at(i);
                for (unsigned int j = 0; j < allele_references->at(column_index).size(); j++) {
                    int allele_2 = allele_references->at(column_index).at(j);
                    unsigned int r_index = current_indexer->reference_allele_to_index(i, j);
                    unsigned int index = current_indexer->get_index(b_index, r_index);
                    current_projection_column->at(index) = tr_prb * emission_probability_computer.at(allele_1, allele_2);
                }
            }
        }
        else {
            vector<unsigned int> compatible_bipartitions = hmm_columns[column_index-1]->get_backward_compatible_bipartitions(b_index);
            for (unsigned int b : compatible_bipartitions) {
                // Get the alpha values from the compatible bipartition b
                vector<long double> alpha_values;
                for (int r_index = 0; r_index < pow(n_references, 2); r_index++) {
                    alpha_values.push_back(previous_projection_column->at((n_references*n_references*b)+r_index));
                }
                // Iterating over the groups to calculate the alpha values.
                for (int allele_1 = 0; allele_1 < n_alleles; allele_1++) {
                    for (int allele_2 = 0; allele_2 < n_alleles; allele_2++) {
                        Vector2D<double> * transition_matrix = current_transition_table->get_transition_matrix(allele_1, allele_2);
                        vector<unsigned int> * reordering_map = current_transition_table->get_reordering_map(allele_1, allele_2);
                        vector<long double> result = MatMul(alpha_values, *transition_matrix);
                        // Taking the alpha values from the result array and then mapping them to the current_projection_column based on reordering_map.
                        for (unsigned int pos_index = 0; pos_index < reordering_map->size(); pos_index++) {
                            unsigned int r_index = reordering_map->at(pos_index);
                            current_projection_column->at((n_references * n_references * b_index) + r_index) += result.at(pos_index)*emission_probability_computer.at(allele_1, allele_2);
                        }
                    }
                }
            }
        }
        long double forward_backward;
        long double forward;
        assert (current_projection_column->size() == backward_probabilities->size());

        for (unsigned int i = 0; i < backward_probabilities->size(); i++) {
            vector<unsigned int> ref = current_indexer->index_to_reference_allele(i, 0);
            vector<unsigned int> alleles;
            alleles.push_back(allele_references->at(column_index).at(ref.at(0)));
            alleles.push_back(allele_references->at(column_index).at(ref.at(1)));
            // Get the genotype index
            unsigned int g_index = 0;
            for (int allele = 0; allele < 2; allele++) {
                g_index += binomial_coefficient(allele + alleles.at(allele), alleles.at(allele) - 1);
            }
            forward = current_projection_column->at(i) / scaling_parameters[column_index];
            forward_backward = forward * backward_probabilities->at(i);
            
            // HARDCODED FOR A PEDIGREE SIZE OF 1.
            genotype_likelihood_table.at(0, column_index).likelihoods[g_index] += forward_backward;
        }
    }
    // store the computed projection column (in case there is one)
    if(current_projection_column != 0){
        delete forward_pass_column_table[0];
        forward_pass_column_table[0] = current_projection_column;
    }

    // we can remove the backward-probability column
    if(backward_pass_column_table[column_index] != nullptr){
        delete backward_pass_column_table[column_index];
        backward_pass_column_table[column_index] = nullptr;
    }

    // scale the likelihoods
    for(size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index){
        genotype_likelihood_table.at(individuals_index,column_index).divide_likelihoods_by(normalization);
    }
}

vector<long double> GenotypeHMM::get_genotype_likelihoods(unsigned int individual_id, unsigned int position)
{
    assert(pedigree->id_to_index(individual_id) < genotype_likelihood_table.get_size0());
    assert(position < input_column_iterator.get_column_count());

    return genotype_likelihood_table.at(pedigree->id_to_index(individual_id),position).likelihoods;

}

void GenotypeHMM::update_emission_probability(Vector2D<long double>* em_prob, const int bit_changed, ColumnIndexingIterator iterator, vector<const Entry *> entries) {
    int n_alleles = em_prob->get_size0();
    if (bit_changed >= 0) {
        int newBit = iterator.get_binary_vector()[bit_changed];
        if (entries.at(bit_changed)->get_allele_type() == -1) {
            return;
        }
        if (newBit == 1) {
            for (int i = 0; i < n_alleles; i++) {
                for (int j = 0; j < n_alleles; j++) {
                    em_prob->set(i, j, em_prob->at(i, j) * (entries.at(bit_changed)->get_emission_score()[j]/entries.at(bit_changed)->get_emission_score()[i]));
                }
            }
        }
        else {
            for (int i = 0; i < n_alleles; i++) {
                for (int j = 0; j < n_alleles; j++) {
                    em_prob->set(i, j, em_prob->at(i, j) * (entries.at(bit_changed)->get_emission_score()[i]/entries.at(bit_changed)->get_emission_score()[j]));
                }
            }
        }
    }
    else {
        for (int i = 0; i < n_alleles; i++) {
            long double value = 1.0;
            for (int entry_index = 0; entry_index < entries.size(); entry_index++) {
                if (entries.at(entry_index)->get_allele_type() == -1) {
                    continue;
                }
                value = value * (entries.at(entry_index)->get_emission_score()[i]);
            }
            for (int j = 0; j < n_alleles; j++) {
                em_prob->set(i, j, value);
            }
        }
    }
}