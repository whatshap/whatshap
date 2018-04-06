#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "genotypecolumncostcomputer.h"
#include "genotypedptable.h"
#include "columnindexingiterator.h"
#include "transitionprobabilitycomputer.h"

using namespace std;

GenotypeDPTable::GenotypeDPTable(ReadSet* read_set, const vector<unsigned int>& recombcost, const Pedigree* pedigree, const vector<unsigned int>* positions)
    :read_set(read_set),
     recombcost(recombcost),
     pedigree(pedigree),
     input_column_iterator(*read_set, positions),
     backward_input_column_iterator(*read_set, positions),
     transition_probability_table(input_column_iterator.get_column_count(),nullptr),
     scaling_parameters(input_column_iterator.get_column_count(),-1.0L)
{
   genotype_likelihood_table = Vector2D<genotype_likelihood_t>(pedigree->size(),input_column_iterator.get_column_count(),genotype_likelihood_t());
   read_set->reassignReadIds();

   assert(input_column_iterator.get_column_count() == backward_input_column_iterator.get_column_count());

   // create all pedigree partitions
   for(size_t i = 0; i < std::pow(4,pedigree->triple_count()); ++i)
   {
       pedigree_partitions.push_back(new PedigreePartitions(*pedigree,i));
   }

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

GenotypeDPTable::~GenotypeDPTable()
{
    init(forward_projection_column_table,0);
    init(backward_projection_column_table, 0);
    init(indexers,0);
    init(pedigree_partitions,0);
    init(transition_probability_table,0);
}

void GenotypeDPTable::clear_forward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(forward_projection_column_table, 1);
}

void GenotypeDPTable::clear_backward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(backward_projection_column_table, column_count);
}

unique_ptr<vector<unsigned int> > GenotypeDPTable::extract_read_ids(const vector<const Entry *>& entries) {
    unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
    for (size_t i=0; i<entries.size(); ++i) {
        read_ids->push_back(entries[i]->get_read_id());
    }
    return read_ids;
}

void GenotypeDPTable::compute_index(){
    size_t column_count = input_column_iterator.get_column_count();
    if(column_count == 0) return;
    init(indexers, column_count);
    // do one forward pass to get the indexers (that are needed in forward and backward pass)
    input_column_iterator.jump_to_column(0);
    unique_ptr<vector<const Entry*> > current_input_column;
    unique_ptr<vector<const Entry*> > next_input_column;
    next_input_column = input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
    ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0, *next_read_ids);
    indexers[0] = next_indexer;
    unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());
    transition_probability_table[0] = new TransitionProbabilityComputer(0, recombcost[0], pedigree, pedigree_partitions);

    for(size_t column_index=0; column_index < input_column_iterator.get_column_count(); ++column_index){
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
        ColumnIndexingScheme* current_indexer = next_indexer;

        if (input_column_iterator.has_next()){
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
            next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);

            current_indexer->set_next_column(next_indexer);
            indexers[column_index+1] = next_indexer;
            transition_probability_table[column_index+1] = new TransitionProbabilityComputer(column_index+1, recombcost[column_index+1], pedigree, pedigree_partitions);
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
            next_indexer = 0;
        }
    }
}

void GenotypeDPTable::compute_backward_prob()
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
    unique_ptr<vector<const Entry*> > next_input_column = backward_input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);

    // backward pass: create sparse table
    size_t k = (size_t)sqrt(column_count);
    for(int column_index=column_count-1; column_index >= 0; --column_index){
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
        // peek ahead and get the next column
        if (backward_input_column_iterator.has_next()){
            next_input_column = backward_input_column_iterator.get_next();
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
        }

        // compute the backward probabilities
        compute_backward_column(column_index,std::move(current_input_column));

        // check whether to delete the previous column
        if ((k>1) && (column_index < column_count-1) && (((column_index+1)%k) != 0)) {
            delete backward_projection_column_table[column_index+1];
            backward_projection_column_table[column_index+1] = nullptr;
        }
    }
}

void GenotypeDPTable::compute_forward_prob()
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
    size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
    for (size_t column_index=0; column_index<input_column_iterator.get_column_count(); ++column_index) {
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

void GenotypeDPTable::compute_backward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{
   assert(column_index < backward_input_column_iterator.get_column_count());

   // check if column already exists
   if(column_index > 0){
       if (backward_projection_column_table[column_index-1] != nullptr) return;
   }

   ColumnIndexingScheme* current_indexer = indexers[column_index];
   assert(current_indexer != nullptr);

   // number of transmission values
   unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());

   // if current input column was not provided, create it
   if(current_input_column.get() == nullptr) {
       backward_input_column_iterator.jump_to_column(column_index);
       current_input_column = backward_input_column_iterator.get_next();
   }

   // obtain previous projection column (same index as current column!)
   Vector2D<long double>* previous_projection_column = nullptr;
   // check if there is a projection column
   if(column_index < backward_input_column_iterator.get_column_count()-1){
       previous_projection_column = backward_projection_column_table[column_index];
   }

   // initialize the new projection column (= current index -1)
   Vector2D<long double>* current_projection_column = nullptr;
   if(column_index > 0){
       current_projection_column = new Vector2D<long double>(indexers[column_index-1]->forward_projection_size(),transmission_configurations,0.0L);
   }

   // create column cost computer for each transmission vector
   vector<GenotypeColumnCostComputer> cost_computers;
   cost_computers.reserve(transmission_configurations);
   for(unsigned int i = 0; i < transmission_configurations; ++i){
       cost_computers.emplace_back(*current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[i]);
   }

   // for scaled version of forward backward alg, keep track of the sum of backward
   long double scaling_sum = 0.0L;

   // iterate over all bipartitions
   unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
   while (iterator->has_next()){
       int bit_changed = -1;
       iterator->advance(&bit_changed);
       if (bit_changed >= 0) {
           // update bipartition in the cost computers
           for(auto& cost_computer : cost_computers) {
               cost_computer.update_partitioning(bit_changed);
           }
       } else {
           // set bipartition in the cost computers
           for(auto& cost_computer : cost_computers) {
              cost_computer.set_partitioning(iterator->get_partition());
           }
       }

       // Determine index in forward projection column from where to fetch the current cost
       long double backward_prob = 1.0L;

       // iterate over all transmission configurations
       for(size_t i = 0; i < transmission_configurations; ++i){
           // number of allele assignments
           unsigned int number_of_allele_assignments = 1<<pedigree_partitions[i]->count();

           // get entry from forward projection column (which is equal to current backward prob. for all genotypes)
           size_t forward_projection_index = 0;
           if (column_index + 1 < backward_input_column_iterator.get_column_count()) {
               forward_projection_index = iterator->get_forward_projection();
               backward_prob = previous_projection_column->at(forward_projection_index,i);
           }

           size_t backward_projection_index = iterator->get_backward_projection();
           // sum up entries in backward projection column
           for(unsigned int a = 0; a < number_of_allele_assignments; ++a){
               if(column_index > 0){
                   long double local_cost = cost_computers[i].get_cost(a);
                   for(size_t j = 0; j < transmission_configurations; ++j){
                       long double transition_prob = transition_probability_table[column_index]->get_prob_transmission(j,i) * transition_probability_table[column_index]->get_prob_allele_assignment(i,a);
                       current_projection_column->at(backward_projection_index, j) += backward_prob * local_cost * transition_prob;
                   }
               }
               scaling_sum += backward_prob;
           }
       }
   }

   // go through (old) projection column to scale the values -> when we lookup betas later, they will sum up to 1
   if(previous_projection_column != 0){
       previous_projection_column->divide_entries_by(scaling_sum);
   }
   if(current_projection_column != 0){
       current_projection_column->divide_entries_by(scaling_sum);
       backward_projection_column_table[column_index-1] = current_projection_column;
   }
   scaling_parameters[column_index] = scaling_sum;
}

// given the current matrix column, compute the forward probability table
void GenotypeDPTable::compute_forward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{
    assert(column_index < input_column_iterator.get_column_count());

    ColumnIndexingScheme* current_indexer = indexers[column_index];
    assert(current_indexer != nullptr);

    // compute the number of different transmission vectors
    unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());

    // if the current input column was not provided, then create it
    if(current_input_column.get() == nullptr) {
        input_column_iterator.jump_to_column(column_index);
        current_input_column = input_column_iterator.get_next();
    }

    // obtain previous projection column (which is assumed to have already been computed)
    Vector2D<long double>* previous_projection_column = nullptr;
    if (column_index > 0) {
        previous_projection_column = forward_projection_column_table[0];
        assert(previous_projection_column != nullptr);
    }

    // obtain the backward projection table, from where to get the backward probabilities
    size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
    Vector2D<long double>* backward_probabilities = nullptr;
    if(column_index + 1 < input_column_iterator.get_column_count()){
        backward_probabilities = backward_projection_column_table[column_index];
        // if column is not stored, recompute it
        if(backward_probabilities == nullptr){
            // compute index of next column that has been stored
            size_t next = std::min((unsigned int) ( ((column_index + k) / k) * k ), input_column_iterator.get_column_count()-1);
            for(size_t i = next; i > column_index; --i){
                compute_backward_column(i);
            }

            // last column just computed still needs to be scaled
            backward_projection_column_table[column_index]->divide_entries_by(scaling_parameters[column_index]);
        }
        backward_probabilities = backward_projection_column_table[column_index];
        assert(backward_probabilities != nullptr);
    }

    // initialize the new projection column (2D: has entry for every bipartition and transmission value)
    Vector2D<long double>* current_projection_column = nullptr;
    if(column_index + 1 < input_column_iterator.get_column_count()){
        current_projection_column = new Vector2D<long double>(current_indexer->forward_projection_size(),transmission_configurations,0.0L);
    }

    // create column cost computer for each transmission vector
    vector<GenotypeColumnCostComputer> cost_computers;
    cost_computers.reserve(transmission_configurations);
    for(unsigned int i = 0; i < transmission_configurations; ++i){
        cost_computers.emplace_back(*current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[i]);
    }

    // sum of alpha*beta, used to normalize the likelihoods
    long double normalization = 0.0L;

    // iterate over all bipartitions
    unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
    while (iterator->has_next()) {
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        if (bit_changed >= 0) {
            // update bipartition in the cost computers
            for(auto& cost_computer : cost_computers) {
                cost_computer.update_partitioning(bit_changed);
            }
        } else {
            // set bipartition in the cost computers
            for(auto& cost_computer : cost_computers) {
               cost_computer.set_partitioning(iterator->get_partition());
            }
        }

        // Determine index in backward projection column from where to fetch the previous cost
        size_t backward_projection_index = 0;
        if (column_index > 0) {
            backward_projection_index = iterator->get_backward_projection();
        }

        // iterate over all transmission vectors
        for(size_t i = 0; i < transmission_configurations; ++i){
            // keep track of sum of previous values (alpha_i-1 * transition_prob)
            long double sum_prev_values = 0.0L;
            unsigned int number_of_allele_assignments = 1<<pedigree_partitions[i]->count();
            if(column_index > 0){
                for(size_t j = 0; j < transmission_configurations; ++j){
                    // add product of previous cost * transition_probability
                    sum_prev_values += ((previous_projection_column->at(backward_projection_index,j)) * (transition_probability_table[column_index]->get_prob_transmission(j,i)) );
                }
            } else {
                sum_prev_values = 1.0L;
            }

            // iterate over all allele assignments
            for(unsigned int a = 0; a < number_of_allele_assignments; ++a){
                // get already computed backward probability from table
                long double backward_probability = 1.0L;
                if(backward_probabilities != nullptr){
                    size_t forward_projection_index = iterator->get_forward_projection();
                    backward_probability = backward_probabilities->at(forward_projection_index,i);
                }

                long double forward_probability = ( sum_prev_values * cost_computers[i].get_cost(a) * (transition_probability_table[column_index]->get_prob_allele_assignment(i,a)) ) / scaling_parameters[column_index];
                long double forward_backward = forward_probability * backward_probability;
                normalization += forward_backward;

                // marginalize over all genotypes
                for (size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index) {
                    unsigned int partition0 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,0);
                    unsigned int partition1 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,1);
                    unsigned int allele0 = (a >> partition0) & 1;
                    unsigned int allele1 = (a >> partition1) & 1;
                    genotype_likelihood_table.at(individuals_index,column_index).likelihoods[allele0 + allele1] += forward_backward;
                }

                // set forward projections
                if(current_projection_column != 0){
                    current_projection_column->at( iterator->get_forward_projection() ,i) += forward_probability;
                }
            }
        }
    }

    // store the computed projection column (in case there is one)
    if(current_projection_column != 0){
        delete forward_projection_column_table[0];
        forward_projection_column_table[0] = current_projection_column;
    }

    // we can remove the backward-probability column
    if(backward_projection_column_table[column_index] != nullptr){
        delete backward_projection_column_table[column_index];
        backward_projection_column_table[column_index] = nullptr;
    }

    // scale the likelihoods
    for(size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index){
        genotype_likelihood_table.at(individuals_index,column_index).divide_likelihoods_by(normalization);
    }
}

vector<long double> GenotypeDPTable::get_genotype_likelihoods(unsigned int individual_id, unsigned int position)
{
    assert(pedigree->id_to_index(individual_id) < genotype_likelihood_table.get_size0());
    assert(position < input_column_iterator.get_column_count());

    return genotype_likelihood_table.at(pedigree->id_to_index(individual_id),position).likelihoods;

}
