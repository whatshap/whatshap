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

using namespace std;

GenotypeDPTable::GenotypeDPTable(ReadSet* read_set, const vector<unsigned int>& recombcost, const Pedigree* pedigree, const vector<unsigned int>* positions)
    :read_set(read_set),
     recombcost(recombcost),
     pedigree(pedigree),
     input_column_iterator(*read_set, positions)
{
   genotype_likelihood_table = Vector2D<genotype_likelihood_t>(pedigree->size(),input_column_iterator.get_column_count(),genotype_likelihood_t());
   read_set->reassignReadIds();

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

   // compute forward and backward probabilities
   //compute_backward_prob();
   compute_forward_prob();

}

GenotypeDPTable::~GenotypeDPTable()
{
    init(forward_projection_column_table,0);
    init(backward_projection_column_table, 0);
    init(indexers,0);
    init(pedigree_partitions,0);
}

//TODO seperate indexers for forward and backward pass?
void GenotypeDPTable::clear_forward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(forward_projection_column_table, column_count);
    init(indexers, column_count);
}

void GenotypeDPTable::clear_backward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(backward_projection_column_table, column_count);
    init(indexers, column_count);
}

unique_ptr<vector<unsigned int> > GenotypeDPTable::extract_read_ids(const vector<const Entry *>& entries) {
    unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
    for (size_t i=0; i<entries.size(); ++i) {
        read_ids->push_back(entries[i]->get_read_id());
    }
    return read_ids;
}

size_t GenotypeDPTable::popcount(size_t x) {
    unsigned int count = 0;
    for (;x; x >>= 1) {
        count += x & 1;
    }
    return count;
}

void GenotypeDPTable::compute_backward_prob()
{
    clear_backward_table();

    // if no reads are in the read set, nothing to do
    if(input_column_iterator.get_column_count() == 0){
        return;
    }

    // start at rightmost column

}

void GenotypeDPTable::compute_forward_prob()
{
    // the same code is used as in PegigreeDPTable::compute_table()...
    // first clear table
    clear_forward_table();

    // if no reads are in read set, nothing to compute
    if (input_column_iterator.get_column_count() == 0) {
        // TODO anything needs to be done here?
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
    ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0, *next_read_ids);
    indexers[0] = next_indexer;

    // forward pass: create a sparse table, storing values at every sqrt(#columns)-th position
    size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
    for (size_t column_index=0; column_index<input_column_iterator.get_column_count(); ++column_index) {
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);

        ColumnIndexingScheme* current_indexer = next_indexer;
        // peek ahead and get the next column TODO: what about last column? Anything needs to be done?
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
            next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
            current_indexer->set_next_column(next_indexer);
            indexers[column_index + 1] = next_indexer;
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
            next_indexer = 0;
        }

        // compute forward probabilities for the current column
        compute_forward_column(column_index,std::move(current_input_column));

        // check whether to delete the previous column
        if ((k>1) && (column_index > 0) && (((column_index-1)%k) != 0)) {
            delete forward_projection_column_table[column_index-1];
            forward_projection_column_table[column_index-1] = nullptr;
        }
    }
}

// given the current matrix column, compute the forward probability table
void GenotypeDPTable::compute_forward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{/** TODO for now, transition probabilities are assumed to be all the same, no recombination costs are
     considered.
       **/

    assert(column_index < input_column_iterator.get_column_count());

    // check whether requested column is already there, if so nothing to do
    if (forward_projection_column_table[column_index] != nullptr) return;

    ColumnIndexingScheme* current_indexer = indexers[column_index];
    assert(current_indexer != nullptr);

    // compute the number of different transmission vectors
    unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());

    // if the current input column was not provided, then create it
    if(current_input_column.get() == nullptr)
    {
        input_column_iterator.jump_to_column(column_index);
        current_input_column = input_column_iterator.get_next();
    }

    // reserve memory for the current DP column to be computed: number of bipartitions, transmission configurations
    // stores sums of previous costs for each bipartition and transmission configuration
    Vector2D<long double> dp_column(current_indexer->column_size(),transmission_configurations,0.0L);
    // keep track of sum of all alpha_i*beta_i
    long double normalization_sum = 0.0L;

    // obtain previous projection column (which is assumed to have already been computed)
    Vector2D<long double>* previous_projection_column = nullptr;
    if (column_index > 0) {
        previous_projection_column = forward_projection_column_table[column_index - 1];
    }

    // initialize the new projection column (2D: has entry for every bipartition and transmission value)
    Vector2D<long double>* current_projection_column = nullptr;
    if(column_index + 1 < input_column_iterator.get_column_count()){
        current_projection_column = new Vector2D<long double>(current_indexer->forward_projection_size(),transmission_configurations,0.0L);
    }

    // create column cost computer for each transmission vector
    vector<GenotypeColumnCostComputer> cost_computers;
    cost_computers.reserve(transmission_configurations);
    for(unsigned int i = 0; i < transmission_configurations; i++){
        cost_computers.emplace_back(*current_input_column, column_index, read_sources, pedigree, *pedigree_partitions[i]);
    }

    // for scaled version of forward backward algorithm, keep track of sum of all forward prob. in that column
    long double scaling_sum = 0.0L;

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

        // Determine index in the current DP column to be written
        size_t current_index = iterator->get_index();

        // iterate over all transmission vectors
        for(size_t i = 0; i < transmission_configurations; i++){
            // keep track of sum of previous values (alpha_i-1 * transition_prob)
            long double sum_prev_values = 0.0L;
            unsigned int number_of_allele_assignments = 1<<pedigree_partitions[i]->count();
            // TODO transition prob.!!!!used to normalize the transision probabilities, s.t. they sum up to 1
            if(column_index > 0){
                for(size_t j = 0; j < transmission_configurations; j++){
                    // TODO!!!
                    long double transition_prob = 1.0L;
                    // get cost from the previous column
                    long double previous_cost = previous_projection_column->at(backward_projection_index,j);
                    
                    // multiply with transition probability
                    previous_cost *= transition_prob;
                    sum_prev_values += previous_cost;
                }
            } else {
                // TODO times the 1/degree(start)
                sum_prev_values = 1.0L;
            }

            // store the sum of previous values in table
            dp_column.set(current_index,i,sum_prev_values);

            // iterate over all allele assignments
            for(unsigned int a = 0; a < number_of_allele_assignments; a++){
                // TODO get already computed backward probability from table
                long double backward_probability = 1.0L;
                long double forward_probability = dp_column.at(current_index,i)*cost_computers[i].get_cost(a);
                long double forward_backward = forward_probability * backward_probability;
                scaling_sum += forward_probability;
                normalization_sum += forward_backward;

                // marginalize over all genotypes
                for (size_t individuals_index = 0; individuals_index < pedigree->size(); individuals_index++) {
                    unsigned int partition0 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,0);
                    unsigned int partition1 = pedigree_partitions[i]->haplotype_to_partition(individuals_index,1);
                    unsigned int allele0 = (a >> partition0) & 1;
                    unsigned int allele1 = (a >> partition1) & 1;
                    unsigned int genotype = allele0 + allele1;
                    genotype_likelihood_table.at(individuals_index,column_index).likelihoods[genotype] += forward_backward;
                }

                // set forward projections
                if(current_projection_column != 0){
                    unsigned int forward_index = iterator->get_forward_projection();
                    current_projection_column->set(forward_index,i,forward_probability + current_projection_column->at(forward_index,i));
                }
            }
        }
    }

    // store the computed projection column (in case there is one)
    if(current_projection_column != 0){
        forward_projection_column_table[column_index] = current_projection_column;
    }

    // again go through projection column to scale the values (divide them by sum of forward prob.)
    if(current_projection_column != 0){
        for(size_t i = 0; i < current_indexer->forward_projection_size(); i++){
            for(size_t j = 0; j < transmission_configurations; j++){
                current_projection_column->at(i,j) /= scaling_sum;
            }
        }
    }
}

void GenotypeDPTable::compute_backward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{
    throw std::runtime_error("compute_backward_column() not yet implemented.");
}

vector<long double> GenotypeDPTable::get_genotype_likelihoods(unsigned int individual, unsigned int position)
{
    assert(individual < pedigree->size());
    assert(position < input_column_iterator.get_column_count());

    vector<long double> result;
    for(unsigned int i = 0; i < 3; i++){
        result.push_back(genotype_likelihood_table.at(individual,position).likelihoods[i]);
    }

    return result;
}

// TODO !!
long double GenotypeDPTable::compute_transition_prob(size_t t1, size_t t2, size_t length, unsigned int r){
    size_t x = t1 ^ t2;
    // count how many bits are set
    x = popcount(x);
    long double recomb_prob = pow(10,-(long double)(r)/10.0L);
    return pow(recomb_prob,x)*pow(1-recomb_prob,length-x);
}
