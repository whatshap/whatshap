#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "pedigreecolumncostcomputer.h"
#include "pedigreedptable.h"

using namespace std;

PedigreeDPTable::PedigreeDPTable(ReadSet* read_set, const vector<unsigned int>& recombcost, const Pedigree* pedigree)
	:
	read_set(read_set),
	recombcost(recombcost),
	pedigree(pedigree),
	optimal_score(0u),
	optimal_score_index(0u),
	index_backtrace_table(),
	transmission_backtrace_table()
{
	read_set->reassignReadIds();

	for (size_t i=0; i<std::pow(4, pedigree->triple_count()); ++i) {
		pedigree_partitions.push_back(new PedigreePartitions(*pedigree, i));
	}

	// translate all individual ids to individual indices
	for (size_t i=0; i<read_set->size(); ++i) {
		read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
	}

	compute_table();
}

PedigreeDPTable::~PedigreeDPTable() {
	for(size_t i=0; i<pedigree_partitions.size(); ++i) {
		delete pedigree_partitions[i];
	}

	for(size_t i=0; i<indexers.size(); ++i) {
		delete indexers[i];
	}

	for(size_t i=0; i<index_backtrace_table.size(); ++i) {
		delete index_backtrace_table[i];
	}

	for(size_t i=0; i<transmission_backtrace_table.size(); ++i) {
		delete transmission_backtrace_table[i];
	}
}

unique_ptr<vector<unsigned int> > PedigreeDPTable::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
	for (size_t i=0; i<entries.size(); ++i) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}


size_t PedigreeDPTable::popcount(size_t x) {
	unsigned int count = 0;
	for (;x; x >>= 1) {
		count += x & 1;
	}
	return count;
}


void PedigreeDPTable::clear_table() {
	if(!indexers.empty()) { // clear indexers, if present
		for(size_t i=0; i<indexers.size(); ++i) {
			delete indexers[i];
		}
		indexers.resize(0);
	}

	if(!index_backtrace_table.empty()) { // clear backtrace_table, if present
		for(size_t i=0; i<index_backtrace_table.size(); ++i) {
			delete index_backtrace_table[i];
		}
		index_backtrace_table.resize(0);
	}

	if(!transmission_backtrace_table.empty()) { // clear backtrace_table, if present
		for(size_t i=0; i<transmission_backtrace_table.size(); ++i) {
			delete transmission_backtrace_table[i];
		}
		transmission_backtrace_table.resize(0);
	}
	optimal_score = 0;
	optimal_score_index = 0;
	optimal_transmission_value = 0;
	previous_transmission_value = 0;
}


void PedigreeDPTable::compute_table() {
	// Compute the number of different transmission vectors
	unsigned int transmission_configurations = std::pow(4, pedigree->triple_count());
	ColumnIterator input_column_iterator(*read_set);

	clear_table();

	// empty read-set, nothing to phase, so MEC score is 0
	if (!input_column_iterator.has_next()) {
		optimal_score = 0;
		optimal_score_index = 0;
		return;
	}
	
	unsigned int n = 0;
	unique_ptr<vector<const Entry *> > current_input_column;
	unique_ptr<vector<const Entry *> > next_input_column;
	// get the next column ahead of time
	next_input_column = input_column_iterator.get_next();
	unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
	ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0,*next_read_ids);
	indexers.push_back(next_indexer);
	unique_ptr<Vector2D<unsigned int>> previous_projection_column;
	unique_ptr<Vector2D<unsigned int>> current_projection_column;
	optimal_score = numeric_limits<unsigned int>::max();
	optimal_score_index = 0;
	optimal_transmission_value = 0;
	previous_transmission_value = 0;

	while(next_indexer != 0) {
		// move on projection column
		previous_projection_column = std::move(current_projection_column);
		// make former next column the current one
		current_input_column = std::move(next_input_column);
		unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);

		ColumnIndexingScheme* current_indexer = next_indexer;
		// peek ahead and get the next column
		if (input_column_iterator.has_next()) {
			next_input_column = input_column_iterator.get_next();
			next_read_ids = extract_read_ids(*next_input_column);
			next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
			current_indexer->set_next_column(next_indexer);
			indexers.push_back(next_indexer);
		} else {
			assert(next_input_column.get() == 0);
			assert(next_read_ids.get() == 0);
			next_indexer = 0;
		}
		// reserve memory for the DP column
		Vector2D<unsigned int> dp_column(current_indexer->column_size(), transmission_configurations, 0);
		Vector2D<unsigned int>* transmission_backtrace_column = nullptr;
		Vector2D<unsigned int>* index_backtrace_column = nullptr;
		// if not last column, reserve memory for forward projections column
		if (next_input_column.get() != 0) {
			current_projection_column = unique_ptr<Vector2D<unsigned int>>(
				new Vector2D<unsigned int>(current_indexer->forward_projection_size(), transmission_configurations, numeric_limits<unsigned int>::max())
			);
			// NOTE: forward projection size will always be even
			transmission_backtrace_column = new Vector2D<unsigned int>(current_indexer->forward_projection_size(), transmission_configurations, numeric_limits<unsigned int>::max());
			index_backtrace_column = new Vector2D<unsigned int>(current_indexer->forward_projection_size(), transmission_configurations, numeric_limits<unsigned int>::max());
		}

		// do the actual compution on current column
		vector<PedigreeColumnCostComputer> cost_computers;
		cost_computers.reserve(transmission_configurations);
		for(unsigned int i = 0; i < transmission_configurations; ++i) {
			cost_computers.emplace_back(*current_input_column, n, read_sources, pedigree, *pedigree_partitions[i]);
		}

		unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();

		while (iterator->has_next()) {
			int bit_changed = -1;
			iterator->advance(&bit_changed);
			if (bit_changed >= 0) {
				for(auto& cost_computer : cost_computers) {
					cost_computer.update_partitioning(bit_changed);
				}
			} else {
				for(auto& cost_computer : cost_computers) {
					cost_computer.set_partitioning(iterator->get_partition());
				}
			}

			// Determine index in backward projection column from were to fetch the previous cost
			size_t backward_projection_index = iterator->get_backward_projection();
			// Determine index in the current DP column to be written
			size_t current_index = iterator->get_index();

			// Compute aggregate cost based on cost in previous and cost in current column
			vector<unsigned int> min_recomb_index(transmission_configurations);
			bool found_valid_transmission_vector = false;
			for (size_t i = 0; i < transmission_configurations; ++i) {
				// Compute cost incurred by current cell of DP table
				unsigned int current_cost = cost_computers[i].get_cost();
				unsigned int min = numeric_limits<unsigned int>::max();
				size_t min_index = 0;
				if (current_cost < numeric_limits<unsigned int>::max()) {
					found_valid_transmission_vector = true;
				}
				for (size_t j = 0; j < transmission_configurations; ++j) {
					// Step 1: add up cost from current_cost column and previous columns
					unsigned int val;
					unsigned int previous_cost = 0;
					if (previous_projection_column.get() != nullptr) {
						previous_cost = previous_projection_column->at(backward_projection_index,j);
					}
					if ((current_cost < numeric_limits<unsigned int>::max()) && (previous_cost < numeric_limits<unsigned int>::max())) {
						val = current_cost + previous_cost;
					} else {
						val = numeric_limits<unsigned int>::max();
					}
					// Step 2: add further cost incurred by recombination
					// change in bit 0 --> recombination in mother
					size_t x = i ^ j; // count the number of bits set in x

					if (val < numeric_limits<unsigned int>::max()) {
						val += popcount(x) * recombcost[n];
					}
					
					// check for new minimum
					if (val < min) {
						min = val;
						min_index = j;
					}
				}
				dp_column.set(current_index, i, min);
				min_recomb_index[i] = min_index;
			}
			if (!found_valid_transmission_vector) {
				throw std::runtime_error("Error: Mendelian conflict");
			}

			// if not last DP column, then update forward projection column and backtrace column
			if (next_input_column.get() == 0) {
				// update running optimal score index
				for (size_t i = 0; i < transmission_configurations; ++i) {
					if (dp_column.at(current_index, i) < optimal_score) {
						optimal_score = dp_column.at(current_index, i);
						optimal_score_index = iterator->get_index();
						optimal_transmission_value = i;
						previous_transmission_value = min_recomb_index[i];
					}
				}
			} else {
				unsigned int forward_index = iterator->get_forward_projection();
				unsigned int it_idx = iterator->get_index();
				for (unsigned int i = 0; i < transmission_configurations; ++i) {
					if(dp_column.at(current_index, i) < current_projection_column->at(forward_index,i)) {
						current_projection_column->set(forward_index, i, dp_column.at(current_index, i));
						index_backtrace_column->set(forward_index, i, it_idx);
						transmission_backtrace_column->set(forward_index,i, min_recomb_index[i]);
					}
				}
			}
		}

		// add newly computed backtrace_table column
		index_backtrace_table.push_back(index_backtrace_column);
		transmission_backtrace_table.push_back(transmission_backtrace_column);

		++n;
	} // end of main loop over columns
}

unsigned int PedigreeDPTable::get_optimal_score() {
	//if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
	return optimal_score;
}

unique_ptr<vector<index_and_inheritance_t> > PedigreeDPTable::get_index_path() {
	unique_ptr<vector<index_and_inheritance_t> > index_path = unique_ptr<vector<index_and_inheritance_t> >(new vector<index_and_inheritance_t>(indexers.size()));
	if (indexers.size() == 0) {
		return index_path;
	}
	index_and_inheritance_t v;
	unsigned int prev_inheritance_value = previous_transmission_value;
	v.index = optimal_score_index;
	v.inheritance_value = optimal_transmission_value;
	index_path->at(indexers.size()-1) = v;
	for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
		unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
		unsigned int backtrace_index = iterator->index_backward_projection(v.index);
		v.index = index_backtrace_table[i-1]->at(backtrace_index, v.inheritance_value);
		v.inheritance_value = prev_inheritance_value;
		prev_inheritance_value = transmission_backtrace_table[i-1]->at(backtrace_index, v.inheritance_value);
		index_path->at(i-1) = v;
	}

	return index_path;
}

void PedigreeDPTable::get_super_reads(std::vector<ReadSet*>* output_read_set, vector<unsigned int>* transmission_vector) {
	assert(output_read_set != nullptr);
	assert(output_read_set->size() == pedigree->size());
	assert(transmission_vector != nullptr);
	transmission_vector->clear();

	ColumnIterator column_iterator(*read_set);
	const vector<unsigned int>* positions = column_iterator.get_positions();

	std::vector<std::pair<Read*,Read*>> superreads;
	for (unsigned int i=0; i<pedigree->size(); i++) {
		 superreads.emplace_back(new Read("superread_0_"+std::to_string(i), -1, 0, -1),new Read("superread_1_"+std::to_string(i), -1, 0, -1));
	}

	if (index_backtrace_table.empty()) {
		assert(!column_iterator.has_next());
	} else {
		// run through the file again with the column_iterator
		unsigned int i = 0; // column index
		unique_ptr<vector<index_and_inheritance_t> > index_path = get_index_path();
		while (column_iterator.has_next()) {
			index_and_inheritance_t v = index_path->at(i);
			unique_ptr<vector<const Entry *> > column = column_iterator.get_next();
			PedigreeColumnCostComputer cost_computer(*column, i, read_sources, pedigree, *pedigree_partitions[v.inheritance_value]);
			cost_computer.set_partitioning(v.index);

			auto population_alleles = cost_computer.get_alleles();
			
			// TODO: compute proper weights based on likelihoods.
			for (unsigned int k=0; k<pedigree->size(); k++) {
				superreads[k].first->addVariant(positions->at(i), population_alleles[k].first, 0);
				superreads[k].second->addVariant(positions->at(i), population_alleles[k].second, 0);
			}
			transmission_vector->push_back(v.inheritance_value);
			++i; // next column
		}
	}
	for(unsigned int k=0;k<pedigree->size();k++) {
		assert(output_read_set->at(k) != nullptr);
		output_read_set->at(k)->add(superreads[k].first);
		output_read_set->at(k)->add(superreads[k].second);
	}
}

vector<bool>* PedigreeDPTable::get_optimal_partitioning() {
	unique_ptr<vector<index_and_inheritance_t> > index_path = get_index_path();
	vector<bool>* partitioning = new vector<bool>(read_set->size(),false);

	for(size_t i=0; i< index_path->size(); ++i) {
		unsigned int mask = 1; // mask to pass over the partitioning (i.e., index)
		for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
			unsigned int index = index_path->at(i).index;
			if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
				partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
			}
			mask = mask << 1;
		}
	}
	
	return partitioning;
}
