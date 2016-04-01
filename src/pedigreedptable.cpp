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

//#define DB // db

using namespace std;

// TODO read_marks is unnecessary because every read in the read_set already
// knows its id
PedigreeDPTable::PedigreeDPTable(ReadSet* read_set, vector<unsigned int> read_marks, vector<unsigned int> recombcost, Pedigree* pedigree)
	:
	read_set(read_set),
	read_marks(std::move(read_marks)),
	recombcost(std::move(recombcost)),
	triples(pedigree->triples),
	genotypes(pedigree->genotypes_map),
	id_of_individuals(pedigree->individuals),
	optimal_score(0u),
	optimal_score_index(0u),
	index_backtrace_table(),
	transmission_backtrace_table(),
	read_count(0u)
{
//   std::set<unsigned int> temp;
//   for(auto& triple: triples)
//       {
// 	for(auto& ind:triple)
// 	{
// 	  if(temp.insert(ind).second){
// 	    id_of_individuals.push_back(ind);
// 	  }
// 	}
//       }
  read_set->reassignReadIds();
  //for(unsigned int i=0;i<3;i++)
  //std::cout<<"compute table"<<genotypes[0][0][i]<<endl;
//    for(unsigned int i=0;i<3;i++)
//   std::cout<<"compute table"<<genotypes[0][1][i]<<endl;
//    for(unsigned int i=0;i<3;i++)
//   std::cout<<"compute table"<<genotypes[0][2][i]<<endl;
  compute_table();

}

PedigreeDPTable::~PedigreeDPTable() {
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

// db
#ifdef DB
// helper function to output the bit representation of an unsigned int
string bit_rep(unsigned int a, unsigned int len) {

  string s;
  for(int i=0; i< len; ++i) {
    s = ((a&1)?"1":"0") + s;
    a = a >> 1;
  }
  
  return s;
}

void output_vector(const vector<unsigned int> * v) {
  for(int j=v->size()-1; j>= 0; --j) {
    if(v->at(j) == -1) cout << "_ ";
    else cout << v->at(j) << " ";
  }
}

void output_vector_enum(const vector<unsigned int> * v, unsigned int len) {
  for(int j = v->size()-1; j >= 0; --j) {
    cout << j << " [" << bit_rep(j,len) << "] : ";
    if(v->at(j) == -1) cout << "_";
    else cout << v->at(j);
    cout << endl;
  }
}
#endif

void compute_final_cost(const num_of_recomb_uints_t& prev, const num_of_recomb_uints_t& current, unsigned int penalty, num_of_recomb_uints_t* min_costs, num_of_recomb_uints_t* min_cost_indices) {
  for (size_t i = 0; i < current.size(); ++i) {
    unsigned int min = numeric_limits<unsigned int>::max();
    size_t min_index = 0;
    for (size_t j = 0; j < prev.size(); ++j) {
      // Step 1: add up cost from current column and previous columns
      unsigned int val = current[i] + prev[j];
      // Step 2: add further cost incurred by recombination
      // change in bit 0 --> recombination in mother
      auto x = i ^ j; // count the number of bits set in x
      unsigned int count = 0; 
 
      for (;x; x >>= 1) {
	count += x & 1;
      }
      val += count * penalty;
      
      // check for new minimum
      if (val < min) {
        min = val;
        min_index = j;
      }
    }
    min_costs->at(i) = min;
    min_cost_indices->at(i) = min_index;
  }
}

void PedigreeDPTable::compute_table() {
  unsigned int num_recombs = std::pow(4, triples.size());
  ColumnIterator column_iterator(*read_set);
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

  // empty read-set, nothing to phase, so MEC score is 0
  if (!column_iterator.has_next()) {
    optimal_score = 0;
    optimal_score_index = 0;
    return;
  }  
  
  unsigned int n = 0;
  unique_ptr<vector<const Entry *> > current_column;
  unique_ptr<vector<const Entry *> > next_column;
  // get the next column ahead of time
  next_column = column_iterator.get_next();
  unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_column);
  ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0,*next_read_ids);
  indexers.push_back(next_indexer);
  unique_ptr<vector<num_of_recomb_uints_t> > previous_projection_column;
  unique_ptr<vector<num_of_recomb_uints_t> > current_projection_column;
  optimal_score = numeric_limits<unsigned int>::max();
  optimal_score_index = 0;
  optimal_transmission_value = 0;
  previous_transmission_value = 0;
  unsigned int nc = column_iterator.get_column_count();
//  if ((genotypesm.size() != nc) || (genotypesf.size() != nc) || (genotypesc.size() != nc)) {
//    throw std::runtime_error("Genotype vector length mismatch");
// }
#ifdef DB
  int i = 0;
#endif
  while(next_indexer != 0) {
   // std::cout<<"column"<<n<<endl;
    // move on projection column
    previous_projection_column = std::move(current_projection_column);
    // make former next column the current one
    current_column = std::move(next_column);
    unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
   
    ColumnIndexingScheme* current_indexer = next_indexer;
    // peek ahead and get the next column
    if (column_iterator.has_next()) {
      next_column = column_iterator.get_next();
      next_read_ids = extract_read_ids(*next_column);
      next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
      current_indexer->set_next_column(next_indexer);
      indexers.push_back(next_indexer);
    } else {
      assert(next_column.get() == 0);
      assert(next_read_ids.get() == 0);
      read_count = column_iterator.get_read_count();
      next_indexer = 0;
    }
    // reserve memory for the DP column
    num_of_recomb_uints_t null_array(num_recombs);
    vector<num_of_recomb_uints_t> dp_column(current_indexer->column_size(), null_array);
    vector<num_of_recomb_uints_t>* transmission_backtrace_column = nullptr;
    vector<num_of_recomb_uints_t>* index_backtrace_column = nullptr;
    // if not last column, reserve memory for forward projections column
    if (next_column.get() != 0) {
#ifdef DB
      cout << i << " : " << endl;
      ++i;
      cout << "allocate current projection / backtrace columns of size : " << current_indexer->forward_projection_size() << endl;
      cout << "forward projection width : " << current_indexer->get_forward_projection_width() << endl << endl;
#endif

      num_of_recomb_uints_t dummy_max_arr(num_recombs, numeric_limits<unsigned int>::max());
      current_projection_column = unique_ptr<vector<num_of_recomb_uints_t> >(
        new vector<num_of_recomb_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr)
      );
      // NOTE: forward projection size will always be even
      transmission_backtrace_column = new vector<num_of_recomb_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr);
      index_backtrace_column = new vector<num_of_recomb_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr);
    }

    // do the actual compution on current column
    vector<PedigreeColumnCostComputer> cost_computers;
    cost_computers.reserve(num_recombs);
    for(unsigned int i = 0; i < num_recombs; ++i) {
      cost_computers.emplace_back(*current_column, read_marks, i, triples);
    }

    unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();

    // db
#ifdef DB
    cout << "previous projection column (costs) :" << endl << endl;
    if(previous_projection_column.get()!=0) {
      output_vector_enum(previous_projection_column.get(),current_indexer->get_backward_projection_width());
    }
    cout << endl;

    cout << "row ids : ";
    output_vector(current_indexer->get_read_ids());

    cout << " .. column size : " << current_indexer->column_size() << endl;
    
    cout << "forward projection mask : ";
    if(next_column.get()!=0) {
      output_vector(current_indexer->get_forward_projection_mask());
      cout << " .. width : " << current_indexer->get_forward_projection_width();
    }
    cout << endl;

    cout << "------------------" << endl;
#endif
    
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

      // Fetch cost from previous column
      num_of_recomb_uints_t cost(num_recombs);
      if (previous_projection_column.get() != nullptr) {
        cost = previous_projection_column->at(iterator->get_backward_projection());
      }

#ifdef DB
      cout << iterator->get_backward_projection() << " [" << bit_rep(iterator->get_backward_projection(), current_indexer->get_backward_projection_width()) << "] -> " << cost;
#endif

      // TODO adapt to general case
      // Compute cost incurred by current cell of DP table
      num_of_recomb_uints_t current_cost;
      current_cost.reserve(num_recombs);
      std::vector<Pedigree::triple_entry_t> tempgeno;
    //  std::vector<std::vector<unsigned int>> tempgeno (genotypes.size(), std::vector<unsigned int>(3) );
      if (triples.empty())
      {
	tempgeno.emplace_back(Pedigree::triple_entry_t{(genotypes.begin()->second)[n],numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max()});
      }
      else{
      for(auto& triple: triples)
	{
	  auto t0 = genotypes.find(triple[0])->second;
	  auto t1 = genotypes.find(triple[1])->second;
	  auto t2 = genotypes.find(triple[2])->second;
	  tempgeno.emplace_back(Pedigree::triple_entry_t{t0[n], t1[n], t2[n]});
	}
      }
     // for (auto& genotriple: genotypes)
     // for(unsigned int h=0;h<genotypes.size();h++)
     // {
	 // tempgeno.emplace_back(std::vector<unsigned int>{genotriple[0][n],genotriple[1][n],genotriple[2][n]});
     // }
     // std::cout<<"in dptable"<<tempgeno[0][0]<< tempgeno[0][1]<<tempgeno[0][2]<<endl;
      for(auto& cost_computer : cost_computers) {
	  current_cost.push_back(cost_computer.get_cost(tempgeno));
      }

      

#ifdef DB
      cout << " + " << cost_computer.get_cost(genotypesm[n], genotypesf[n], genotypesc[n]) << " = " << cost << " -> " << iterator->get_index() << " [" << bit_rep(iterator->get_index(), current_indexer->get_read_ids()->size()) << "]";
      if(next_column.get()!=0) {
        cout << " -> " << iterator->get_forward_projection() << " [" << bit_rep(iterator->get_forward_projection(), current_indexer->get_forward_projection_width()) << "]";// fpw = " << current_indexer->get_forward_projection_width();
      }
      cout << endl;
#endif

      // Compute aggregate cost based on cost in previous and cost in current column
      num_of_recomb_uints_t final_col_cost(num_recombs);
      num_of_recomb_uints_t min_recomb_index(num_recombs);
      
      compute_final_cost(cost, current_cost, recombcost[n], &final_col_cost, &min_recomb_index);
      // ... and store it in current DP column
      dp_column[iterator->get_index()] = final_col_cost;
      // if not last DP column, then update forward projection column and backtrace column
      if (next_column.get() == 0) {
        // update running optimal score index
        for (size_t i = 0; i < final_col_cost.size(); ++i) {
          if (final_col_cost[i] < optimal_score) {
            optimal_score = final_col_cost[i];
            optimal_score_index = iterator->get_index();
            optimal_transmission_value = i;
            previous_transmission_value = min_recomb_index[i];
          }
        }
      } else {
        unsigned int forward_index = iterator->get_forward_projection();
        num_of_recomb_uints_t& current_proj_entry = current_projection_column->at(forward_index);
        num_of_recomb_uints_t& index_backtrace_column_entry = index_backtrace_column->at(forward_index);
        num_of_recomb_uints_t& transmission_backtrace_column_entry = transmission_backtrace_column->at(forward_index);
        unsigned int it_idx = iterator->get_index();
        for (unsigned int i = 0; i < current_proj_entry.size(); ++i) {
          if(final_col_cost[i] < current_proj_entry[i]) {
            current_proj_entry[i] = final_col_cost[i];
            index_backtrace_column_entry[i] = it_idx;
            transmission_backtrace_column_entry[i] = min_recomb_index[i];
          }
        }
      }
    }

#ifdef DB
    cout << endl;
#endif

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
  //cout<<"columns"<<indexers.size()<<endl;
  for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
    unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
    unsigned int backtrace_index = iterator->index_backward_projection(v.index);
    v.index = index_backtrace_table[i-1]->at(backtrace_index)[v.inheritance_value];
    v.inheritance_value = prev_inheritance_value;
    prev_inheritance_value = transmission_backtrace_table[i-1]->at(backtrace_index)[v.inheritance_value];
    index_path->at(i-1) = v;
  }

  //db
#ifdef DB
  cout << "index path : " << endl;
  output_vector(index_path.get());
  cout << endl;
#endif

  return index_path;
}

void PedigreeDPTable::get_super_reads(std::vector<ReadSet*>* output_read_set, vector<unsigned int>* transmission_vector) {
 // assert(std::all(output_read_set->begin(), output_read_set->end(), [](ReadSet* r){ return r != nullptr}));
 // assert(output_read_setf != 0u);
 // assert(output_read_setc != 0u);
  assert(transmission_vector != 0u);
  transmission_vector->clear();

  ColumnIterator column_iterator(*read_set);
  const vector<unsigned int>* positions = column_iterator.get_positions();
  
  std::vector<std::pair<Read*,Read*>> superreads;
  for(unsigned int i=0;i<id_of_individuals.size();i++)
  {
     superreads.emplace_back(new Read("superread_0_"+std::to_string(i), -1, 0),new Read("superread_1_"+std::to_string(i), -1, 0)) ;     
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
      PedigreeColumnCostComputer cost_computer(*column, read_marks, v.inheritance_value, triples);
      cost_computer.set_partitioning(v.index);
      std::vector<Pedigree::triple_entry_t> tempgeno;
    //  std::vector<std::vector<unsigned int>> tempgeno (genotypes.size(), std::vector<unsigned int>(3) );
      if (triples.empty())
      {
	tempgeno.emplace_back(Pedigree::triple_entry_t{(genotypes.begin()->second)[i],numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max()});
      }
      else{
      for(auto& triple: triples)
	{
	  auto t0 = genotypes.find(triple[0])->second;
	  auto t1 = genotypes.find(triple[1])->second;
	  auto t2 = genotypes.find(triple[2])->second;
	  tempgeno.emplace_back(Pedigree::triple_entry_t{t0[i], t1[i], t2[i]});
	}
      }

      auto population_alleles = cost_computer.get_alleles(tempgeno);
      
      // TODO: compute proper weights based on likelihoods.
       for(unsigned int k=0;k<id_of_individuals.size();k++)
       {
	 superreads[k].first->addVariant(positions->at(i), population_alleles[id_of_individuals[k]].first, 0);
	 superreads[k].second->addVariant(positions->at(i), population_alleles[id_of_individuals[k]].second, 0);
       }
      transmission_vector->push_back(v.inheritance_value);
      ++i; // next column
    }
  }
  for(unsigned int k=0;k<id_of_individuals.size();k++){
    (*output_read_set).at(k)->add(superreads[k].first);
    (*output_read_set).at(k)->add(superreads[k].second);
  }
}

vector<bool>* PedigreeDPTable::get_optimal_partitioning() {
  unique_ptr<vector<index_and_inheritance_t> > index_path = get_index_path();
  vector<bool>* partitioning = new vector<bool>(read_count,false);

  for(size_t i=0; i< index_path->size(); ++i) {

#ifdef DB
    cout << "index : " << index_path->at(i) << endl;
#endif

    unsigned int mask = 1; // mask to pass over the partitioning (i.e., index)
    for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {

#ifdef DB
      cout << indexers[i]->get_read_ids()->at(j) << " : ";
#endif

      unsigned int index = index_path->at(i).index;

#ifdef DB
      cout << index << " & " << mask << " = " << (index & mask);
#endif

      if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
        partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
#ifdef DB
        cout << " : true";
#endif
      }

#ifdef DB
      cout << endl;
#endif
      mask = mask << 1;
    }
  }
  
  return partitioning;
}
