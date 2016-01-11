#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>

#include "columncostcomputer.h"
#include "dptable.h"

//#define DB // db

using namespace std;

DPTable::DPTable(ReadSet* read_set, vector<unsigned int> read_marks, vector<unsigned int> recombcost, vector<unsigned int> genotypesm, vector<unsigned int> genotypesf, vector<unsigned int> genotypesc)
: read_set(read_set), read_marks(std::move(read_marks)), recombcost(std::move(recombcost)), genotypesm(std::move(genotypesm)), genotypesf(std::move(genotypesf)), genotypesc(std::move(genotypesc)),  indexers(), optimal_score(0u), optimal_score_index(0u),
 backtrace_table(), forrecomb(),read_count(0u)
{
  read_set->reassignReadIds();
  compute_table();

}

DPTable::~DPTable() {
  for(size_t i=0; i<indexers.size(); ++i) {
    delete indexers[i];
  }

  for(size_t i=0; i<backtrace_table.size(); ++i) {
    delete backtrace_table[i];
  }
  
  for(size_t i=0; i<forrecomb.size(); ++i) {
    delete forrecomb[i];
  }
}

unique_ptr<vector<unsigned int> > DPTable::extract_read_ids(const vector<const Entry *>& entries) {
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

void compute_final_cost(const four_uints_t& prev, const four_uints_t& current, unsigned int penalty, four_uints_t* min_costs, four_uints_t* min_cost_indices) {
  for (size_t i = 0; i < 4; ++i) {
    unsigned int min = numeric_limits<unsigned int>::max();
    size_t min_index = 0;
    for (size_t j = 0; j < 4; ++j) {
      // Step 1: add up cost from current column and previous columns
      unsigned int val = current[i] + prev[j];
      // Step 2: add further cost incurred by recombination
      // change in bit 0 --> recombination in mother
      if (i%2 != j%2) {
        val += penalty;
      }
      // change in bit 1 --> recombination in father
      if (i/2 != j/2) {
        val += penalty;
      }
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

void DPTable::compute_table() {
  ColumnIterator column_iterator(*read_set);
  if(!indexers.empty()) { // clear indexers, if present
    for(size_t i=0; i<indexers.size(); ++i) {
      delete indexers[i];
    }
    indexers.resize(0);
  }

  if(!backtrace_table.empty()) { // clear backtrace_table, if present
    for(size_t i=0; i<backtrace_table.size(); ++i) {
      delete backtrace_table[i];
    }
    backtrace_table.resize(0);
  }
  
  if(!forrecomb.empty()) { // clear backtrace_table, if present
    for(size_t i=0; i<forrecomb.size(); ++i) {
      delete forrecomb[i];
    }
    forrecomb.resize(0);
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
  unique_ptr<vector<four_uints_t> > previous_projection_column;
  unique_ptr<vector<four_uints_t> > current_projection_column;
  four_uints_t running_optimal_score;
  unsigned int running_optimal_score_index = 0; // optimal score and its index
  unsigned int nc = column_iterator.get_column_count();
  if ((genotypesm.size() != nc) || (genotypesf.size() != nc) || (genotypesc.size() != nc)) {
    throw std::runtime_error("Genotype vector length mismatch");
  }
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
    four_uints_t null_array = {{0, 0, 0, 0}};
    vector<four_uints_t> dp_column(current_indexer->column_size(), null_array);
    vector<four_uints_t>* forrecomb_col = nullptr;
    vector<four_uints_t>* backtrace_column = nullptr;
    // if not last column, reserve memory for forward projections column
    if (next_column.get() != 0) {
#ifdef DB
      cout << i << " : " << endl;
      ++i;
      cout << "allocate current projection / backtrace columns of size : " << current_indexer->forward_projection_size() << endl;
      cout << "forward projection width : " << current_indexer->get_forward_projection_width() << endl << endl;
#endif

      four_uints_t dummy_max_arr = {{numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max() }};
      current_projection_column = unique_ptr<vector<four_uints_t> >(
        new vector<four_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr)
      );
      // NOTE: forward projection size will always be even
      forrecomb_col = new vector<four_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr);
      backtrace_column = new vector<four_uints_t>(current_indexer->forward_projection_size(), dummy_max_arr);
    }

    // do the actual compution on current column
    ColumnCostComputer cost_computer_0(*current_column, read_marks, 0);
    ColumnCostComputer cost_computer_1(*current_column, read_marks, 1);
    ColumnCostComputer cost_computer_2(*current_column, read_marks, 2);
    ColumnCostComputer cost_computer_3(*current_column, read_marks, 3);
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
        cost_computer_0.update_partitioning(bit_changed);
        cost_computer_1.update_partitioning(bit_changed);
        cost_computer_2.update_partitioning(bit_changed);
        cost_computer_3.update_partitioning(bit_changed);
      } else {
        cost_computer_0.set_partitioning(iterator->get_partition());
        cost_computer_1.set_partitioning(iterator->get_partition());
        cost_computer_2.set_partitioning(iterator->get_partition());
        cost_computer_3.set_partitioning(iterator->get_partition());
        if(next_column.get() == 0) { // only if we're at the last column
          running_optimal_score_index = iterator->get_index(); // default to first
        }
      }

      four_uints_t cost = {{0, 0, 0, 0}};
      if (previous_projection_column.get() != nullptr) {
        cost = previous_projection_column->at(iterator->get_backward_projection());
      }

#ifdef DB
      cout << iterator->get_backward_projection() << " [" << bit_rep(iterator->get_backward_projection(), current_indexer->get_backward_projection_width()) << "] -> " << cost;
#endif

      four_uints_t current_cost = {{ 
        cost_computer_0.get_cost(genotypesm[n], genotypesf[n], genotypesc[n]), 
        cost_computer_1.get_cost(genotypesm[n], genotypesf[n], genotypesc[n]),
        cost_computer_2.get_cost(genotypesm[n], genotypesf[n], genotypesc[n]),
        cost_computer_3.get_cost(genotypesm[n], genotypesf[n], genotypesc[n])
      }};
      

#ifdef DB
      cout << " + " << cost_computer.get_cost(genotypesm[n], genotypesf[n], genotypesc[n]) << " = " << cost << " -> " << iterator->get_index() << " [" << bit_rep(iterator->get_index(), current_indexer->get_read_ids()->size()) << "]";
      if(next_column.get()!=0) {
        cout << " -> " << iterator->get_forward_projection() << " [" << bit_rep(iterator->get_forward_projection(), current_indexer->get_forward_projection_width()) << "]";// fpw = " << current_indexer->get_forward_projection_width();
      }
      cout << endl;
#endif
       
      four_uints_t final_col_cost = {0,0,0,0};
      four_uints_t min_recomb_index = {0,0,0,0};
      compute_final_cost(cost, current_cost, recombcost[n], &final_col_cost, &min_recomb_index);
      dp_column[iterator->get_index()] = final_col_cost;
      // if not last DP column, then update forward projection column and backtrace column
      if (next_column.get() == 0) {
        // update running optimal score index
        auto& current_optimal_cost = dp_column[running_optimal_score_index];
        if(*min_element(begin(final_col_cost), end(final_col_cost)) < *min_element(begin(current_optimal_cost), end(current_optimal_cost))){
          running_optimal_score_index = iterator->get_index();
        }
      } else {
        unsigned int forward_index = iterator->get_forward_projection();
        auto& current_proj_entry = (*current_projection_column)[forward_index];
        auto& backtrace_column_entry = (*backtrace_column)[forward_index];
        auto& recombminindex_column_entry = (*forrecomb_col)[forward_index];
        auto it_idx = iterator->get_index();
        for (unsigned int i = 0; i < 4; ++i) {
          if(final_col_cost[i] < current_proj_entry[i]) {
            current_proj_entry[i] = final_col_cost[i];
            backtrace_column_entry[i] = it_idx;
            recombminindex_column_entry[i] = min_recomb_index[i];
          }
        }
      }
    }

#ifdef DB
    cout << endl;
#endif

    if(next_column.get() == 0) { // record optimal score
      running_optimal_score = dp_column[running_optimal_score_index];
    }

    // add newly computed backtrace_table column
    backtrace_table.push_back(backtrace_column);
    forrecomb.push_back(forrecomb_col);

    ++n;

  } // end of main loop over columns

  // store optimal score for table at end of computation
  auto min = numeric_limits<unsigned int>::max();
  for(unsigned int i = 0; i < 4; ++i) {
    if(running_optimal_score[i] < min) {
      min = running_optimal_score[i];
      optimal_score_array_index = i;
    }
  }
  optimal_score = min;
  optimal_score_index = running_optimal_score_index;
}

unsigned int DPTable::get_optimal_score() {
  //if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
  return optimal_score;
}

unique_ptr<vector<index_and_inheritance_t> > DPTable::get_index_path() {
  unique_ptr<vector<index_and_inheritance_t> > index_path = unique_ptr<vector<index_and_inheritance_t> >(new vector<index_and_inheritance_t>(indexers.size()));
  if (indexers.size() == 0) {
    return index_path;
  }
  index_and_inheritance_t v;
  v.index = optimal_score_index;
  v.inheritance_value = optimal_score_array_index;
  index_path->at(indexers.size()-1) = v;
  //cout<<"columns"<<indexers.size()<<endl;
  for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
    unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
    unsigned int backtrace_index = iterator->index_backward_projection(v.index);
    v.inheritance_value = forrecomb[i-1]->at(backtrace_index)[v.inheritance_value];
    v.index = backtrace_table[i-1]->at(backtrace_index)[v.inheritance_value];
     
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

void DPTable::get_super_reads(ReadSet* output_read_setm, ReadSet* output_read_setf, ReadSet* output_read_setc, vector<unsigned int>* transmission_vector) {
  assert(output_read_setm != 0u);
  assert(output_read_setf != 0u);
  assert(output_read_setc != 0u);
  assert(transmission_vector != 0u);
  transmission_vector->clear();

  ColumnIterator column_iterator(*read_set);
  const vector<unsigned int>* positions = column_iterator.get_positions();
  
  Read* r0m = new Read("superread_mother0", -1, 0);
  Read* r1m = new Read("superread_mother1", -1, 0);
  Read* r0f = new Read("superread_father0", -1, 0);
  Read* r1f = new Read("superread_father1", -1, 0);
  Read* r0c = new Read("superread_child0", -1, 0);
  Read* r1c = new Read("superread_child1", -1, 0);
  
  if (backtrace_table.empty()) {
    assert(!column_iterator.has_next());
  } else {
    // run through the file again with the column_iterator
    unsigned int i = 0; // column index
    unique_ptr<vector<index_and_inheritance_t> > index_path = get_index_path();
    while (column_iterator.has_next()) {
      index_and_inheritance_t v = index_path->at(i);
      unique_ptr<vector<const Entry *> > column = column_iterator.get_next();
      ColumnCostComputer cost_computer(*column, read_marks, v.inheritance_value);
      cost_computer.set_partitioning(v.index);
      ColumnCostComputer::trio_alleles_t trio_alleles = cost_computer.get_alleles(genotypesm[i], genotypesf[i], genotypesc[i]);
      // TODO: compute proper weights based on likelihoods.
      r0m->addVariant(positions->at(i), trio_alleles.mother.first, 0);
      r1m->addVariant(positions->at(i), trio_alleles.mother.second, 0);
      r0f->addVariant(positions->at(i), trio_alleles.father.first, 0);
      r1f->addVariant(positions->at(i), trio_alleles.father.second, 0);
      r0c->addVariant(positions->at(i), trio_alleles.child.first, 0);
      r1c->addVariant(positions->at(i), trio_alleles.child.second, 0);
      transmission_vector->push_back(v.inheritance_value);
      ++i; // next column
    }
  }

  output_read_setm->add(r0m);
  output_read_setm->add(r1m);
  output_read_setf->add(r0f);
  output_read_setf->add(r1f);
  output_read_setc->add(r0c);
  output_read_setc->add(r1c);
}

vector<bool>* DPTable::get_optimal_partitioning() {
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
