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

DPTable::DPTable(ReadSet* read_set, vector<unsigned int> read_marks, vector<unsigned int> recombcost,bool all_heterozygous)
: read_set(read_set), read_marks(std::move(read_marks)), recombcost(std::move(recombcost)),indexers(), optimal_score(0u), optimal_score_index(0u),
 backtrace_table(), forrecomb(),read_count(0u), all_heterozygous(all_heterozygous)
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

namespace {
     template<typename T, std::size_t s>
     array<T, s> compute_final_cost(array<T,s>& prev, array<T,s>& current, unsigned int penalty) {
       
       array<T, s> result;
       for(typename array<T,s>::size_type i = 0; i < s; ++i) {
         unsigned int min = numeric_limits<unsigned int>::max();
         for(typename array<T,s>::size_type j = 0; j < s; ++j) {
             auto val = current[i] + prev[j] + (i != j ? penalty : 0u);
             if(val < min) {
               min = val;
             }
         }
         result[i] = min;
       }
      
       return result;
     }
     
     template<typename T, std::size_t s>
     array<T, s> compute_min_index(array<T,s>& prev, array<T,s>& current, unsigned int penalty) {
       
       array<T, s> result;
       
       for(typename array<T,s>::size_type i = 0; i < s; ++i) {
         unsigned int min = numeric_limits<unsigned int>::max();
         for(typename array<T,s>::size_type j = 0; j < s; ++j) {
             auto val = current[i] + prev[j] + (i != j ? penalty : 0u);
             if(val < min) {
               min = val;
               result[i] = j;
             }
         }
         
       }
      
       return result;
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
  unique_ptr<vector<array<unsigned int, 4> > > previous_projection_column;
  unique_ptr<vector<array<unsigned int, 4> > > current_projection_column;
  array<unsigned int, 4> running_optimal_score;
  unsigned int running_optimal_score_index; // optimal score and its index
  unsigned int temp;
  double pi = 0.05; // percentage of columns processed
  double pc = pi;
  unsigned int nc = column_iterator.get_column_count();
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
    array<unsigned int, 4> null_array = {{0, 0, 0, 0}};
    vector<array<unsigned int, 4>> dp_column(current_indexer->column_size(), null_array);
    vector<array<unsigned int, 4>>* forrecomb_col = nullptr;
    vector<array<unsigned int, 4>>* backtrace_column = nullptr;
    // if not last column, reserve memory for forward projections column
    if (next_column.get() != 0) {
#ifdef DB
      cout << i << " : " << endl;
      ++i;
      cout << "allocate current projection / backtrace columns of size : " << current_indexer->forward_projection_size() << endl;
      cout << "forward projection width : " << current_indexer->get_forward_projection_width() << endl << endl;
#endif

    array<unsigned int, 4> dummy_max_arr = {{numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max(),
                                             numeric_limits<unsigned int>::max(), numeric_limits<unsigned int>::max()
                                   }};
      current_projection_column = unique_ptr<vector<array<unsigned int, 4>> >(
        new vector<array<unsigned int, 4>>(current_indexer->forward_projection_size(), dummy_max_arr));
      // NOTE: forward projection size will always be even
      
      forrecomb_col = new vector<array<unsigned int, 4>>(current_indexer->forward_projection_size(), dummy_max_arr);

      backtrace_column = new vector<array<unsigned int, 4>>(current_indexer->forward_projection_size(), dummy_max_arr);

    }

    // do the actual compution on current column
    ColumnCostComputer cost_computer_0(*current_column, read_marks, 0, all_heterozygous);
    ColumnCostComputer cost_computer_1(*current_column, read_marks, 1, all_heterozygous);
    ColumnCostComputer cost_computer_2(*current_column, read_marks, 2, all_heterozygous);
    ColumnCostComputer cost_computer_3(*current_column, read_marks, 3, all_heterozygous);
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

      array<unsigned int, 4> cost = {{0, 0, 0, 0}};
      auto* prev_col_ptr = previous_projection_column.get();
      if (prev_col_ptr != nullptr) {
        auto& prev_col_entry = (*prev_col_ptr)[iterator->get_backward_projection()];
        for(unsigned int i = 0; i < 4; ++i) {
          cost[i] += prev_col_entry[i];
        }
      }

      // db
#ifdef DB
      cout << iterator->get_backward_projection() << " [" << bit_rep(iterator->get_backward_projection(), current_indexer->get_backward_projection_width()) << "] -> " << cost;
#endif
    //  std::cout<<"prevcost"<<cost[0]<<cost[1]<<cost[2]<<cost[3];
     // std:: cout <<"bla bla"<<cost_computer_3.get_cost()<< cost_computer_1.get_cost()<<
        cost_computer_2.get_cost()<< cost_computer_3.get_cost();

      array<unsigned int, 4> current_cost = {{ 
        cost_computer_0.get_cost(), cost_computer_1.get_cost(),
        cost_computer_2.get_cost(), cost_computer_3.get_cost()
      }};
      

      // db
#ifdef DB
      cout << " + " << cost_computer.get_cost() << " = " << cost << " -> " << iterator->get_index() << " [" << bit_rep(iterator->get_index(), current_indexer->get_read_ids()->size()) << "]";
      if(next_column.get()!=0) {
        cout << " -> " << iterator->get_forward_projection() << " [" << bit_rep(iterator->get_forward_projection(), current_indexer->get_forward_projection_width()) << "]";// fpw = " << current_indexer->get_forward_projection_width();
      }
      cout << endl;
#endif
       
      auto final_col_cost = compute_final_cost(cost, current_cost, recombcost[n]);
      auto min_recomb_index=compute_min_index(cost, current_cost, recombcost[n]);
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
          }
        }
        for (unsigned int i = 0; i < 4; ++i) {
         recombminindex_column_entry[i]=min_recomb_index[i];
          
        }
      }
    }

    // db
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

    // db
    /*
    for(size_t j=0;j<current_column->size(); ++j) {
      cout << dp_column->at(j) << endl;
    }
    cout << endl;
    */

    // completion percentage output
    /*
    if((double(n)/double(nc)) > pc) {
      cout << int(pc*100.0) << " %" << endl;
      pc += pi;
    }
    */
  }

  // store optimal score for table at end of computation
  auto min = numeric_limits<unsigned int>::max();
  for(unsigned int i = 0; i < 4; ++i) {
    if(running_optimal_score[i] < min) {
      min = running_optimal_score[i];
      optimal_score_array_index = i;
    }
  }
  optimal_score_index = running_optimal_score_index;
}

unsigned int DPTable::get_optimal_score() {
  //if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
  return optimal_score;
}

unique_ptr<vector<unsigned int> > DPTable::get_index_path() {

  unique_ptr<vector<unsigned int> > index_path = unique_ptr<vector<unsigned int> >(new vector<unsigned int>(indexers.size()));
  if (indexers.size() == 0) {
    return index_path;
  }
  unsigned int index = optimal_score_index;
  unsigned int temp;
  index_path->at(indexers.size()-1) = index;
  //cout<<"columns"<<indexers.size()<<endl;
  for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
    unique_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
    unsigned int backtrace_index = iterator->index_backward_projection(index);
    if (i==indexers.size()-1){
    index = backtrace_table[i-1]->at(backtrace_index)[optimal_score_array_index];
    temp= forrecomb[i-1]->at(backtrace_index)[optimal_score_array_index];
    }
    else{
      index = backtrace_table[i-1]->at(backtrace_index)[temp];
    temp= forrecomb[i-1]->at(backtrace_index)[temp]; 
    }
    
    index_path->at(i-1) = index;
  }

  //db
#ifdef DB
  cout << "index path : " << endl;
  output_vector(index_path.get());
  cout << endl;
#endif

  return index_path;
}

void DPTable::get_super_readsm(ReadSet* output_read_set) {
  assert(output_read_set != 0u);

  ColumnIterator column_iterator(*read_set);
  const vector<unsigned int>* positions = column_iterator.get_positions();
// std::cout<<positions<<endl;
 // for(int i=0;i<positions->size();i++) std::cout<<positions->at(i)<<endl;
 //std::cout<<positions<<endl;

  Read* r0m = new Read("superread0", -1, 0);
  Read* r1m = new Read("superread1", -1, 0);
  
  if (backtrace_table.empty()) {
    assert(!column_iterator.has_next());
  } else {
    // run through the file again with the column_iterator
    unsigned int i = 0; // column index
    unique_ptr<vector<unsigned int> > index_path = get_index_path();
int count=0;
    while (column_iterator.has_next()) {
      unsigned int index = index_path->at(i);
      unique_ptr<vector<const Entry *> > column = column_iterator.get_next();
    const std::vector<const Entry*>& columnx=  *column;
    const std::vector<unsigned int>& read_marksx=read_marks;
     // unsigned int index_m=0;
      
       bool flag= false;
    
     //  int j=0;
      for (vector<const Entry*>::const_iterator it = columnx.begin(); it != columnx.end(); ++it) {
       // unsigned int temp=index%2;
        auto& entry = **it;
        if((read_marksx[entry.get_read_id()]==1u) && ((*it)->get_allele_type())!= Entry::BLANK) {
          flag=true;
          break;
         // index_m+=pow(2,j)*temp;
         // j++;
        }      
        //index/=2;
     }      
      
     if(flag){
    count++;
        ColumnCostComputer cost_computer(*column, read_marks, 0, all_heterozygous);
        cost_computer.set_partitioning_m(index);
        r0m->addVariant(positions->at(i), cost_computer.get_allele(0), cost_computer.get_weight(0));
        r1m->addVariant(positions->at(i), cost_computer.get_allele(1), cost_computer.get_weight(1));
     }
      ++i; // next column
    }
  }
 
  output_read_set->add(r0m);
  output_read_set->add(r1m);
}


//---------------------------

void DPTable::get_super_readsf(ReadSet* output_read_set) {
  assert(output_read_set != 0);

  ColumnIterator column_iterator(*read_set);
  const vector<unsigned int>* positions = column_iterator.get_positions();

  Read* r0f = new Read("superread0", -1, 0);
  Read* r1f = new Read("superread1", -1, 0);
  
  if (backtrace_table.empty()) {
    assert(!column_iterator.has_next());
  } else {
    // run through the file again with the column_iterator
    unsigned int i = 0; // column index
    unique_ptr<vector<unsigned int> > index_path = get_index_path();
int count=0;
    while (column_iterator.has_next()) {
      unsigned int index = index_path->at(i);
      unique_ptr<vector<const Entry *> > column = column_iterator.get_next();
    const std::vector<const Entry*>& columnx=  *column;
    const std::vector<unsigned int>& read_marksx=read_marks;
     // unsigned int index_m=0;
      
       bool flag= false;
    
     //  int j=0;
       for (vector<const Entry*>::const_iterator it = columnx.begin(); it != columnx.end(); ++it) {
       // unsigned int temp=index%2;
        auto& entry = **it;
        if((read_marksx[entry.get_read_id()]==2u) && ((*it)->get_allele_type())!= Entry::BLANK) {
          flag=true;
          break;
         // index_m+=pow(2,j)*temp;
         // j++;
        }      
        //index/=2;
     }      
      
     if(flag){
    count++;
        ColumnCostComputer cost_computer(*column, read_marks, 0, all_heterozygous);
        cost_computer.set_partitioning_f(index);

        r0f->addVariant(positions->at(i), cost_computer.get_allele(0), cost_computer.get_weight(0));
        r1f->addVariant(positions->at(i), cost_computer.get_allele(1), cost_computer.get_weight(1));
     }
      ++i; // next column
    }
  }
 
  output_read_set->add(r0f);
  output_read_set->add(r1f);
}


//-----------------------

void DPTable::get_super_readsc(ReadSet* output_read_set) {
  assert(output_read_set != 0);

  ColumnIterator column_iterator(*read_set);
  const vector<unsigned int>* positions = column_iterator.get_positions();

  Read* r0c = new Read("superread0", -1, 0);
  Read* r1c = new Read("superread1", -1, 0);
  
  if (backtrace_table.empty()) {
    assert(!column_iterator.has_next());
  } else {
    // run through the file again with the column_iterator
    unsigned int i = 0; // column index
int count=0;
    unique_ptr<vector<unsigned int> > index_path = get_index_path();
    while (column_iterator.has_next()) {
      unsigned int index = index_path->at(i);
      unique_ptr<vector<const Entry *> > column = column_iterator.get_next();
    const std::vector<const Entry*>& columnx=  *column;
    const std::vector<unsigned int>& read_marksx=read_marks;
     // unsigned int index_m=0;
      
       bool flag= false;
    
     //  int j=0;
      for (vector<const Entry*>::const_iterator it = columnx.begin(); it != columnx.end(); ++it) {
       // unsigned int temp=index%2;
        auto& entry = **it;
        if((read_marksx[entry.get_read_id()]==0u) && ((*it)->get_allele_type())!= Entry::BLANK) {
          flag=true;
          break;
         // index_m+=pow(2,j)*temp;
         // j++;
        }      
        //index/=2;
     }      
      
     if(flag){
    count++;
        ColumnCostComputer cost_computer(*column, read_marks, 0, all_heterozygous);
        cost_computer.set_partitioning_c(index);

        r0c->addVariant(positions->at(i), cost_computer.get_allele(0), cost_computer.get_weight(0));
        r1c->addVariant(positions->at(i), cost_computer.get_allele(1), cost_computer.get_weight(1));
     }
      ++i; // next column
    }
  }

  output_read_set->add(r0c);
  output_read_set->add(r1c);
}

//---------

vector<bool>* DPTable::get_optimal_partitioning() {

  unique_ptr<vector<unsigned int> > index_path = get_index_path();
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

      unsigned int index = index_path->at(i);

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
