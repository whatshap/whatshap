#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>

#include "columncostcomputer.h"
#include "dptable.h"

//#define HALF_TABLE // switch on to use the trick that uses half the memory (when this is complete, it will be this by default)
//#define DB // db

using namespace std;

DPTable::DPTable(bool all_heterozygous) {
  this->all_heterozygous = all_heterozygous;
  this->read_count = 0;
}

auto_ptr<vector<unsigned int> > DPTable::extract_read_ids(const vector<const Entry *>& entries) {
  auto_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
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

void DPTable::compute_table(ColumnIterator * column_iterator) {
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
  
  if (!column_iterator->has_next()) return;
  
  unsigned int n = 0;
  auto_ptr<vector<const Entry *> > current_column(0);
  auto_ptr<vector<const Entry *> > next_column(0);
  // get the next column ahead of time
  next_column = column_iterator->get_next();
  auto_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_column);
  ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0,*next_read_ids);
  indexers.push_back(next_indexer);
  auto_ptr<vector<unsigned int> > previous_projection_column(0);
  auto_ptr<vector<unsigned int> > current_projection_column(0);
  unsigned int running_optimal_score;
  unsigned int running_optimal_score_index; // optimal score and its index
  double pi = 0.05; // percentage of columns processed
  double pc = pi;
  unsigned int nc = column_iterator->get_column_count();
  while(next_indexer != 0) {
    // move on projection column
    previous_projection_column = current_projection_column;
    // make former next column the current one
    current_column = next_column;
    auto_ptr<vector<unsigned int> > current_read_ids = next_read_ids;
    ColumnIndexingScheme* current_indexer = next_indexer;
    // peek ahead and get the next column
    if (column_iterator->has_next()) {
      next_column = column_iterator->get_next();
      next_read_ids = extract_read_ids(*next_column);
      next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
      current_indexer->set_next_column(next_indexer);
      indexers.push_back(next_indexer);
    } else {
      assert(next_column.get() == 0);
      assert(next_read_ids.get() == 0);
      read_count = column_iterator->get_read_count();
      next_indexer = 0;
    }
    // reserve memory for the DP column
    vector<unsigned int> dp_column(current_indexer->column_size(),0);
    vector<unsigned int>* backtrace_column = 0;
    // if not last column, reserve memory for forward projections column
    if (next_column.get() != 0) {
      current_projection_column = auto_ptr<vector<unsigned int> >(new vector<unsigned int>(current_indexer->forward_projection_size(), numeric_limits<unsigned int>::max()));
      // NOTE: forward projection size will always be even
#ifdef HALF_TABLE
      backtrace_column = new vector<unsigned int>((current_indexer->forward_projection_size()>>1), numeric_limits<unsigned int>::max());
#else
      backtrace_column = new vector<unsigned int>(current_indexer->forward_projection_size(), numeric_limits<unsigned int>::max());
#endif
    }

    // do the actual compution on current column
    ColumnCostComputer cost_computer(*current_column, all_heterozygous);
    auto_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();

    // db
#ifdef DB
    cout << "previous projection column (costs) :" << endl << endl;
    if(previous_projection_column.get()!=0) {
      output_vector_enum(previous_projection_column.get(),current_indexer->get_backward_projection_width());
    }
    cout << endl;

    cout << "row ids : ";
    output_vector(current_indexer->get_read_ids());
    cout << "(pivot)" << endl;

    cout << "column size : " << current_indexer->column_size() << endl;
    
    cout << "forward projection mask : ";
    if(next_column.get()!=0) {
      output_vector(current_indexer->get_forward_projection_mask());
    }
    cout << endl;

    cout << "------------------" << endl;
#endif
    
    while (iterator->has_next()) {
      int bit_changed = -1;
      iterator->advance(&bit_changed);
      if (bit_changed >= 0) {
        cost_computer.update_partitioning(bit_changed);
      } else {
        cost_computer.set_partitioning(iterator->get_partition());
        if(next_column.get() == 0) { // only if we're at the last column
          running_optimal_score_index = iterator->get_index(); // default to first
        }
      }

      unsigned int cost = 0;
      if (previous_projection_column.get() != 0) {
        cost += previous_projection_column->at(iterator->get_backward_projection());
      }

      // db
#ifdef DB
      cout << iterator->get_backward_projection() << " [" << bit_rep(iterator->get_backward_projection(), current_indexer->get_backward_projection_width()) << "] -> " << cost;
#endif

      cost += cost_computer.get_cost();

      // db
#ifdef DB
      cout << " + " << cost_computer.get_cost() << " = " << cost << " -> " << iterator->get_index() << " [" << bit_rep(iterator->get_index(), current_indexer->get_read_ids()->size()) << "]";
      if(next_column.get()!=0) {
        cout << " -> " << iterator->get_forward_projection() << " [" << bit_rep(iterator->get_forward_projection(), current_indexer->get_forward_projection_width()) << "]";
      }
      cout << endl;
#endif

      dp_column[iterator->get_index()] = cost;
      // if not last DP column, then update forward projection column and backtrace column
      if (next_column.get() == 0) {
        // update running optimal score index
        if (cost < dp_column[running_optimal_score_index]) {
          running_optimal_score_index = iterator->get_index();
        }
      } else {
        unsigned int forward_index = iterator->get_forward_projection();
        if (current_projection_column->at(forward_index) > cost) {
          current_projection_column->at(forward_index) = cost;
#ifdef HALF_TABLE
          forward_index = min(forward_index, iterator->get_forward_dual_projection());
          // ensure that we pass the smaller of the two indexes,
          // because backtrace_column is half the size of
          // current_projection_column (which is symmetric, hence this
          // dual computation)
#endif
          backtrace_column->at(forward_index) = iterator->get_index();
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
  optimal_score = running_optimal_score;
  optimal_score_index = running_optimal_score_index;
}

unsigned int DPTable::get_optimal_score() {
  if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
  return optimal_score;
}

auto_ptr<vector<unsigned int> > DPTable::get_index_path() {

  auto_ptr<vector<unsigned int> > index_path = auto_ptr<vector<unsigned int> >(new vector<unsigned int>(indexers.size()));
  unsigned int index = optimal_score_index;
  index_path->at(indexers.size()-1) = index;
  for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
    if(i>0) {
      auto_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
      unsigned int backtrace_index = iterator->index_backward_projection(index);
      index = backtrace_table[i-1]->at(backtrace_index);
      index_path->at(i-1) = index;
    }
  }

  return index_path;
}

auto_ptr<DPTable::haplotype_t> DPTable::get_haplotype(ColumnIterator * column_iterator) {

  auto_ptr<DPTable::read_t> s = get_super_reads(column_iterator);
  unsigned int c_len = s->at(0).size();
  vector<Entry::allele_t> h0(c_len);
  vector<Entry::allele_t> h1(c_len);

  for(unsigned int i=0; i< c_len; ++i) { // loop over the positions in super-read

    h0[i] = s->at(0)[i]->get_allele_type();
    h1[i] = s->at(1)[i]->get_allele_type();
  }

  auto_ptr<DPTable::haplotype_t> h = auto_ptr<DPTable::haplotype_t>(new DPTable::haplotype_t);
  h->push_back(h0);
  h->push_back(h1);

  return h;
}

// what used to be 'get_haplotype', but this is a more accurate name
auto_ptr<DPTable::haplotype_t> DPTable::get_raw_haplotype(ColumnIterator * column_iterator) {

  auto_ptr<DPTable::read_t> s = get_super_reads(column_iterator);
  unsigned int c_len = s->at(0).size();
  unsigned int p_len = s->at(0)[c_len-1]->get_read_id();
  vector<Entry::allele_t> h0(p_len);
  vector<Entry::allele_t> h1(p_len);

  unsigned int offset = s->at(0)[0]->get_read_id(); // first position
  unsigned int j = 0;
  for(unsigned int i=0; i< c_len; ++i) { // loop over the positions

    while((j+offset) < s->at(0)[i]->get_read_id()) { // fill gaps with '-'s
      h0[j] = Entry::BLANK;
      h1[j] = Entry::BLANK;
      ++j;
    }

    h0[j] = s->at(0)[i]->get_allele_type();
    h1[j] = s->at(1)[i]->get_allele_type();
    ++j; // next position
  }

  auto_ptr<DPTable::haplotype_t> h = auto_ptr<DPTable::haplotype_t>(new DPTable::haplotype_t);
  h->push_back(h0);
  h->push_back(h1);
  
  return h;
}

auto_ptr<DPTable::read_t> DPTable::get_super_reads(ColumnIterator * column_iterator) {

  auto_ptr<DPTable::read_t> r = auto_ptr<DPTable::read_t>(new DPTable::read_t);

  if(!column_iterator->has_next()) return r;

  const vector<unsigned int> * positions = column_iterator->get_positions();

  unsigned int c_len = positions->size();
  vector<Entry*> r0(c_len);
  vector<Entry*> r1(c_len);

  // run through the file again with the column_iterator
  unsigned int i = 0; // column index
  auto_ptr<vector<unsigned int> > index_path = get_index_path();
  while(column_iterator->has_next()) {
    unsigned int index = index_path->at(i);
    auto_ptr<vector<const Entry *> > column = column_iterator->get_next();
    ColumnCostComputer cost_computer(*column, all_heterozygous);
    cost_computer.set_partitioning(index);

    // yes, I abuse slightly Entry in using read_id for a position,
    // but it's just the perfect object to represent a snp position of
    // a .wif file.  Perhaps we can generalize Entry by naming the
    // read_id field 'info' or 'tag' or something
    r0[i] = new Entry(positions->at(i), cost_computer.get_allele(0), cost_computer.get_weight(0));
    r1[i] = new Entry(positions->at(i), cost_computer.get_allele(1), cost_computer.get_weight(1));
    ++i; // next column
  }

  r->push_back(r0);
  r->push_back(r1);

  return r;
}

auto_ptr<vector<bool> > DPTable::get_optimal_partitioning() {

  auto_ptr<vector<unsigned int> > index_path = get_index_path();

  // db
  /*
  for(size_t i=0; i< index_path.size(); ++i) {
    cout << "index : " << index_path[i] << " " << endl;
    for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
      cout << indexers[i]->get_read_ids()->at(j) << endl;
    }
  }
  */

  auto_ptr<vector<bool> > partitioning = auto_ptr<vector<bool> >(new vector<bool>(read_count,false));
  for(size_t i=0; i< index_path->size(); ++i) {
    unsigned int mask = 1; // mask to pass over the partitioning (i.e., index)
    for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
      unsigned int index = index_path->at(i);
      if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
        partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
      }
      mask = mask << 1;
    }
  }
  
  return partitioning;
}

// HALF_TABLE stuff ...

auto_ptr<vector<unsigned int> > DPTable::get_index_path_HALF_TABLE() {

  auto_ptr<vector<unsigned int> > index_path = auto_ptr<vector<unsigned int> >(new vector<unsigned int>(indexers.size()));
  unsigned int index = optimal_score_index;
  index_path->at(indexers.size()-1) = index;
  for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
    if(i>0) {
      auto_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
      unsigned int backtrace_index = min(iterator->index_backward_projection(index), iterator->dual_index_backward_projection(index));
      index = backtrace_table[i-1]->at(backtrace_index);
      index_path->at(i-1) = index; // index is re-used again!! keep this in mind!!
    }
  }

  return index_path;
}

auto_ptr<DPTable::haplotype_t> DPTable::get_haplotype_HALF_TABLE(ColumnIterator & column_iterator) {

  auto_ptr<DPTable::haplotype_t> h = auto_ptr<DPTable::haplotype_t>(new DPTable::haplotype_t);

  if(!column_iterator.has_next()) return h;

  const vector<unsigned int> * positions = column_iterator.get_positions();

  unsigned int p_len = positions->size();
  unsigned int c_len = positions->at(p_len-1);
  vector<Entry::allele_t> h0(c_len);
  vector<Entry::allele_t> h1(c_len);

  // run through the file again with the column_iterator
  unsigned int j = positions->at(0)-1; // position index
  auto_ptr<vector<unsigned int> > index_path = get_index_path();

  auto_ptr<ColumnIndexingIterator> previous_iterator(0);
  auto_ptr<vector<const Entry *> > previous_column(0);
  unsigned int previous_index = -1;

  auto_ptr<ColumnIndexingIterator> iterator = indexers[0]->get_iterator();
  auto_ptr<vector<const Entry *> > column = column_iterator.get_next();
  unsigned int index = index_path->at(0);
  
  // start off the haplotype
  ColumnCostComputer first_cost_computer(*(column.get()), all_heterozygous);
  first_cost_computer.set_partitioning(index);
  h0[j] = first_cost_computer.get_allele(0);
  h1[j] = first_cost_computer.get_allele(1);
  ++j; // next position
  
  for(unsigned int i= 1; i < p_len; ++i) { // for each column i
    assert(column_iterator.has_next());

    while(j < positions->at(i)-1) { // fill gaps with '-'s
      h0[j] = Entry::BLANK;
      h1[j] = Entry::BLANK;
      ++j;
    }

    // pass over previous stuff
    previous_iterator = iterator;
    previous_column = column;
    previous_index = index;
    
    // get new stuff
    iterator = indexers[i]->get_iterator();
    column = column_iterator.get_next();    
    index = index_path->at(i);

    // db
#ifdef DB

    cout << "i : " << i << endl;
    cout << "---------------------------" << endl << endl;

    cout << "previous index : " << previous_index << endl;
    cout << "index : " << index << endl;
    cout << "dual index : " << iterator->dual_index(index) << endl << endl;

    cout << "previous forward : " << previous_iterator->index_forward_projection(previous_index) << " [" << bit_rep(previous_iterator->index_forward_projection(previous_index),indexers[i-1]->get_forward_projection_width()) << "]" << endl;

    cout << "backward : " << iterator->index_backward_projection(index) << " [" << bit_rep(iterator->index_backward_projection(index),indexers[i]->get_backward_projection_width()) << "]" << endl;
    
    cout << "dual backward : " << iterator->dual_index_backward_projection(index) << " [" << bit_rep(iterator->dual_index_backward_projection(index),indexers[i]->get_backward_projection_width()) << "]" << endl;
    cout << endl;
#endif

    // now the fun part
    if(iterator->index_backward_projection(index) != previous_iterator->index_forward_projection(previous_index)) {
      // index of this column either matches the previous column, or its dual
      index = iterator->dual_index(index);
    }

    ColumnCostComputer cost_computer(*(column.get()), all_heterozygous);
    cost_computer.set_partitioning(index);

    h0[j] = cost_computer.get_allele(0);
    h1[j] = cost_computer.get_allele(1);
    ++j; // next position
  }

  h->push_back(h0);
  h->push_back(h1);

  return h;
}
