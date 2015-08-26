#include <unistd.h>
#include <cassert>
#include "columnindexingscheme.h"
#include "columnindexingiterator.h"

using namespace std;

ColumnIndexingIterator::ColumnIndexingIterator(const ColumnIndexingScheme* parent) {
  assert(parent != 0);
  this->parent = parent;
  this->graycodes = new GrayCodes(parent->read_ids.size());
  this->index = -1;
  this->forward_projection = -1;
  this->forward_dual_projection = -1; // used for HALF_TABLE
}

ColumnIndexingIterator::~ColumnIndexingIterator() {
  delete graycodes;
}

bool ColumnIndexingIterator::has_next() {
  return graycodes->has_next();
}

// mt
int ColumnIndexingIterator::get_length() const { 
  return graycodes->get_length();
}
unsigned int ColumnIndexingIterator::compute_forward_projection() {
    unsigned int fp =0;
    //if (parent->forward_projection_mask == NULL) return -1;
    for(size_t i=0;i<parent->forward_projection_mask->size(); ++i) {
	if ((int)(parent->forward_projection_mask->at(i)) >= 0) 
	    fp += ((((1<<i)&index)?1:0)*(1<< parent->forward_projection_mask->at(i)));
    }
    return fp;
}

void ColumnIndexingIterator::reset(const ColumnIndexingScheme* parent) {
  assert(parent != 0);
  this->parent = parent;
  if (graycodes == NULL)
      graycodes = new GrayCodes(parent->read_ids.size());
  else
      graycodes->reset(parent->read_ids.size());
  index = -1;
  forward_projection = -1;
  forward_dual_projection = -1; // used for HALF_TABLE
}


void ColumnIndexingIterator::advance_idx(int *bit_changed, const int idx) {
    //assert(graycodes->has_next());  COMMENTATA

  int graycode_bit_changed = -1;
  index = graycodes->get_next(&graycode_bit_changed, idx);
  // first iteration?
  if (graycode_bit_changed == -1) {
      if (parent->forward_projection_mask != 0) {
	  forward_projection = compute_forward_projection();
      }
  } else {
    if (parent->forward_projection_mask != 0) {

      // index of bit in the forward_projection
      int bit_index = parent->forward_projection_mask->at(graycode_bit_changed);
      if (bit_index >= 0) {
	  forward_projection = forward_projection ^ (((unsigned int)1) << bit_index);
	  //forward_dual_projection = forward_dual_projection ^ (((unsigned int)1) << bit_index);
      }
    }
  }
  if (bit_changed != 0) {
    *bit_changed = graycode_bit_changed;
  }
}

void ColumnIndexingIterator::advance(int* bit_changed) {
  assert(graycodes->has_next());

  int graycode_bit_changed = -1;
  index = graycodes->get_next(&graycode_bit_changed);
  // first iteration?
  if (graycode_bit_changed == -1) {
    assert(index == 0);
    if (parent->forward_projection_mask != 0) {
	forward_projection = 0;
    }
  } else {
    if (parent->forward_projection_mask != 0) {

      // index of bit in the forward_projection
      int bit_index = parent->forward_projection_mask->at(graycode_bit_changed);
      if (bit_index >= 0) {
	  forward_projection = forward_projection ^ (((unsigned int)1) << bit_index);
	  //forward_dual_projection = forward_dual_projection ^ (((unsigned int)1) << bit_index);
      }
    }
  }
  if (bit_changed != 0) {
    *bit_changed = graycode_bit_changed;
  }
}

unsigned int ColumnIndexingIterator::get_forward_projection() {
  assert(index >= 0);
  return forward_projection;
}

unsigned int ColumnIndexingIterator::get_forward_dual_projection() {
  assert(index >= 0);
  return forward_dual_projection;
}

unsigned int ColumnIndexingIterator::get_backward_projection() {
  assert(index >= 0);
  return index & ((((unsigned int)1)<<parent->backward_projection_width) - 1);
}

unsigned int ColumnIndexingIterator::get_index() {
  assert(index >= 0);
  return index;
}

unsigned int ColumnIndexingIterator::get_partition() {
  assert(index >= 0);
  return index;
}

// for backtracking

unsigned int ColumnIndexingIterator::index_backward_projection(unsigned int i) {
  assert(i >= 0); // assert the proper boundaries
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  return i & ((((unsigned int)1) << parent->backward_projection_width) -1);
}

// in the HALF_TABLE case

unsigned int ColumnIndexingIterator::dual_index_backward_projection(unsigned int i) {
  assert(i >= 0); // assert the proper boundaries
  assert(i < (((unsigned int)1) << parent->read_ids.size()));
  
  // return backward projection of the dual
  return dual_index(i) & ((((unsigned int)1) << parent->backward_projection_width) -1);
}

unsigned int ColumnIndexingIterator::dual_index(unsigned int i) {
  assert(i >= 0);
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  unsigned int dual_i = i;
  unsigned int s = 1;
  for(int j=0; j< parent->read_ids.size(); ++j) { // compute the dual of i
    dual_i = dual_i ^ s;
    s = s << 1;
  }

  return dual_i;
}

unsigned int ColumnIndexingIterator::index_forward_projection(unsigned int i) {
  assert(i >= 0);
  assert(i < (((unsigned int)1) << parent->read_ids.size()));

  unsigned int i_forward_projection = 0;
  for(size_t j=0; j< parent->read_ids.size(); ++j) {
    unsigned int m = parent->forward_projection_mask->at(j);
    if((int)m != -1) {
      unsigned int s = (((unsigned int)1) << m);
      i_forward_projection += (s&i);
    }
  }
  
  return i_forward_projection;
}
