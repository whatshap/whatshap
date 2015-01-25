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
}

ColumnIndexingIterator::~ColumnIndexingIterator() {
  delete graycodes;
}

bool ColumnIndexingIterator::has_next() {
  return graycodes->has_next();
}

// auxilliary function that returns the complement of unsigned int i
// with respect to a width, should its most significant bit be 1
// (i.e., it is 1 at bit index width), otherwise it returns i
unsigned int wrapped(unsigned int i, unsigned int width) {

  if((i >> (width-1)) == 1)
    return ~i & ((((unsigned int)1) << width) -1);

  return i;
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
  }
  else {

    if (parent->forward_projection_mask != 0) {
      // index of bit in the forward_projection
      int bit_index = parent->forward_projection_mask->at(graycode_bit_changed);
      if (bit_index >= 0) {
	forward_projection = forward_projection ^ (((unsigned int)1) << bit_index);
      }
      forward_projection = wrapped(forward_projection, parent->forward_projection_width);
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

unsigned int ColumnIndexingIterator::get_backward_projection() {
  assert(index >= 0);
  return wrapped(index & ((((unsigned int)1)<<parent->backward_projection_width) - 1), parent->backward_projection_width);
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

  return wrapped(i & ((((unsigned int)1) << parent->backward_projection_width) -1), parent->backward_projection_width);
}
