#include <cassert>
#include "column.h"
#include "columnindexingiterator.h"

using namespace std;

ColumnIndexingIterator::ColumnIndexingIterator(Column* parent) {
	assert(parent != 0);
	this->parent = parent;
	this->graycodes = new GrayCodes(parent->get_read_ids()->size());
	this->b_index = -1;

}


ColumnIndexingIterator::~ColumnIndexingIterator() {
	delete graycodes;
}


bool ColumnIndexingIterator::has_next() {
	return graycodes->has_next();
}


void ColumnIndexingIterator::advance(int* bit_changed) {
	assert(graycodes->has_next());
	int graycode_bit_changed = -1;
	b_index = graycodes->get_next(&graycode_bit_changed);
	binaryVector = graycodes->toBinary(b_index);
	if (bit_changed != 0) {
		*bit_changed = graycode_bit_changed;
	}
}

unsigned int ColumnIndexingIterator::get_b_index() {
	return this->b_index;
}

vector<int> ColumnIndexingIterator::get_binary_vector() const {
	return this->binaryVector;
}
