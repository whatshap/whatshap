#include <cassert>
#include "columnindexingscheme.h"
#include "columnindexingiterator.h"

#include<iostream>
#include<math.h>

using namespace std;

ColumnIndexingIterator::ColumnIndexingIterator(const ColumnIndexingScheme* parent, unsigned int number_of_partitions) {
	assert(parent != 0);
	this->parent = parent;
	this->graycodes = new GeneralizedGrayCodes(parent->read_ids.size(), number_of_partitions);
	this->index = -1;
	this->forward_projection = -1;
	this->number_of_partitions = number_of_partitions;
}


ColumnIndexingIterator::~ColumnIndexingIterator() {
	delete graycodes;
}


bool ColumnIndexingIterator::has_next() {
	return graycodes->has_next();
}


void ColumnIndexingIterator::advance(int* position_changed, int* partition_changed) {
	assert(graycodes->has_next());
	
	int graycode_position_changed = -1;
	int graycode_partition_changed = -1;
	index = graycodes->get_next(&graycode_position_changed, &graycode_partition_changed);

	if (graycode_position_changed == -1) {
		assert(index == 0);
		if (parent->forward_projection_mask != 0) {
			forward_projection = 0;
		}
	} else {
		assert(graycode_position_changed >= 0);
		if (parent->forward_projection_mask != 0) {
			// get index (wrt to all reads in intersection) of read that switched partitions
			int read_index = parent->forward_projection_mask->at(parent->read_ids.size() - graycode_position_changed - 1);
			if (read_index >= 0) {
				// change corresponding position in forward projection and compute decimal number
				int old = forward_projection;
				forward_projection = switch_read(forward_projection, read_index, graycode_partition_changed, parent->forward_projection_size());
			}
		}
	}
}


unsigned int ColumnIndexingIterator::get_forward_projection() {
	assert(index >= 0);
	return forward_projection;
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

// TODO implement
unsigned int ColumnIndexingIterator::index_backward_projection(unsigned int i) {
	throw std::runtime_error("Not yet implemented.");
/**
	assert(i >= 0); // assert the proper boundaries
	assert(i < (((unsigned int)1) << parent->read_ids.size()));

	return i & ((((unsigned int)1) << parent->backward_projection_width) -1);
**/
}

// TODO implement
unsigned int ColumnIndexingIterator::index_forward_projection(unsigned int i) {
	throw std::runtime_error("Not yet implemented.");
/**	assert(i >= 0);
	assert(i < (((unsigned int)1) << parent->read_ids.size()));

	unsigned int i_forward_projection = 0;
	unsigned int s = 1;
	for(int j=0; j< parent->read_ids.size(); ++j) {
		unsigned int m = parent->forward_projection_mask->at(j);
		if(m != -1) {
			unsigned int s = (((unsigned int)1) << m);
			i_forward_projection += (s&i);
		}
	}

	return i_forward_projection;
**/
}


unsigned int ColumnIndexingIterator::switch_read(unsigned int old_index, unsigned int read_to_switch, unsigned int new_partition, unsigned int used_bits) {
	unsigned int result = 0;
	unsigned int factor = 1;
	unsigned int i = 0;

	while (used_bits > 0) {
		unsigned int digit = old_index % number_of_partitions;
		if (i == read_to_switch) digit = new_partition;
		result += digit * factor;
		factor *= number_of_partitions;
		i += 1;
		old_index /= number_of_partitions;
		used_bits -= 1;
	}
	return result;
}
