#include <cassert>
#include "columnindexingscheme.h"
#include "columnindexingiterator.h"

#include<iostream>
#include<math.h>

using namespace std;

ColumnIndexingIterator::ColumnIndexingIterator(const ColumnIndexingScheme* parent, unsigned int number_of_partitions) {
	assert(parent != 0);
	this->parent = parent;
	this->graycodes = new GrayCodes(parent->read_ids.size(), number_of_partitions);
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
			int read_index = parent->forward_projection_mask->at(graycode_position_changed);
			if (read_index >= 0) {
				// change corresponding position in forward projection and compute decimal number
				forward_projection = switch_read(forward_projection, read_index, graycode_partition_changed, parent->forward_projection_size());
			}
		}
	}
	if (position_changed != 0){
		*position_changed = graycode_position_changed;
	}
	if (partition_changed != 0){
		*partition_changed = graycode_partition_changed;
	}
}


unsigned int ColumnIndexingIterator::get_forward_projection() {
	assert(index >= 0);
	return forward_projection;
}


unsigned int ColumnIndexingIterator::get_backward_projection() {
	assert(index >= 0);
	unsigned int steps = parent->backward_projection_width;
	unsigned int digits = index;
	unsigned int result = 0;
	unsigned int factor = 1;
	while(steps > 0){
		result += (digits % number_of_partitions) * factor;
		digits /= number_of_partitions;
		factor *= number_of_partitions;
		steps -= 1;	 
	}
	return result;
}


unsigned int ColumnIndexingIterator::get_index() {
	assert(index >= 0);
	return index;
}


unsigned int ColumnIndexingIterator::get_partition() {
	assert(index >= 0);
	return index;
}


unsigned int ColumnIndexingIterator::index_backward_projection(unsigned int i) {
	assert(i >= 0);
	assert(i < pow(number_of_partitions, parent->read_ids.size()));

	unsigned int steps = parent->backward_projection_width;
	unsigned int digits = i;
	unsigned int result = 0;
	unsigned int factor = 1;

	while(steps > 0){
		result += (digits % number_of_partitions) * factor;
		digits /= number_of_partitions;
		factor *= number_of_partitions;
		steps -= 1;
	}
	return result;
}


unsigned int ColumnIndexingIterator::index_forward_projection(unsigned int i) {
	assert(i >= 0);
	assert(i < pow(number_of_partitions, parent->read_ids.size()));
	
	unsigned int result = 0;
	unsigned int factor = 1;
	for(unsigned int j = 0; j < parent->read_ids.size(); j++){
		unsigned int m = parent->forward_projection_mask->at(j);
		if(m != -1){
			// get the partition the read is assigned to
			unsigned int partition =  (i / (unsigned int) pow(number_of_partitions, j)) % number_of_partitions;
			result += partition * factor;
			factor *= number_of_partitions;
		}
	}
	return result;
}


unsigned int ColumnIndexingIterator::switch_read(unsigned int old_index, unsigned int read_to_switch, unsigned int new_partition, unsigned int used_bits) {
	unsigned int factor = pow(number_of_partitions, read_to_switch);

	// get the partition the read is assigned to currently
	unsigned int old_partition = (old_index / factor) % number_of_partitions;
	unsigned int tmp = old_index - (old_partition * factor);

	return tmp + new_partition * factor;
}
