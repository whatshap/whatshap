#include <cassert>
#include "activelist.h"

using namespace std;

ActiveList::ActiveList(unsigned int capacity) : list(capacity,0) {
  this->capacity = capacity;
  this->active_element = 0; // default to 0
  this->size = 0;
  this->push_index = -1;
  this->smallest_index = 0;
}

void ActiveList::set_active_element(unsigned int active_element) {
  this->active_element = active_element;
}

bool ActiveList::can_push() {

  if(push_index >= 0) return true; // in case of repeated calls to can_push()

  if(size < capacity) return true;
  else {
    if(active_element <= list[smallest_index]) return false;
    else {
      for(unsigned int i=0; i< size; ++i) {
	if(list[i] < active_element) {
	  push_index = i;
	  return true;
	}
      }
    }
  }
  return false;
}

void ActiveList::push(unsigned int element) {
  assert(element >= active_element);

  if(size < capacity) {
    list[size] = element;
    if(list[size] < list[smallest_index]) smallest_index = size;
    ++size;
  }
  else {
    // asserting something to get its sideaffect ... hmm
    if(push_index < 0) { assert(this->can_push()); }
    list[push_index] = element;

    if(push_index == smallest_index) { // update smallest element
      for(unsigned int i = push_index; i< size; ++i) {
	if(list[i] < list[smallest_index]) smallest_index = i;
      }
    }
    push_index = -1; // reset
  }
}

const std::vector<unsigned int> * ActiveList::get_list() {
  return(&list);
}

unsigned int ActiveList::get_active_element() {
  return active_element;
}

unsigned int ActiveList::get_size() {
  return size;
}

unsigned int ActiveList::get_capacity() {
  return capacity;
}
