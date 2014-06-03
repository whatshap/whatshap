#include <cassert>
#include "activelistdelegator.h"

using namespace std;

ActiveListDelegator::ActiveListDelegator(unsigned int capacity) {
  this->capacity = capacity;
  this->active_element = 0; // default to 0
  this->size = 0;
  this->index = -1;
}

void ActiveListDelegator::set_active_element(unsigned int active_element) {
  this->active_element = active_element;

  for(unsigned int i=0; i< size; ++i) { // update all active lists
    lists[i].set_active_element(active_element);
  }
}

unsigned int ActiveListDelegator::delegate() {

  if(index >= 0) return index; // in case of repeated calls to delegate()

  for(unsigned int i=0; i< size; ++i) {
    if(lists[i].can_push()) {
      index = i;
      break;
    }
  }

  if(index < 0) index = size; // all lists are full (or lists is empty)
  return index;
}

void ActiveListDelegator::push(unsigned int element) {
  assert(element >= active_element);
  assert(index <= size); // just to be sure delegate didn't fall off the end

  // in case user does not call delegate(), get it
  if(index < 0) this->delegate();

  if(index == size) { // add new list, and push element to new list
    lists.push_back(ActiveList(capacity));
    lists[index].set_active_element(active_element); // update new list
    ++size;
  }

  assert(lists[index].can_push()); // sanity check
  lists[index].push(element); // push the element
  index = -1; // reset
}
