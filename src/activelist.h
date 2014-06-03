#ifndef ACTIVE_LIST
#define ACTIVE_LIST

#include <vector>

class ActiveList {
private:
  std::vector<unsigned int> list;
  unsigned int capacity; // the capacity of the list
  unsigned int active_element;
  unsigned int size; // current size of active list
  int push_index; // internal index used by can_push and push
  unsigned int smallest_index; // used for more efficient updating
public:
  /** Constructor.
   *  @param capacity Capacity of the active list.
   */
  ActiveList(unsigned int capacity);

  /** Set the active element
   *
   * NOTE : cannot push until this is set
   */
  void set_active_element(unsigned int active_element);

  /** True iff we can push an element to the list (without exceeding capacity).
   *
   *  The idea is that if the size has already reached capacity, but
   *  there is an element in the list that is smaller than the active
   *  element, we can overwrite it with a new element (which should be
   *  larger than or equal to the active element)
   */
  bool can_push();

  /** Push an element to the list */
  void push(unsigned int element);

  /* Return a const pointer to the list (needed only for debugging) */
  const std::vector<unsigned int> * get_list();

  /* Get the active element (needed only for debugging) */
  unsigned int get_active_element();

  /* Get the current size of the list (needed only for debugging) */
  unsigned int get_size();

  /* Get the  capacity of the list (needed only for debugging) */
  unsigned int get_capacity();

};

#endif
