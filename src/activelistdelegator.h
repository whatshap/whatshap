#ifndef ACTIVE_LIST_DELEGATOR
#define ACTIVE_LIST_DELEGATOR

#include <vector>
#include "activelist.h"

class ActiveListDelegator {
private:
  std::vector<ActiveList> lists; // active lists to which it delegates
  unsigned int capacity;
  unsigned int active_element;
  unsigned int size; // number of active lists
  int index; // internal index used by delegate and push
public:
  /** Constructor.
   *  @param capacity Capacities of the active lists to which it delegates
   */
  ActiveListDelegator(unsigned int capacity);

  /** Set the active element */
  void set_active_element(unsigned int active_element);

  /** The id \in (0,size-1) of the smallest list that has availability
   *  (i.e., for which can_push() is true).  If no such list exists,
   *  it returns size (the id of the list that will be created in
   *  order to push this new element).
   */
  unsigned int delegate();

  /** Push an element to the list of id given by delegate(), creating a new
   *  list at the end of lists and pushing to that list, in the case
   *  that id is size.
   */
  void push(unsigned int element);
};

#endif
