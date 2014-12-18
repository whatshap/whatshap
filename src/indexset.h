#ifndef INDEXSET_H
#define INDEXSET_H

#include <set>

/** Class to represent a subset of indices. */
class IndexSet {
public:
	/** Constructor
	 *  @param set_size Gives the size of the full set, that is, the size of the set of which a subset is to be represented.
	 */
	IndexSet();
	virtual ~IndexSet();

	bool contains(std::size_t index) const;

	void add(std::size_t index);

	/** Returns the number of elements in the subset. */
	std::size_t size() const;

	std::string toString();
private:
	std::set<int> set;
};

#endif
