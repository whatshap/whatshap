#ifndef COLUMNITERATOR_H
#define COLUMNITERATOR_H

#include <vector>
#include <memory>
#include <list>

#include "entry.h"
#include "readset.h"

class ColumnIterator {
public:
	ColumnIterator(const ReadSet& set);
	~ColumnIterator();
	/** Returns the total number of columns, i.e. the number of columns
	 *  that will be returned by get_next. */
	unsigned int get_column_count(); 
	/** Returns the total number of reads. */
	unsigned int get_read_count(); 
	bool has_next();
	/** Ownership of Entry objects remains with the ColumnIterator. Pointers
	 *  remain valid only until iterator is destructed.
	 */
	std::unique_ptr<std::vector<const Entry*> > get_next();
	const std::vector<unsigned int>* get_positions();

private:
	typedef struct active_read_t {
		size_t read_index;
		size_t active_entry;
		active_read_t(size_t read_index) : read_index(read_index), active_entry(0) {}
	} active_read_t;
	
	const ReadSet& set;
	/** The number of columns already written. */
	int n;
	/** Index of the read that is to be examined next. */
	int next_read_index;
	std::list<active_read_t> active_reads;
	std::vector<Entry*> blank_entries;
	std::vector<unsigned int>* positions;
};

#endif
