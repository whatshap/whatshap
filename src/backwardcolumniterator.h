#ifndef BACKWARDCOLUMNITERATOR_H
#define BACKWARDCOLUMNITERATOR_H

#include <vector>
#include <memory>
#include <list>

#include "entry.h"
#include "readset.h"

// TODO create column iterator superclass and subclasses to iterate forwards and backwards

class BackwardColumnIterator {
public:
	BackwardColumnIterator(const ReadSet& set, const std::vector<unsigned int>* positions = nullptr);
	~BackwardColumnIterator();
	/** Returns the total number of columns, i.e. the number of columns
	 *  that will be returned by get_next. */
	unsigned int get_column_count();
	/** Returns the total number of reads. */
	unsigned int get_read_count();
	bool has_next();
	/** Ownership of Entry objects remains with the ColumnIterator. Pointers
	 *  remain valid only until iterator is destructed. */
	std::unique_ptr<std::vector<const Entry*> > get_next();
	const std::vector<unsigned int>* get_positions();
	/** Moves iterator such that next call to get_next() will return
	 *  column k. */
	void jump_to_column(int k);

private:
	typedef struct active_read_t {
		size_t read_index;
		size_t active_entry;
		active_read_t(size_t read_index) : read_index(read_index), active_entry(0) {}
		active_read_t(size_t read_index, size_t active_entry) : read_index(read_index), active_entry(active_entry) {}
	} active_read_t;

	const ReadSet& set;
	/** The number of columns already written. */
	int n;
	std::list<active_read_t> active_reads;
	std::vector<Entry*> blank_entries;
	std::vector<unsigned int>* positions;
	// first_reads[k] is the index of the first read (i.e. lowest index) active at column k,
	// in case no read is active in column k, then first_reads[k] is the index of the first read
	// that will become active after column k.
	std::vector<size_t> first_reads;
};

#endif
