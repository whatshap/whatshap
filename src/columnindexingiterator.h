#ifndef COLUMN_INDEXING_ITERATOR_H
#define COLUMN_INDEXING_ITERATOR_H

#include "graycodes.h"

class ColumnIndexingScheme;

class ColumnIndexingIterator {
private:
	const ColumnIndexingScheme* parent;
	GrayCodes* graycodes;
	unsigned int index;
	unsigned int forward_projection;

public:
	ColumnIndexingIterator(const ColumnIndexingScheme* parent);
	virtual ~ColumnIndexingIterator();

	bool has_next();

	/** Move to next index (i.e. DP table row).
	  *
	  *  @param bit_changed If not null, and only one bit in the
	  *  partitioning (as retrieved by get_partition) is changed by this
	  *  call to advance, then the index of this bit is written to the
	  *  referenced variable; if not, -1 is written.
	  */
	void advance(int* bit_changed = 0);

	/** Index of the projection of the current read set onto the intersection between current and next read set. */
	unsigned int get_forward_projection();

	/** Index of the projection of the current read set onto the intersection between previous and the current read set. */
	unsigned int get_backward_projection();

	/** Row index in the DP table (within the current column). */
	unsigned int get_index();

	/** Bit-wise representation of the partitioning corresponding to the current index. */
	unsigned int get_partition();

	/** get index's backward projection (given index i), so that we don't have to iterate up to it, just to get it */
	unsigned int index_backward_projection(unsigned int i);

	/** get index's forward projection */
	unsigned int index_forward_projection(unsigned int i);
};

#endif
