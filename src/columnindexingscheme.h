#ifndef COLUMN_INDEXING_SCHEME_H
#define COLUMN_INDEXING_SCHEME_H

#include <vector>
#include <memory>

class ColumnIndexingIterator;

class ColumnIndexingScheme {
private:
	std::vector<unsigned int> read_ids;
	const ColumnIndexingScheme* previous_column;
	const ColumnIndexingScheme* next_column;
	unsigned int backward_projection_width;
	unsigned int forward_projection_width;
	std::vector<unsigned int>* forward_projection_mask;

public:

	/** Constructor.
	 * @param previousReadIDs IDs of reads active
	 */
	ColumnIndexingScheme(const ColumnIndexingScheme* previous_column, const std::vector<unsigned int>& read_ids);

	~ColumnIndexingScheme();

	/** Set pointer to indexing scheme of next column. MUST be called before the getIterator
	 *  method is called. */
	void set_next_column(const ColumnIndexingScheme* next_column);

	std::unique_ptr<ColumnIndexingIterator> get_iterator();

	unsigned int column_size();

	unsigned int forward_projection_size();
 
	unsigned int get_forward_projection_width();

	unsigned int get_backward_projection_width();

	// return a const pointer to the read ids
	const std::vector<unsigned int> * get_read_ids();

	// return const forward projection mask (for debugging)
	const std::vector<unsigned int> * get_forward_projection_mask();

	friend class ColumnIndexingIterator;
};

#endif
