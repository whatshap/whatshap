#ifndef COLUMN_H
#define COLUMN_H

#include <vector>
#include <memory>
#include "columnindexingiterator.h"

class ColumnIndexingIterator;

class Column {
private:
	unsigned int column_size;
	std::vector<unsigned int> read_ids;
	std::vector<unsigned int> next_read_ids;
	std::vector<unsigned int> act_nonterminating_read_ids;
	std::vector<unsigned int> act_terminating_read_ids;
	unsigned int n_references;
	const Column* previous_column;
	const Column* next_column;
	
public:

	/** Constructor.
	 * @param previousReadIDs IDs of reads active
	 */
	Column(const unsigned int index, const unsigned int* n_ref, const std::vector<unsigned int>& read_ids, const std::vector<unsigned int>& next_read_ids);
	
	// returns the column size
	unsigned int get_column_size();

	// return a pointer to the read ids
	std::vector<unsigned int> * get_read_ids();

	// return a pointer to the active non terminating read ids of the column
	std::vector<unsigned int> * get_active_nonterminating_read_ids();

	// return a pointer to the active terminating read ids of the column
	std::vector<unsigned int> * get_active_terminating_read_ids();

	// return a pointer to the active non terminating read ids of the column
	std::vector<unsigned int> * get_next_read_ids();

	// returns a pointer to the bipartition iterator which uses graycode.
	std::unique_ptr<ColumnIndexingIterator> get_iterator();

	// return the bipartition defined by the index
	std::vector<std::vector<unsigned int>> index_to_bipartition(unsigned int& index, int column_type);

	// return the reference alleles defined by the index
	std::vector<unsigned int> index_to_reference_allele(unsigned int& index, int column_type);

	// gets the index value using the bipartition and the reference alleles given.
	unsigned int get_index(std::vector<unsigned int>& b1, std::vector<unsigned int>& b2, unsigned int& r1, unsigned int& r2);

	// gets the index value using the bipartition index and reference index;
	unsigned int get_index(unsigned int b_index, unsigned int r_index);

	// returns the index from bipartition
	unsigned int bipartition_to_index(std::vector<unsigned int>& b1, std::vector<unsigned int>& b2);

	// returns the index from reference alleles
	unsigned int reference_allele_to_index(unsigned int& r1, unsigned int& r2);

	// returns the compatible bipartitions in backward pass
	std::vector<unsigned int> get_backward_compatible_bipartitions(int b_index);

	// returns the compatible bipartitions in forward pass
	std::vector<unsigned int> get_forward_compatible_bipartitions(const int b_index);
};

#endif
