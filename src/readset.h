#ifndef READSET_H
#define READSET_H

#include <string>
#include <vector>

#include "read.h"

class ColumnIterator;

class ReadSet {
public:
	ReadSet();
	virtual ~ReadSet();
	/** Ownership of pointer is transferred from caller to the ReadSet. */
	void add(Read* read);
	/** Sorts reads and variants within reads and assigns unique read identifiers. 
	 *  After calling finalize(), the read set becomes "frozen" and cannot be modified. **/
	void finalize();
	/** Returns the set of SNP positions. Can only be called after finalization. */
	const std::vector<unsigned int>* get_positions() const;
	unsigned int size() const;
	std::string toString();
private:
	typedef struct read_comparator_t {
		read_comparator_t() {}
		bool operator()(const Read* r1, const Read* r2) {
			return r1->firstPosition() < r2->firstPosition();
		}
	} read_comparator_t;

	std::vector<Read*> reads;
	std::vector<unsigned int>* positions;
	bool finalized;
	friend class ColumnIterator;
};

#endif