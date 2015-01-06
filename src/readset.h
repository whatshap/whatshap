#ifndef READSET_H
#define READSET_H

#include <string>
#include <vector>
#include <unordered_map>

#include "read.h"
#include "indexset.h"

class ColumnIterator;

class ReadSet {
public:
	ReadSet();
	virtual ~ReadSet();
	/** Ownership of pointer is transferred from caller to the ReadSet. */
	void add(Read* read);
	/** Only retains reads with at least two variants, sorts reads and variants within reads
	 *  and assigns unique read identifiers. After calling finalize(), the read set becomes 
	 * "frozen" and cannot be modified. 
	 **/
	void finalize();
	/** Sort reads by first variant position. */
	void sort();
	bool isFinalized() const;
	/** Returns the set of SNP positions. To create this set,
	 *  this method iterates over all contained reads.
	 *  Caller owns the returned pointer. */
	std::vector<unsigned int>* get_positions() const;
	unsigned int size() const;
	std::string toString();
	/** Access a read in the set. Ownership stays with the ReadSet. */
	Read* get(int i) const;
	/** Access a read in the set by its name. Ownership stays with the ReadSet. */
	Read* getByName(std::string name) const;
	/** Creates a subset of reads as given by the set of indices. Note that this
	 *  creates a COPY of each read.
	 */
	ReadSet* subset(const IndexSet* indices) const;
private:
	typedef struct read_comparator_t {
		read_comparator_t() {}
		bool operator()(const Read* r1, const Read* r2) {
			// if both reads don't have variants, then ressort to compare names
			if ((r1->getVariantCount() == 0) && (r2->getVariantCount() == 0)) {
				return r1->getName() < r2->getName();
			}
			// put reads with no variants first in the set
			if (r1->getVariantCount() == 0) return true;
			if (r2->getVariantCount() == 0) return false;
			// standard case: sort by positions
			if (r1->firstPosition() != r2->firstPosition()) {
				return r1->firstPosition() < r2->firstPosition();
			} else {
				// otherwise ressort to comparing name
				return r1->getName() < r2->getName();
			}
		}
	} read_comparator_t;

	std::vector<Read*> reads;
	// Maps names of reads it their index in the "reads" vector
	typedef std::unordered_map<std::string,size_t> read_name_map_t;
	read_name_map_t read_name_map;
	bool finalized;
};

#endif