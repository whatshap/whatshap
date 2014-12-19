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
	bool isFinalized() const;
	/** Returns the set of SNP positions. Can only be called after finalization. */
	const std::vector<unsigned int>* get_positions() const;
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
			return r1->firstPosition() < r2->firstPosition();
		}
	} read_comparator_t;

	std::vector<Read*> reads;
	// Maps names of reads it their index in the "reads" vector
	typedef std::unordered_map<std::string,size_t> read_name_map_t;
	read_name_map_t read_name_map;
	std::vector<unsigned int>* positions;
	bool finalized;
	friend class ColumnIterator;
};

#endif