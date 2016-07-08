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
	/** Sort reads by first variant position. */
	void sort();
	/** Returns the set of SNP positions. To create this set,
	 *  this method iterates over all contained reads.
	 *  Caller owns the returned pointer. */
	std::vector<unsigned int>* get_positions() const;
	unsigned int size() const;
	std::string toString();
	/** Access a read in the set. Ownership stays with the ReadSet. */
	Read* get(int i) const;
	/** Access a read in the set by its name. Ownership stays with the ReadSet. */
	Read* getByName(std::string name, int source_id) const;
	/** Creates a subset of reads as given by the set of indices. Note that this
	 *  creates a COPY of each read.
	 */
	ReadSet* subset(const IndexSet* indices) const;
	/** Assigns read_ids to all instances of Entry stored in the reads such that
	 *  each read_id matches the index of the corresponding read in the ReadSet. */
	void reassignReadIds();
private:
	typedef struct read_comparator_t {
		read_comparator_t() {}
		bool operator()(const Read* r1, const Read* r2) {
			if ((r1->getVariantCount() > 0) || (r2->getVariantCount() > 0)) {
				// put reads with no variants first in the set
				if (r1->getVariantCount() == 0) return true;
				if (r2->getVariantCount() == 0) return false;
				// standard case: sort by positions
				if (r1->firstPosition() != r2->firstPosition()) {
					return r1->firstPosition() < r2->firstPosition();
				}
			}
			// break ties by using hash value
			name_and_source_id_hasher_t hasher;
			std::size_t hash1 = hasher(name_and_source_id_t(r1->getName(), r1->getSourceID()));
			std::size_t hash2 = hasher(name_and_source_id_t(r2->getName(), r2->getSourceID()));
			if (hash1 != hash2) {
				return hash1 < hash2;
			}
			// this is the extremely unlikely case of a hash collision
			// ressort to comparing names and source_ids.
			int name_cmp = r1->getName().compare(r2->getName());
			if (name_cmp != 0) {
				return name_cmp < 0;
			}
			return r1->getSourceID() < r2->getSourceID();
		}
	} read_comparator_t;

	typedef struct name_and_source_id_t {
		name_and_source_id_t(std::string name, int source_id) : name(name), source_id(source_id) {}
		bool operator==(const name_and_source_id_t& other) const {
			return (name.compare(other.name) == 0) && (source_id == other.source_id);
		}
		std::string name;
		int source_id;
	} name_and_source_id_t;

	typedef struct name_and_source_id_hasher_t {
		std::size_t operator()(const name_and_source_id_t& x) const {
			return (std::hash<std::string>()(x.name)) ^ (std::hash<int>()(x.source_id));
		}
	} name_and_source_id_hasher_t;

	std::vector<Read*> reads;
	// Maps names of reads it their index in the "reads" vector
	typedef std::unordered_map<name_and_source_id_t,size_t,name_and_source_id_hasher_t> read_name_map_t;
	read_name_map_t read_name_map;
};

#endif