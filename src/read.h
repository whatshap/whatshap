#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>

#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include "entry.h"

class Read {
public:
	Read() : source_id(-1), sample_id(-1), id(-1) {}
	Read(const std::string& name, int mapq, int source_id, int sample_id);
	virtual ~Read() {}
	std::string toString();
	void addVariant(int position, int allele, int quality);
	void sortVariants();
	/** Returns the position of the first variant. **/
	int firstPosition() const;
	/** Returns the position of the last variant. **/
	int lastPosition() const;
	void setID(int id);
	int getID() const;
	/** Add all positions contained in this read to the given set. */
	void addPositionsToSet(std::unordered_set<unsigned int>* set);
	int getPosition(size_t variant_idx) const;
	void setPosition(size_t variant_idx, int position);
	int getAllele(size_t variant_idx) const;
	void setAllele(size_t variant_idx, int allele);
	int getVariantQuality(size_t variant_idx) const;
	void setVariantQuality(size_t variant_idx, int quality);
	const Entry* getEntry(size_t variant_idx) const;
	int getVariantCount() const;
	const std::string& getName() const;
	const std::vector<int>& getMapqs() const;
	void addMapq(int mapq);
	int getSourceID() const;
	int getSampleID() const;
	bool isSorted() const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(name, mapqs, source_id, sample_id, variants);
	}

private:
	typedef struct enriched_entry_t {
		Entry entry;
		int position;
		enriched_entry_t() : position(-1) {}
		enriched_entry_t(int position, int allele, int quality) :
			entry(0,Entry::allele_t(allele),quality), position(position) {}
		template<class Archive>
		void serialize(Archive& ar) {
			ar(position, entry);
		}
	} enriched_entry_t;

	typedef struct entry_comparator_t {
		entry_comparator_t() {}
		bool operator()(const enriched_entry_t& e1, const enriched_entry_t& e2) {
			return e1.position < e2.position;
		}
	} entry_comparator_t;

	std::string name;
	std::vector<int> mapqs;
	int source_id;
	int sample_id;
	int id;
	std::vector<enriched_entry_t> variants;
};

#endif
