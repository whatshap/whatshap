#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>

#include "entry.h"

class Read {
public:
	Read(const std::string& name, int mapq, int source_id, int sample_id, int reference_start = -1, const std::string& BX_tag = "");
	virtual ~Read() {}
	std::string toString();
	void addVariant(int position, std::vector<int> alleles, std::vector<int> qualities);
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
	std::vector<int> getAllele(size_t variant_idx) const;
	void setAllele(size_t variant_idx, std::vector<int> allele);
	std::vector<int> getVariantQuality(size_t variant_idx) const;
	void setVariantQuality(size_t variant_idx, std::vector<int> quality);
	const Entry* getEntry(size_t variant_idx) const;
	int getVariantCount() const;
	const std::string& getName() const;
	const std::vector<int>& getMapqs() const;
	void addMapq(int mapq);
	int getSourceID() const;
	int getSampleID() const;
	int getReferenceStart() const;
	const std::string& getBXTag() const;
	bool isSorted() const;
	bool hasBXTag() const;
private:
	typedef struct enriched_entry_t {
		Entry entry;
		int position;
		enriched_entry_t(int position, std::vector<int> allele, std::vector<int> quality) :
			entry(0,std::vector<Entry::allele_t>(),quality), position(position) {
				std::vector<Entry::allele_t> alleles;
				for (auto a : allele) alleles.push_back(Entry::allele_t(a));
				entry.set_allele_type(alleles);
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
	int reference_start;
	std::string BX_tag;
	std::vector<enriched_entry_t> variants;
};

#endif
