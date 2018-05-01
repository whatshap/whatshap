#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>


#include "read.h"

using namespace std;

Read::Read(const std::string& name, int mapq, int source_id, int sample_id, int reference_start, const std::string& BX_tag) : name(name), mapqs(1, mapq), source_id(source_id), sample_id(sample_id), reference_start(reference_start), BX_tag(BX_tag) {
	this->id = -1;
}


string Read::toString() {
	ostringstream oss;
	oss << name << " mapq:(";
	for (size_t i=0; i<mapqs.size(); ++i) {
		if (i>0) oss << ",";
		oss << mapqs[i];
	}
	oss << ") source:" << source_id << " sample:" << sample_id << " (";
	for (size_t i=0; i<variants.size(); ++i) {
		if (i>0) oss << ";";
		oss << "[" << variants[i].position << "," << variants[i].entry << "]";
	}
	oss << ")";
	return oss.str();
}


void Read::addVariant(int position, int allele, int quality) {
	variants.push_back(enriched_entry_t(position, allele, quality));
}


void Read::sortVariants() {
	sort(variants.begin(), variants.end(), entry_comparator_t());
	for (size_t i=1; i<variants.size(); ++i) {
		if (variants[i-1].position == variants[i].position) {
			ostringstream oss;
			oss << "Duplicate variant in read " << name << " at position " << variants[i].position;
			throw std::runtime_error(oss.str());
		}
	}
}


int Read::firstPosition() const {
	if (variants.size() == 0) throw std::runtime_error("No variants present");
	return variants[0].position;
}


int Read::lastPosition() const {
	if (variants.size() == 0) throw std::runtime_error("No variants present");
	return variants[variants.size()-1].position;
}


void Read::setID(int id) {
	this->id = id;
	for (size_t i=0; i<variants.size(); ++i) {
		variants[i].entry.set_read_id(id);
	}
}


int Read::getID() const {
	return id;
}


void Read::addPositionsToSet(std::unordered_set<unsigned int>* set) {
	assert(set != 0);
	for (size_t i=0; i<variants.size(); ++i) {
		set->insert(variants[i].position);
	}
}


int Read::getPosition(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].position;
}


void Read::setPosition(size_t variant_idx, int position) {
	assert(variant_idx < variants.size());
	variants[variant_idx].position = position;
}


int Read::getAllele(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_allele_type();
}


void Read::setAllele(size_t variant_idx, int allele) {
	assert(variant_idx < variants.size());
	variants[variant_idx].entry.set_allele_type((Entry::allele_t)allele);
}


int Read::getVariantQuality(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_phred_score();
}


void Read::setVariantQuality(size_t variant_idx, int quality) {
	assert(variant_idx < variants.size());
	variants[variant_idx].entry.set_phred_score(quality);
}


const Entry* Read::getEntry(size_t variant_idx) const {
	return &(variants[variant_idx].entry);
}


int Read::getVariantCount() const {
	return variants.size();
}


const string& Read::getName() const {
	return name;
}


const vector<int>& Read::getMapqs() const {
	return mapqs;
}


void Read::addMapq(int mapq) {
	mapqs.push_back(mapq);
}


int Read::getSourceID() const {
	return source_id;
}


int Read::getSampleID() const {
	return sample_id;
}

int Read::getReferenceStart() const {
	return reference_start;
}

const std::string& Read::getBXTag() const {
	return BX_tag; 
}

bool Read::isSorted() const {
	entry_comparator_t comparator;
	for (size_t i=1; i<variants.size(); ++i) {
		if (!comparator(variants[i-1],variants[i])) {
			return false;
		}
	}
	return true;
}

bool Read::hasBXTag() const {
	return (BX_tag != "");
}
