#include <sstream>
#include <algorithm> 
#include <stdexcept>
#include <cassert>

#include "read.h"

using namespace std;

Read::Read(const std::string& name, int mapq) : name(name), mapqs(1, mapq) {
	this->id = -1;
}

string Read::toString() {
	ostringstream oss;
	oss << name << " (";
	for (size_t i=0; i<mapqs.size(); ++i) {
		if (i>0) oss << ",";
		oss << mapqs[i];
	}
	oss << ") (";
	for (size_t i=0; i<variants.size(); ++i) {
		if (i>0) oss << ";";
		oss << "[" << variants[i].position << "," << variants[i].base << "," << variants[i].entry << "]";
	}
	oss << ")";
	return oss.str();
}

void Read::addVariant(int position, char base, int allele, int quality) {
	variants.push_back(enriched_entry_t(position, base, allele, quality));
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

char Read::getBase(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].base;
}

int Read::getAllele(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_allele_type();
}

int Read::getBaseQuality(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_phred_score();
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
