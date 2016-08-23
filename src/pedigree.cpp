#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include "pedigree.h"

using namespace std;

Pedigree::Pedigree() {
	variant_count = -1;
}


Pedigree::~Pedigree() {
	for (size_t i=0; i<individual_ids.size(); ++i) {
		for (size_t j=0; j<genotype_likelihoods[i].size(); ++j) {
			if (genotype_likelihoods[i][j] != nullptr) {
				delete genotype_likelihoods[i][j];
			}
		}
	}
}


void Pedigree::addIndividual(unsigned int id, std::vector<unsigned int> genotypes, std::vector<PhredGenotypeLikelihoods*> genotype_likelihoods) {
	if (variant_count == -1) {
		variant_count = genotypes.size();
	}
	assert(genotypes.size() == variant_count);
	assert(genotype_likelihoods.size() == variant_count);
	this->genotypes.push_back(genotypes);
	this->genotype_likelihoods.push_back(genotype_likelihoods);
	individual_ids.push_back(id);
	id_to_index_map[id] = individual_ids.size() - 1;
}


void Pedigree::addRelationship(unsigned int mother_id, unsigned int father_id, unsigned int child_id) {
	triple_entry_t triple_entry = {id_to_index(mother_id), id_to_index(father_id), id_to_index(child_id)};
	triples.push_back(triple_entry);
}


size_t Pedigree::id_to_index(unsigned int individual_id) const {
	auto it = id_to_index_map.find(individual_id);
	if (it == id_to_index_map.end()) {
	    ostringstream oss;
	    oss << "Individual with ID " << individual_id << " not present in pedigree.";
		throw std::runtime_error(oss.str());
	}
	return it->second;
}


unsigned int Pedigree::index_to_id(size_t individual_index) const {
	assert(individual_index < individual_ids.size());
	return individual_ids[individual_index];
}


unsigned int Pedigree::get_genotype(size_t individual_index, size_t variant_index) const {
	return genotypes[individual_index][variant_index];
}


unsigned int Pedigree::get_genotype_by_id(unsigned int individual_id, unsigned int variant_index) const {
	assert(variant_index < variant_count);
	return get_genotype(id_to_index(individual_id), variant_index);
}


const PhredGenotypeLikelihoods* Pedigree::get_genotype_likelihoods(size_t individual_index, size_t variant_index) const {
	return genotype_likelihoods[individual_index][variant_index];
}


const PhredGenotypeLikelihoods* Pedigree::get_genotype_likelihoods_by_id(unsigned int individual_id, unsigned int variant_index) const {
	assert(variant_index < variant_count);
	return get_genotype_likelihoods(id_to_index(individual_id), variant_index);
}


const std::vector<Pedigree::triple_entry_t>& Pedigree::get_triples() const {
	return triples;
}


std::string Pedigree::toString() const {
	ostringstream oss;
	oss << "Pedigree:" << endl;
	oss << "  individuals (index,id):";
	for (size_t i=0; i<individual_ids.size(); ++i) {
		oss << " " << i << "," << individual_ids[i];
	}
	oss << endl;
	oss << "  triples by index (mother,father,child):";
	for (size_t i=0; i<triples.size(); ++i) {
		oss << " (" << triples[i][0] << ","  << triples[i][1] << "," << triples[i][2] << ")";
	}
	oss << endl;
	oss << "  triples by id (mother,father,child):";
	for (size_t i=0; i<triples.size(); ++i) {
		oss << " (" << individual_ids[triples[i][0]] << ","  << individual_ids[triples[i][1]] << "," << individual_ids[triples[i][2]] << ")";
	}
	oss << endl;
	oss << "  genotypes (and likelihoods):" << endl;
	for (size_t i=0; i<individual_ids.size(); ++i) {
		oss << "    individual index:" << i << " / id:" << individual_ids[i] << ":" << endl;
		for (size_t j=0; j<variant_count; ++j) {
			oss << "      " << genotypes[i][j] << "(GL:";
			if (genotype_likelihoods[i][j] == nullptr) {
				oss << "None)" << endl;;
			} else {
				oss << genotype_likelihoods[i][j]->toString();
			}
		}
	}
	
	return oss.str();

}
