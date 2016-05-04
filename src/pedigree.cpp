#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include "pedigree.h"

using namespace std;

Pedigree::Pedigree() {
}


void Pedigree::addIndividual(unsigned int id, std::vector<unsigned int> genotypes) {
	this->genotypes.push_back(genotypes);
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
		throw std::runtime_error("Individual with given ID not present in pedigree.");
	}
	return it->second;
}


unsigned int Pedigree::get_genotype(size_t individual_index, size_t column) const {
	return genotypes[individual_index][column];
}

const std::vector<Pedigree::triple_entry_t>& Pedigree::get_triples() const {
	return triples;
}