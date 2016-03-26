#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include "pedigree.h"

using namespace std;

Pedigree::Pedigree() {
}


void Pedigree::addIndividual(unsigned int id, std::vector<unsigned int> genotypes) {
	genotypes_map[id] = genotypes;
}


void Pedigree::addRelationship(unsigned int m, unsigned int f, unsigned int c) {
	triple_entry_t triple_entry = {m, f, c};
	triples.push_back(triple_entry);
	genotype_entry_t genotype_entry = {
		genotypes_map.at(m),
		genotypes_map.at(f),
		genotypes_map.at(c)
	};
	triple_genotypes.push_back(genotype_entry);
}
