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
	/* TODO: Append an item to triple_genotypes.
	To get the genotypes, use
	genotypes_map.at(m), genotypes_map.at(f), and genotypes_map.at(c).
	*/
}
