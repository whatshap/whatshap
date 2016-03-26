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
	triple_entry.clear();
	triple_entry.push_back(m);
	triple_entry.push_back(f);
	triple_entry.push_back(c);
	triples.push_back(triple_entry);
	
	//add genotypes now
	genotype_entry[0]=genotypes_map.at(m);
	genotype_entry[1]=genotypes_map.at(f);
	genotype_entry[2]=genotypes_map.at(c);
	triple_genotypes.push_back(genotype_entry);
}
