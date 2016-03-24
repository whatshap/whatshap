#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include "pedigree.h"

using namespace std;

Pedigree::Pedigree(const std::vector<triple_entry_t> triples, const std::vector<genotype_entry_t> genotypes) : triples(std::move(triples)), genotypes(std::move(genotypes)) {
}

void Pedigree::addTriple(unsigned int m, unsigned int f, unsigned int c) {
	triple_entry_t.clear();
	triple_entry_t.push_back(m);
	triple_entry_t.push_back(f);
	triple_entry_t.push_back(c);
	triples.push_back(triple_entry_t);
}

void Pedigree::addGenoTriple(std::vector<unsigned int>& m, std::vector<unsigned int>& f, std::vector<unsigned int>& c) {
	genotype_entry_t[0]= m;
	genotype_entry_t[1]=f;
	genotype_entry_t[2]=c;
	genotypes.push_back(genotype_entry_t);
}