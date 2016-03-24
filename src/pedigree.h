#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>

class Pedigree {

    public:
      typedef std::vector<unsigned int> triple_entry_t;
      typedef std::array<std::vector<unsigned int> , 3> genotype_entry_t;
      triple_entry_t triple_entry_t;
      genotype_entry_t genotype_entry_t;

      Pedigree( std::vector<triple_entry_t> triples, std::vector<genotype_entry_t> genotypes);
      std::vector<triple_entry_t> triples;
      std::vector<genotype_entry_t> genotypes;

};



