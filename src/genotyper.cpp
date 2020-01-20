
#include <iostream>
#include <cassert>
#include <math.h>

#include "columniterator.h"
#include "genotypedistribution.h"

#include "genotyper.h"

using namespace std;

void compute_genotypes(const ReadSet& readset, std::vector<Genotype>* genotypes, std::vector<GenotypeDistribution>* genotype_likelihoods, std::vector<unsigned int>* positions) {
	assert(genotypes != nullptr);
	assert(genotype_likelihoods != nullptr);
	genotypes->clear();
	genotype_likelihoods->clear();
	if (positions == nullptr) {
		positions = readset.get_positions();
		assert(positions != nullptr);
	}
	ColumnIterator column_iterator(readset, positions);
	size_t column_index = 0;
	while (column_iterator.has_next()) {
// 		cerr << "  working in column " << column_index << ", position " << positions->at(column_index) << endl;
		std::unique_ptr<std::vector<const Entry*> > column = column_iterator.get_next();
		GenotypeDistribution distribution;
		for (const Entry* e : *column) {
// 			cerr << "    " << (*e) << endl;
			double p_wrong = max(0.05, pow(10.0,-((double)e->get_phred_score())/10.0));
			switch (e->get_allele_type()) {
				case Entry::REF_ALLELE:
					distribution = distribution * GenotypeDistribution(2.0/3.0-1.0/3.0*p_wrong, 1.0/3.0, 1.0/3.0*p_wrong);
					break;
				case Entry::ALT_ALLELE:
					distribution = distribution * GenotypeDistribution(1.0/3.0*p_wrong, 1.0/3.0, 2.0/3.0-1.0/3.0*p_wrong);
					break;
				default:
					break;
			}
			
		}
		distribution.normalize();
		Genotype genotype;
		if (distribution.errorProbability() < 0.1) {
			genotype = Genotype(distribution.likeliestGenotype(), Genotype::DIPLOID);
		}
		genotypes->push_back(genotype);
		genotype_likelihoods->push_back(distribution);
// 		cerr << "  --> GT=" << genotype << ", GL=" << distribution << endl;
		++column_index;
	}
	assert(genotypes->size() == positions->size());
	delete positions;
}

void compute_polyploid_genotypes(const ReadSet& readset, size_t ploidy, std::vector<Genotype>* genotypes, std::vector<unsigned int>* positions){
	assert(genotypes != nullptr);
	genotypes->clear();
	if (positions == nullptr) {
		positions = readset.get_positions();
		assert(positions != nullptr);
	}
	ColumnIterator column_iterator(readset, positions);
	size_t column_index = 0;
	while (column_iterator.has_next()) {
		column_index += 1;
		std::unique_ptr<std::vector<const Entry*> > column = column_iterator.get_next();
		size_t ref_count = 0;
		size_t alt_count = 0;
		for (const Entry* e : *column) {
			switch (e->get_allele_type()){
				case Entry::REF_ALLELE:
					ref_count += 1;
					break;
				case Entry::ALT_ALLELE:
					alt_count += 1;
					break;
				default:
					break;
			}
		}
		// TODO: Respect genotype likelihoods
		
		// determine the genotype
		size_t total_alleles = ref_count + alt_count;
		if (total_alleles == 0) {
			genotypes->push_back(Genotype());
		} else {
			double alt_frac = alt_count / (double) (total_alleles);
			uint32_t num_alts = (uint32_t)(ploidy*alt_frac+1/(2*ploidy));
			
			std::vector<uint32_t> geno_alleles;
			for (size_t j = 0; j < num_alts; ++j){
				geno_alleles.push_back(1);
			}
			for (size_t j = num_alts; j < ploidy; ++j){
				geno_alleles.push_back(0);
			}
			
			assert(geno_alleles.size() == ploidy);
			Genotype genotype(geno_alleles);
			genotypes->push_back(genotype);
		}
	}
	cout << genotypes->size() << " " << positions->size() << column_index << endl;
	assert(genotypes->size() == positions->size());
	delete positions;
}
