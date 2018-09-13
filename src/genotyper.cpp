
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
			std::vector<unsigned int> phred_scores = e->get_phred_score();
			double p_wrong = 1.0;
			switch (e->get_allele_type()) {
				case 0:
					assert(phred_scores.size() == 2);
					p_wrong = max(0.05, pow(10.0,-((double)phred_scores[1])/10.0));
					distribution = distribution * GenotypeDistribution(2.0/3.0-1.0/3.0*p_wrong, 1.0/3.0, 1.0/3.0*p_wrong);
					break;
				case 1:
					assert(phred_scores.size() == 2);
					p_wrong = max(0.05, pow(10.0,-((double)phred_scores[0])/10.0));
					distribution = distribution * GenotypeDistribution(1.0/3.0*p_wrong, 1.0/3.0, 2.0/3.0-1.0/3.0*p_wrong);
					break;
				default:
					break;
			}
			
		}
		distribution.normalize();
		Genotype genotype;
		if (distribution.errorProbability() < 0.1) {
			genotype = distribution.likeliestGenotype();
		}
		genotypes->push_back(genotype);
		genotype_likelihoods->push_back(distribution);
// 		cerr << "  --> GT=" << genotype << ", GL=" << distribution << endl;
		++column_index;
	}
	assert(genotypes->size() == positions->size());
	delete positions;
}
