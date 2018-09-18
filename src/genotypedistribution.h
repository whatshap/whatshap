#ifndef GENOTYPEDISTRIBUTION_H
#define GENOTYPEDISTRIBUTION_H

#include <vector>

#include "genotypelikelihoods.h"
#include "genotype.h"

class GenotypeDistribution {
private:
	std::vector<double> distribution;
public:
	
	GenotypeDistribution() : distribution(3, 1.0/3.0) {}
	GenotypeDistribution(double hom_ref_prob, double het_prob, double hom_alt_prob);

	double probabilityOf(Genotype genotype) const {
		unsigned int index = genotype.get_index(2,2);
		assert(index <= 2);
		return distribution[index];
	}

	/** Normalize distribution such that it sums to one. */
	void normalize();

	Genotype likeliestGenotype() const;

	/** Returns the probability that the likeliest genotype is wrong. */
	double errorProbability() const;

	GenotypeLikelihoods toPhredLikelihoods() const;

	friend GenotypeDistribution operator*(const GenotypeDistribution& d1, const GenotypeDistribution& d2);
};

std::ostream& operator<<(std::ostream& os, const GenotypeDistribution& d);

#endif // GENOTYPEDISTRIBUTION_H
