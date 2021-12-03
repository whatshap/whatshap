#ifndef GENOTYPEDISTRIBUTION_H
#define GENOTYPEDISTRIBUTION_H

#include <vector>

#include "phredgenotypelikelihoods.h"

class GenotypeDistribution {
private:
	std::vector<double> distribution;
	int n_allele;
public:
	
	GenotypeDistribution() : distribution() {}
	GenotypeDistribution(unsigned int nr_allele);
	GenotypeDistribution(std::vector<std::vector<double> > p_wrong_vector, int allele);
	GenotypeDistribution(std::vector<double> d, int nr_allele);

	double probabilityOf(size_t genotype) const {
		return distribution[genotype];
	}

	int getSize() const {
		return distribution.size();
	}

	/** Normalize distribution such that it sums to one. */
	void normalize();

	int likeliestGenotype() const;

	/** Returns the probability that the likeliest genotype is wrong. */
	double errorProbability() const;

	PhredGenotypeLikelihoods toPhredLikelihoods() const;

	friend GenotypeDistribution operator*(const GenotypeDistribution& d1, const GenotypeDistribution& d2);
};

std::ostream& operator<<(std::ostream& os, const GenotypeDistribution& d);

#endif // GENOTYPEDISTRIBUTION_H
