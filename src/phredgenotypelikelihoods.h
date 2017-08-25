#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <array>

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods(double gl0 = 0, double gl1 = 0, double gl2 = 0);

	double get(size_t genotype) const;

	std::string toString() const;

private:
	std::array<double, 3> gl;
};


#endif
