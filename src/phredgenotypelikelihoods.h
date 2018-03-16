#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <vector>

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods();

	PhredGenotypeLikelihoods(const std::vector<unsigned int>&);

	double get(size_t genotype) const;

	const std::vector<unsigned int>& as_vector() const;

	size_t genotype_count() const;

	std::string toString() const;

private:
	std::vector<unsigned int> gl;
};


#endif
