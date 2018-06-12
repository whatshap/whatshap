#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <vector>

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods(const std::vector<double>& gl);

	double get(size_t genotype) const;

	const  std::vector<double>& as_vector() const;

	size_t genotype_count() const;

	std::string toString() const;

private:
	std::vector<double> gl;
};


#endif
