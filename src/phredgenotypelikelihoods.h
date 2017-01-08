#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <vector>

class PhredGenotypeLikelihoods {
public:
	PhredGenotypeLikelihoods(std::vector<unsigned int>);

	unsigned int get(size_t genotype) const;

	std::vector<unsigned int> get_gl() const;

	std::string toString() const;

private:
	std::vector<unsigned int> gl;
};


#endif
