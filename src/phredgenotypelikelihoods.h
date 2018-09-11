#ifndef PHRED_GENOTYPE_LIKELIHOODS_H
#define PHRED_GENOTYPE_LIKELIHOODS_H

#include <vector>
#include <map>

class PhredGenotypeLikelihoods {
public:
	typedef std::vector<unsigned int> genotype;

	PhredGenotypeLikelihoods(size_t ploidy, size_t n_alleles, const std::vector<double>& gl);

	double get(genotype) const;

	const  std::vector<double>& as_vector() const;

	size_t genotype_count() const;

	std::string toString() const;

private:
	size_t ploidy;
	size_t n_alleles;
	std::vector<double> gl;

	size_t genotype_to_index();
};


#endif
