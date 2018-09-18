#ifndef GENOTYPE_LIKELIHOODS_H
#define GENOTYPE_LIKELIHOODS_H

#include <vector>
#include <set>
#include "genotype.h"

class GenotypeLikelihoods {
public:
	GenotypeLikelihoods(unsigned int ploidy, unsigned int n_alleles, const std::vector<double>& gl, bool is_phred_scaled);

	double get(Genotype genotype) const;

	const  std::vector<double>& as_vector() const;

	size_t genotype_count() const;

	std::string toString() const;

	unsigned int get_ploidy() const;

	unsigned int get_n_alleles() const;

	bool is_phred() const;

	Genotype get_likeliest_genotype(double threshold_prob) const;

private:
	unsigned int ploidy;
	unsigned int n_alleles;
	std::vector<double> gl;
	bool is_phred_scaled;
};

#endif
