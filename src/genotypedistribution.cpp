#include <cassert>
#include <sstream>
#include <cmath>

#include "genotypedistribution.h"
#include "genotype.h"
#include "binomial.h"

using namespace std;

GenotypeDistribution::GenotypeDistribution(unsigned int nr_allele) {
	n_allele = nr_allele;
	int ploidy = 2;
	int n_genotypes = binomial_coefficient(ploidy + n_allele - 1, n_allele - 1);
	for (int i = 0; i < n_genotypes; i++) {
		distribution.push_back(1.0/(double)n_genotypes);
	}
}


GenotypeDistribution::GenotypeDistribution(std::vector<std::vector<double> > p_wrong_vector, int allele) : distribution() {
	// Here the ploidy has been hardcoded to 2
	int ploidy = 2;
	double prob_sum = 0.0;
	vector<double> tmp_distribution;
	n_allele = p_wrong_vector.size();
	int n_genotypes = binomial_coefficient(ploidy + n_allele - 1, n_allele - 1);
	for (int i = 0; i < n_genotypes; i++) {
		vector<uint32_t> alleles = convert_index_to_alleles(i, ploidy);
		double prob = 1;
		for (int j = 0; j < ploidy; j ++) {
			if (alleles[j] != allele) {
				prob = prob * p_wrong_vector[allele][alleles[j]];
			}
			else {
				prob = prob * (1.0 - p_wrong_vector[allele][alleles[j]]);
			}
		}
		prob_sum = prob_sum + prob;
		distribution.push_back(prob);
	}
	for (int i = 0; i < n_genotypes; i ++) {
		distribution[i] = distribution[i]/prob_sum;
	}
}

GenotypeDistribution::GenotypeDistribution(std::vector<double> d, int nr_allele) : distribution() {
	n_allele = nr_allele;
	for (int i = 0; i < d.size(); i++) {
		distribution.push_back(d[i]);
	}
}


int GenotypeDistribution::likeliestGenotype() const {
	int best_index = 0;
	double best = 0.0;
	for (size_t i=0; i<distribution.size(); ++i) {
		if (distribution[i] > best) {
			best = distribution[i];
			best_index = i;
		}
	}
	return best_index;
}


void GenotypeDistribution::normalize() {
	double p_sum = 0.0;
	for (size_t i=0; i<distribution.size(); ++i) {
		 p_sum += distribution[i];
	}
	if (p_sum <= 0.0) {
		distribution = vector<double>(distribution.size(), 1.0/(double)distribution.size());
	} else {
		for (size_t i=0; i<distribution.size(); ++i) {
			distribution[i] /= p_sum;
		}
	}
}


double GenotypeDistribution::errorProbability() const {
	int best_index = 0;
	double best = 0.0;
	for (size_t i=0; i<distribution.size(); ++i) {
		if (distribution[i] > best) {
			best = distribution[i];
			best_index = i;
		}
	}
	double p = 0.0;
	for (size_t i=0; i<distribution.size(); ++i) {
		if (i == best_index) continue;
		p += distribution[i];
	}
	return p;
}


GenotypeDistribution operator*(const GenotypeDistribution& d1, const GenotypeDistribution& d2) {
	assert (d1.distribution.size() == d2.distribution.size());
	assert (d1.n_allele == d2.n_allele);
	vector<double> d(d1.distribution);
	double sum = 0.0;
	for (int i=0; i<d.size(); ++i) {
		d[i] *= d2.distribution[i];
		sum += d[i];
	}
	for (int i=0; i<d.size(); ++i) d[i] /= sum;
	return GenotypeDistribution(d, d1.n_allele);
}

PhredGenotypeLikelihoods GenotypeDistribution::toPhredLikelihoods() const {
	double max = 0.0;
	for (int i=0; i<distribution.size(); ++i) {
		if (distribution[i] > max) max = distribution[i];
	}
	if (max == 0.0) {
		vector<double> phred_genotyle_likelihood(distribution.size(), 0.0);
		return PhredGenotypeLikelihoods(phred_genotyle_likelihood, 2, n_allele);
	}
	vector<double> phred_genotyle_likelihood;
	for (int i=0; i<distribution.size(); ++i) {
		phred_genotyle_likelihood.push_back(distribution[i]);
	}
	return PhredGenotypeLikelihoods(phred_genotyle_likelihood, 2, n_allele);
}

ostream& operator<<(ostream& os, const GenotypeDistribution& d) {
	os << "(HOM_REF:" << d.probabilityOf(0) << ", HET:" << d.probabilityOf(1) << ", HOM_ALT:" << d.probabilityOf(2) << ")";
	return os;
}
