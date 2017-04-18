#include <cassert>
#include <sstream>
#include <cmath>

#include "genotypedistribution.h"

using namespace std;

GenotypeDistribution::GenotypeDistribution(double hom_ref_prob, double het_prob, double hom_alt_prob) : distribution() {
	distribution.push_back(hom_ref_prob);
	distribution.push_back(het_prob);
	distribution.push_back(hom_alt_prob);
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
		distribution = vector<double>(3, 1.0/3.0);
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
	vector<double> d(d1.distribution);
	double sum = 0.0;
	for (int i=0; i<3; ++i) {
		d[i] *= d2.distribution[i];
		sum += d[i];
	}
	for (int i=0; i<3; ++i) d[i] /= sum;
	return GenotypeDistribution(d[0], d[1], d[2]);
}

PhredGenotypeLikelihoods GenotypeDistribution::toPhredLikelihoods() const {
	double max = 0.0;
	for (int i=0; i<3; ++i) {
		if (distribution[i] > max) max = distribution[i];
	}
	if (max == 0.0) return PhredGenotypeLikelihoods(0,0,0);
	return PhredGenotypeLikelihoods(
		(int)round(min(255.0,-log10(distribution[0]/max)*10.0)),
		(int)round(min(255.0,-log10(distribution[1]/max)*10.0)),
		(int)round(min(255.0,-log10(distribution[2]/max)*10.0))
	);
}

ostream& operator<<(ostream& os, const GenotypeDistribution& d) {
	os << "(HOM_REF:" << d.probabilityOf(0) << ", HET:" << d.probabilityOf(1) << ", HOM_ALT:" << d.probabilityOf(2) << ")";
	return os;
}
