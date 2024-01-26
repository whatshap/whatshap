#include "binomial.h"
#include <cmath>
#include <limits>

int binomial_coefficient(int n, int k){
	if (k < 0 || n < 0 || n < k) return 0;
	int result = 1;
	if (k > n-k) k = n-k;

	for (int i = 0; i < k; i++){
		result *= (n-i);
		result /= (i+1);
	}
	return result;
}

double binomial_coefficient_log(int n, int k){
	if (k < 0 || n < 0 || n < k) return 0;
	double result = 0.0;
	if (k > n-k) k = n-k;

	double buffer = 1.0; // product to collect coefficients
	for (int i = 0; i < k; i++){
		double addition = (double)(n-i) / (double)(i+1);
		// avoid to many log operations, only log when current product overflows
		if (buffer * addition > std::numeric_limits<double>::max()) {
			result += std::log(buffer);
			buffer = addition;
		} else {
			buffer *= addition;
		}
	}
	return result + std::log(buffer);
}

double binom_pmf(int n, int k, double p) {
    return std::exp(binomial_coefficient_log(n, k)) * pow(p, k) * pow(1 - p, n - k);
}

double log_binom_pmf(int n, int k, double p) {
    return binomial_coefficient_log(n, k) + k * std::log(p) + (n - k) * std::log(1 - p);
}
