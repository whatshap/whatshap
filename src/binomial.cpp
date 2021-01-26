#include "binomial.h"
#include <cmath>

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

double binomial_coefficient_large(int n, int k){
	if (k < 0 || n < 0 || n < k) return 0;
	double result = 1.0;
	if (k > n-k) k = n-k;

	for (int i = 0; i < k; i++){
		result *= (n-i);
		result /= (i+1);
	}
	return result;
}

double binom_pmf(int n, int k, double p) {
    return binomial_coefficient_large(n, k) * pow(p, k) * pow(1-p, n-k);
}
