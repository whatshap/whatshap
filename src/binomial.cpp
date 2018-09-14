#include "binomial.h"

int binomial_coeff(int n, int k){
	if (k < 0 || n < 0 || n < k) return 0;
	int result = 1;
	if (k > n-k) k = n-k;

	for (int i = 0; i < k; i++){
		result *= (n-i);
		result /= (i+1);
	}
	return result;
}
