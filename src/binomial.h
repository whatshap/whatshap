#ifndef BINOMIAL_H
#define BINOMIAL_H

/**
* Computes the Binomial Coefficient. 
* Use implementation here: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient
*/

int binomial_coefficient(int n, int k);

double binomial_coefficient_large(int n, int k);

double binom_pmf(int n, int k, double p);

# endif // BINOMIAL_H
