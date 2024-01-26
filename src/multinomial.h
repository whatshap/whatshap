#ifndef MULTINOMIAL_H
#define MULTINOMIAL_H

#include <vector>
#include <stdint.h>

/**
* Computes the Multinomial Coefficient. 
* Formula taken from https://en.wikipedia.org/wiki/Multinomial_distribution
*/

double log_multinomial_coefficient(std::vector<uint32_t>& n);

bool check_multinom_input(std::vector<uint32_t>& n, std::vector<double>& p);

double multinom_pmf(std::vector<uint32_t>& n, std::vector<double>& p);

double log_multinom_pmf(std::vector<uint32_t>& n, std::vector<double>& p);

# endif // MULTINOMIAL_H
