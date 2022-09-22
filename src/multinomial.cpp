#include "multinomial.h"
#include "binomial.h"
#include <cmath>
#include <stdint.h>
#include <algorithm>

double multinomial_coefficient(std::vector<uint32_t>& n){
    // sort descending. not necessary, but small optimization for product
    std::vector<uint32_t> s(n);
    std::sort(s.begin(), s.end(), [](const uint32_t a, const uint32_t b) {return a > b; });
    
    // determine sum of samples and create factor list
    uint32_t sum = s[0];
    std::vector<uint32_t> factors;
    for (uint32_t i = 1; i < s.size(); i++) {
        sum += s[i];
        for (uint32_t j = 2; j <= s[i]; j++)
            factors.push_back(j);
    }
    
    double result = 1.0;
    // intertwine mult and div to not create too small or too large numbers
    for (uint32_t i = 0; i < factors.size(); i++) {
        result *= (sum - s[0]); // first s[0] factors are left out
        result /= factors[i];   // factors for s[0] have never been added
    }
	return result;
}

double multinom_pmf(std::vector<uint32_t>& n, std::vector<double>& p) {
    if (n.size() != p.size())
        return 0;       // size of n and p must be identical
    double sum = p[0];
    for (uint32_t i = 1; i < p.size(); i++)
        sum += p[i];
    if (sum != 1.0)
        return 0;       // sum of p must be one for valid input
    if (n.size() == 2)
        return binom_pmf(n[0] + n[1], n[0], p[0]);  // faster for binomial case
    double result = multinomial_coefficient(n);
    for (uint32_t i = 0; i < p.size(); i++)
        result *= p[i];
    return result;
}
