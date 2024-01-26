#include "multinomial.h"
#include "binomial.h"
#include <cmath>
#include <stdint.h>
#include <algorithm>
#include <limits>

double log_multinomial_coefficient(std::vector<uint32_t>& n){
    // sort descending. not necessary, but small optimization for product
    std::vector<uint32_t> s(n.begin(), n.end());
    std::sort(s.begin(), s.end(), [](const uint32_t a, const uint32_t b) {return a > b; });
    
    // determine sum of samples and create factor list
    uint32_t sum = s[0];
    std::vector<uint32_t> factors;
    for (uint32_t i = 1; i < s.size(); i++) {
        sum += s[i];
        for (uint32_t j = 2; j <= s[i]; j++)
            factors.push_back(j);
    }
    
    double result = 0.0;
    double buffer = 1.0;
    // intertwine mult and div to avoid over/underflows
    for (uint32_t i = 0; i < factors.size(); i++) {
        double addition = (double)(sum - s[0]) / (double)factors[i];
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

bool check_multinom_input(std::vector<uint32_t>& n, std::vector<double>& p) {
    if (n.size() != p.size())
        return false;       // size of n and p must be identical
    double sum = p[0];
    for (uint32_t i = 1; i < p.size(); i++)
        sum += p[i];
    if (sum != 1.0)
        return false;       // sum of p must be one for valid input
    return true;
}

double multinom_pmf(std::vector<uint32_t>& n, std::vector<double>& p) {
    if (n.size() == 2)
        return binom_pmf(n[0] + n[1], n[0], p[0]);  // faster for binomial case
    if (!check_multinom_input(n, p))
        return 0;
    double result = std::exp(log_multinomial_coefficient(n));
    for (uint32_t i = 0; i < p.size(); i++)
        result *= p[i];
    return result;
}

double log_multinom_pmf(std::vector<uint32_t>& n, std::vector<double>& p) {
    if (n.size() == 2)
        return log_binom_pmf(n[0] + n[1], n[0], p[0]);  // faster for binomial case
    if (!check_multinom_input(n, p))
        return 0;
    double result = log_multinomial_coefficient(n);
    for (uint32_t i = 0; i < p.size(); i++)
        result += std::log(p[i]);
    return result;
}
