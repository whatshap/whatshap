#ifndef MECHEADER_H
#define MECHEADER_H

#include "readset.h"
#include "pedigree.h"
#include <vector>
#include <stdint.h>

typedef uint32_t Position;
typedef float MecScore;
typedef std::vector<bool> Bipartition;
typedef uint32_t ReadId;
typedef int8_t Allele;
typedef uint16_t RowIndex;
typedef uint32_t Transmission;
typedef std::vector<MecScore> Balance;

static const uint64_t m1  = 0x5555555555555555;
static const uint64_t m2  = 0x3333333333333333;
static const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
static const uint64_t h01 = 0x0101010101010101;
static inline uint64_t popcount(uint64_t bitv) {
    // copied from Wikipedia (https://en.wikipedia.org/wiki/Hamming_weight)
    bitv -= (bitv >> 1) & m1;
    bitv = (bitv & m2) + ((bitv >> 2) & m2);
    bitv = (bitv + (bitv >> 4)) & m4;
    return (bitv * h01) >> 56;
}

static const uint32_t MAX_ROW_LIMIT = 65535;

#endif
