#ifndef GLOBALS_H
#define GLOBALS_H

#include <limits>
#include <cstdlib>
#include <string>
// #include <x86intrin.h>

/**
 * The time limit in seconds, can be set globally
 * Note: This is currently only respected by the ILP, the heuristic doesn't care about it and the reduction rules don't take it into account either
 */
extern int time_limit;

extern double threshold;

extern int verbosity;
extern int no_threads;
extern double eps;

/**
 * Returns the number of set bits in a 64bit-word.
 */
const uint64_t m1  = 0x5555555555555555;
const uint64_t m2  = 0x3333333333333333;
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;
const uint64_t h01 = 0x0101010101010101;
inline uint64_t popcount(uint64_t bitv) {
    // copied from Wikipedia (https://en.wikipedia.org/wiki/Hamming_weight)
    bitv -= (bitv >> 1) & m1;
    bitv = (bitv & m2) + ((bitv >> 2) & m2);
    bitv = (bitv + (bitv >> 4)) & m4;
    return (bitv * h01) >> 56;
    //return _mm_popcnt_u64(bitv);
}

#endif /* GLOBALS_H */

