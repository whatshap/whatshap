#ifndef TRIANGLESPARSEMATRIX_H
#define TRIANGLESPARSEMATRIX_H

#include <unordered_map>
#include <set>
#include <vector>
#include <cmath>
#include <cstdint>

/**
 * A simple storage class, to sparsely store float values for pairs of non-negative integers.
 */
class TriangleSparseMatrix {
 
public:
    struct DoubleInt {
        uint16_t u1;
        uint16_t u2;
        DoubleInt() : u1(0), u2(0) {}
        DoubleInt(uint16_t v1, uint16_t v2) : u1(v1), u2(v2) {}
        DoubleInt(uint32_t u) : u1(u/65536), u2(u%65536) {}
    };
	TriangleSparseMatrix();
    uint64_t entryToIndex(uint32_t i, uint32_t j);
    uint64_t size();
    uint32_t getMaxDim();
    float get(uint32_t i, uint32_t j);
    DoubleInt getDoubleInt(uint32_t i, uint32_t j);
    void set(uint32_t i, uint32_t j, float v);
    void setDoubleInt(uint32_t i, uint32_t j, uint16_t u1, uint16_t u2);
    std::vector<uint64_t> getIndices();
    std::vector<std::pair<uint32_t, uint32_t>> getEntries();

private:
    union MatrixItem {
        float v;
        uint32_t i;
    };
	std::unordered_map<uint64_t, MatrixItem> m;
    uint32_t maxDim;
};

#endif

