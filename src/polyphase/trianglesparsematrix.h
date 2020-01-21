#ifndef TRIANGLESPARSEMATRIX_H
#define TRIANGLESPARSEMATRIX_H

#include <unordered_map>
#include <set>
#include <vector>

/**
 * A simple storage class, to sparsely store float values for pairs of non-negative integers.
 */
class TriangleSparseMatrix {
 
public:
	TriangleSparseMatrix();
    uint64_t entryToIndex(uint32_t i, uint32_t j);
    uint64_t size();
    float get(uint32_t i, uint32_t j);
    void set(uint32_t i, uint32_t j, float v);
    std::vector<std::pair<uint32_t, uint32_t>> getEntries();

private:    
    std::set<uint64_t> entries;
	std::unordered_map<uint64_t, float> m;
};

#endif

