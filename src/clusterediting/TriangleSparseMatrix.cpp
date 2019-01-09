#include "TriangleSparseMatrix.h"
#include <cmath>

//using namespace std;

TriangleSparseMatrix::TriangleSparseMatrix() { }

uint64_t TriangleSparseMatrix::entryToIndex(uint32_t i, uint32_t j) {
    if (i < j)
        return entryToIndex(j, i);
    else if (i > j)
        return (i*(i-1)/2)+j;
    else
        return -1;
}

uint64_t TriangleSparseMatrix::size() {
    return m.size();
}

float TriangleSparseMatrix::get(uint32_t i, uint32_t j) {
    uint64_t index = entryToIndex(i, j);
    std::unordered_map<uint64_t, float>::const_iterator it = m.find(index);
    if (it != m.end())
        return it->second;
    else
        return 0;
}

void TriangleSparseMatrix::set(uint32_t i, uint32_t j, float v) {
    uint64_t index = entryToIndex(i, j);
    if (index != -1) {
        m[index] = v;
        entries.insert(index);
    }
}

std::vector<std::pair<uint32_t, uint32_t>> TriangleSparseMatrix::getEntries() {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    for (uint64_t entry : entries) {
        uint64_t u = std::ceil(std::sqrt(2*(entry+1)+0.25) - 0.5);
        uint64_t v = entry - u * (u-1) / 2;
        pairs.push_back(std::pair<uint32_t, uint32_t>((uint32_t)u, (uint32_t)v));
    }
    return pairs;
}

