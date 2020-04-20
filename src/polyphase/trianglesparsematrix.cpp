#include "trianglesparsematrix.h"
#include <cmath>
#include <algorithm>

//using namespace std;

TriangleSparseMatrix::TriangleSparseMatrix() : maxDim(0) {  }

uint64_t TriangleSparseMatrix::entryToIndex(uint32_t i, uint32_t j) {
    if (i < j)
        return entryToIndex(j, i);
    else if (i > j)
        return ((uint64_t)i*(uint64_t)(i-1)/2)+(uint64_t)j+1;
    else
        return 0;
}

uint64_t TriangleSparseMatrix::size() {
    return m.size();
}

uint32_t TriangleSparseMatrix::getMaxDim() {
    return maxDim;
}

float TriangleSparseMatrix::get(uint32_t i, uint32_t j) {
    uint64_t index = entryToIndex(i, j);
    std::unordered_map<uint64_t, MatrixItem>::const_iterator it = m.find(index);
    if (it != m.end())
        return it->second.v;
    else
        return 0;
}

TriangleSparseMatrix::DoubleInt TriangleSparseMatrix::getDoubleInt(uint32_t i, uint32_t j) {
    uint64_t index = entryToIndex(i, j);
    std::unordered_map<uint64_t, MatrixItem>::const_iterator it = m.find(index);
    if (it != m.end())
        return DoubleInt(it->second.i);
    else
        return DoubleInt(0, 0);
}

void TriangleSparseMatrix::set(uint32_t i, uint32_t j, float v) {
    uint64_t index = entryToIndex(i, j);
    if (index != 0) {
        m[index].v = v;
        maxDim = std::max(maxDim, i+1);
        maxDim = std::max(maxDim, j+1);
    }
}

void TriangleSparseMatrix::setDoubleInt(uint32_t i, uint32_t j, uint16_t u1, uint16_t u2) {
    uint64_t index = entryToIndex(i, j);
    if (index != 0) {
        m[index].i = ((uint32_t)u1 << 16) + (uint32_t)u2;
        maxDim = std::max(maxDim, i+1);
        maxDim = std::max(maxDim, j+1);
    }
}

std::vector<uint64_t> TriangleSparseMatrix::getIndices() {
    std::vector<uint64_t> indices;
    for (auto entry : m) {
        indices.push_back(entry.first-1);
    }
    std::sort(indices.begin(), indices.end());
    return indices;
}

std::vector<std::pair<uint32_t, uint32_t>> TriangleSparseMatrix::getEntries() {
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    for (auto entry : m) {
        uint64_t u = std::ceil(std::sqrt(2*(entry.first)+0.25) - 0.5);
        uint64_t v = (entry.first-1) - u * (u-1) / 2;
        pairs.push_back(std::pair<uint32_t, uint32_t>((uint32_t)u, (uint32_t)v));
    }
    return pairs;
}

