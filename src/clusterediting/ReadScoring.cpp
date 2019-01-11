#include "ReadScoring.h"
#include <vector>
#include <unordered_map>
#include <limits>
#include <cmath>

void ReadScoring::scoreReadset(TriangleSparseMatrix *result, ReadSet *readset, const double errorrate, const uint32_t minOverlap, const uint32_t ploidy) const {
	uint32_t num_reads = readset->size();
    
    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlap;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint64_t> begins;
    std::vector<uint64_t> ends;
    std::vector<std::vector<uint64_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    for (uint32_t i = 0; i < num_reads; i++) {
        begins.push_back(readset->get(i)->firstPosition());
        ends.push_back(readset->get(i)->lastPosition());
        std::vector<uint64_t> pos;
        std::vector<uint32_t> all;
        for (int k = 0; k < readset->get(i)->getVariantCount(); k++) {
            pos.push_back(readset->get(i)->getPosition(k));
            all.push_back(readset->get(i)->getAllele(k));
        }
        positions.push_back(pos);
        alleles.push_back(all);
    }
    
    for (uint32_t i = 0; i < num_reads; i++) {
        for (uint32_t j = i+1; j < num_reads; j++) {
            if (ends[i] < begins[j] || ends[j] < begins[i])
                continue;
            uint32_t ov = 0;
            uint32_t di = 0;
            uint32_t k = 0;
            uint32_t l = 0;
            while (k < positions[i].size() && l < positions[j].size()) {
                if (positions[i][k] == positions[j][l]) {
                    if (alleles[i][k] != alleles[j][l])
                        di++;
                    ov++; k++; l++;
                } else if (positions[i][k] < positions[j][l]) {
                    k++;
                } else {
                    l++;
                }
            }
            if (ov >= minOverlap) {
                overlap.set(i, j, ov);
                diffs.set(i, j, di);
            }
        }
    }
    
    double avgDisagr = 0.0;
    size_t numPairs = 0;
    size_t numBases = 0;
    double hammingDistSame = 2*(1.0-errorrate)*errorrate;
    
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlap.getEntries();
    
    for (std::pair<size_t, size_t> p : entries) {
        size_t i = p.first;
        size_t j = p.second;
        numPairs++;
        numBases += overlap.get(i, j);
        avgDisagr += (double)(diffs.get(i, j));
    }
    
    double fracSame = ploidy > 1 ? 1.0/ploidy : 0.0;
    double hammingDistDiff = 1.0;
    
    if (numBases == 0) {
        avgDisagr = 0.0;
        hammingDistDiff = (1.0+hammingDistSame)/2;
    } else {
        avgDisagr = avgDisagr / numBases;
        /* 
         * (frac_same*num_pairs)*errorrate + ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr
		 * => ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate
		 * => x = (num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate) / ((1-frac_same)*num_pairs)
         */
		hammingDistDiff = (numPairs*avgDisagr - (fracSame*numPairs)*hammingDistSame) / ((1.0-fracSame)*numPairs);
		hammingDistDiff = std::max(hammingDistSame, std::min((1.0+hammingDistSame)/2, hammingDistDiff));
    }
    
    // calculate the actual scores
    
//     TriangleSparseMatrix sim;
//     std::unordered_map<std::pair<size_t, size_t>, float> cache;
//     for (std::pair<size_t, size_t> p : entries) {
//         size_t i = p.first;
//         size_t j = p.second;
//         size_t ov = overlap.get(i, j);
//         size_t di = diffs.get(i, j);
//         std::pair<size_t, size_t> ovdi(ov, di);
//         std::unordered_map<std::pair<size_t, size_t>, float>::const_iterator it = cache.find(ovdi);
//         if (it == cache.end()) {
//             cache[ovdi] = logratioSim(ov, di, hammingDistSame, hammingDistDiff);
//         }
//         sim.set(i, j, cache[ovdi]);
//     }
    
    std::unordered_map<uint64_t, float> cache;
    for (std::pair<uint32_t, uint32_t> p : entries) {
        uint32_t i = p.first;
        uint32_t j = p.second;
        uint32_t ov = overlap.get(i, j);
        uint32_t di = diffs.get(i, j);
        uint64_t ovdi = (ov*(ov+1))/2+di;
        std::unordered_map<uint64_t, float>::const_iterator it = cache.find(ovdi);
        if (it == cache.end()) {
            cache[ovdi] = logratioSim(ov, di, hammingDistSame, hammingDistDiff);
        }
//         sim.set(i, j, cache[ovdi]);
        result->set(i, j, cache[ovdi]);
    }
    
//     return sim;
}

float ReadScoring::logratioSim(const uint32_t overlap, const uint32_t diff, const double distSame, const double distDiff) const {
    double pSame = binomPmf(overlap, diff, distSame);
    double pDiff = binomPmf(overlap, diff, distDiff);
    if (pSame == 0)
        return -std::numeric_limits<float>::infinity();
    else if (pDiff == 0)
        return std::numeric_limits<float>::infinity();
    else
        return log(pSame / pDiff);
}

double ReadScoring::binomPmf(const uint32_t n, const uint32_t k, const double p) const {
    double coeff = 1.0;
    for (uint32_t i = 0; i < k; i++) {
        coeff *= ((n-i) / (k-i));
    }
    return coeff * pow(p, k) * pow(1-p, n-k);
}
