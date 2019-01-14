#include "ReadScoring.h"
#include <vector>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>

void ReadScoring::scoreReadset(TriangleSparseMatrix *result, ReadSet *readset, const double errorrate, const uint32_t minOverlap, const uint32_t ploidy) const {

    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint64_t> begins;
    std::vector<uint64_t> ends;
    std::vector<std::vector<uint64_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, minOverlap);
       
    double avgDisagr = 0.0;
    size_t numPairs = 0;
    size_t numBases = 0;
    double hammingDistSame = 2*(1.0-errorrate)*errorrate;
    
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlaps.getEntries();
    
    for (std::pair<size_t, size_t> p : entries) {
        size_t i = p.first;
        size_t j = p.second;
        numPairs++;
        numBases += overlaps.get(i, j);
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
    
    std::unordered_map<uint64_t, float> cache;
    for (std::pair<uint32_t, uint32_t> p : entries) {
        uint32_t i = p.first;
        uint32_t j = p.second;
        uint32_t ov = overlaps.get(i, j);
        uint32_t di = diffs.get(i, j);
        uint64_t ovdi = (ov*(ov+1))/2+di;
        std::unordered_map<uint64_t, float>::const_iterator it = cache.find(ovdi);
        if (it == cache.end()) {
            cache[ovdi] = logratioSim(ov, di, hammingDistSame, hammingDistDiff);
        }
        result->set(i, j, cache[ovdi]);
    }
    
}

void ReadScoring::scoreReadset(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    uint32_t numReads = readset->size();
    
    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint64_t> begins;
    std::vector<uint64_t> ends;
    std::vector<std::vector<uint64_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, minOverlap);
    
    // determine default relative hamming distance for same and different haplotypes
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlaps.getEntries();
    double avgDisagr = 0.0;
    size_t numPairs = 0;
    size_t numBases = 0;
    double defaultSameDist = 2*(1.0-0.03)*0.03;
    for (std::pair<size_t, size_t> p : entries) {
        size_t i = p.first;
        size_t j = p.second;
        numPairs++;
        numBases += overlaps.get(i, j);
        avgDisagr += (double)(diffs.get(i, j));
    }
    
    double fracSame = ploidy > 1 ? 1.0/ploidy : 0.0;
    double defaultDiffDist = 1.0;
    
    if (numBases == 0) {
        avgDisagr = 0.0;
        defaultDiffDist = (1.0+defaultSameDist)/2;
    } else {
        avgDisagr = avgDisagr / numBases;
        /* 
         * (frac_same*num_pairs)*errorrate + ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr
		 * => ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate
		 * => x = (num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate) / ((1-frac_same)*num_pairs)
         */
		defaultDiffDist = (numPairs*avgDisagr - (fracSame*numPairs)*defaultSameDist) / ((1.0-fracSame)*numPairs);
		defaultDiffDist = std::max(defaultSameDist, std::min((1.0+defaultSameDist)/2, defaultDiffDist));
    }
    
    // compute longest read length and average read length (in base pairs) and divide by 2
    uint32_t windowSize = 0;
    uint32_t longestReadLen = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        uint32_t len = readset->get(i)->lastPosition() - readset->get(i)->firstPosition();
        windowSize += len;
        longestReadLen = std::max(longestReadLen, len);
    }
    windowSize /= (2*numReads);

    //divide snp positions by window size
    std::vector<uint32_t> windowStarts;
    std::vector<uint32_t>& snpPositions = *(readset->get_positions());
    std::unordered_map<uint32_t, double> posToSameDist;
    std::unordered_map<uint32_t, double> posToDiffDist;
    
    uint32_t windowStartPosition = 0;
    for (uint32_t current = 0; current < snpPositions.size(); current++) {
        if (snpPositions[current] - windowStartPosition > windowSize || current == 0) {
            windowStarts.push_back(current);
            windowStartPosition = snpPositions[current];
//             std::cout<<"New SNP Window starts at "<<current<<" / "<<snpPositions[current]<<std::endl;
        }
    }
    windowStarts.push_back(snpPositions.size()+1);

    // determine relative hamming distance for same and different haplotypes for each window
    std::vector<double> sameDist;
    std::vector<double> diffDist;
    for (uint32_t windowIdx = 0; windowIdx < windowStarts.size()-1; windowIdx++) {
        // window bounds
        uint32_t start = snpPositions[windowStarts[windowIdx]];
        uint32_t end = snpPositions[windowStarts[windowIdx+1]-1];
        std::vector<uint32_t> coveredReads;
        std::vector<double> dists;
        
        // binary search to find first read, which can theoretically cover this window
        uint32_t firstIndex = std::lower_bound(begins.begin(), begins.end(), start-longestReadLen) - begins.begin();
        // iterate until start position of read is behind required start
        for (uint32_t j = firstIndex; begins[j] <= start && j < begins.size(); j++) {
            if (ends[j] >= end) {
                coveredReads.push_back(j);
            }
        }
        
        // for all candidate reads, determine their first and last defined position, which lies inside this window
        std::vector<uint32_t> firstIdx;
        std::vector<uint32_t> lastIdx;
        for (uint32_t idx = 0; idx < coveredReads.size(); idx++) {
			//if (idx >= coveredReads.size())
			//	continue;
			//if (coveredReads[idx] >= positions.size())
			//	continue;
			bool found = false;
            for (uint32_t k = 0; k < positions[coveredReads[idx]].size() && !found; k++) {
				//if (k >= positions[coveredReads[idx]].size())
				//	continue;
                if (positions[coveredReads[idx]][k] >= start) {
                    firstIdx.push_back(k);
                    found = true;
                }
            }
			if (!found)
				firstIdx.push_back(0);
			found = false;
            for (uint32_t k = 0; k < positions[coveredReads[idx]].size() && !found; k++) {
				//if (k >= positions[coveredReads[idx]].size())
				//	continue;
                if (positions[coveredReads[idx]][k] > end) {
                    lastIdx.push_back(k);
                    found = true;
                }
            }
			if (!found)
				lastIdx.push_back(positions[coveredReads[idx]].size());
        }
        
        // iterate over all read pairs and determine their relative hamming distance
        for (uint32_t idx1 = 0; idx1 < coveredReads.size(); idx1++) {
            for (uint32_t idx2 = idx1 + 1; idx2 < coveredReads.size(); idx2++) {
				//if (idx1 >= firstIdx.size())
				//	continue;
				//if (idx2 >= firstIdx.size())
				//	continue;
				//if (idx1 >= lastIdx.size())
				//	continue;
				//if (idx2 >= lastIdx.size())
				//	continue;
				//if (coveredReads[idx1] >= positions.size())
				//	continue;
				//if (coveredReads[idx2] >= positions.size())
				//	continue;
                uint32_t ov = 0;
                uint32_t di = 0;
                uint32_t k = firstIdx[idx1];
                uint32_t l = firstIdx[idx2];
                while (k < lastIdx[idx1] && l < lastIdx[idx2] /*&& k < positions[coveredReads[idx1]].size() && l < positions[coveredReads[idx2]].size()*/) {
                    if (positions[coveredReads[idx1]][k] == positions[coveredReads[idx2]][l]) {
                        if (alleles[coveredReads[idx1]][k] != alleles[coveredReads[idx2]][l])
                            di++;
                        ov++; k++; l++;
                    } else if (positions[coveredReads[idx1]][k] < positions[coveredReads[idx2]][l]) {
                        k++;
                    } else {
                        l++;
                    }
                }
                if (ov >= minOverlap) {
                    dists.push_back(((double)di)/((double)ov));
                }
            }
        }
        
        // take average of first 1/ploidy fraction as sameDist value and average of rest as diffDist value
        std::sort(dists.begin(), dists.end());
        
        double same = 0;
        double diff = 0;
        if (dists.size() / ploidy == 0) {
            // fallback: if not enough read pairs, use default values
            same = defaultSameDist;
            diff = defaultDiffDist;
        } else {
            for (uint32_t j = 0; j < dists.size(); j++) {
                if (j < dists.size() / ploidy)
                    same += dists[j];
                else
                    diff += dists[j];
            }
            same /= std::floor((double)(dists.size()) / ploidy);
            diff /= (dists.size() - std::floor((double)(dists.size()) / ploidy));
        }
        sameDist.push_back(same);
        sameDist.push_back(diff);
        
        // store values in a map for every snp position
        for (uint32_t j = windowStarts[windowIdx]; j < windowStarts[windowIdx+1]; j++) {
            posToSameDist[snpPositions[j]] = same;
            posToDiffDist[snpPositions[j]] = diff;
        }
    }
            
    // now, iterate over all overlapping read pairs and compute their score
    for (std::pair<uint32_t, uint32_t> p : entries) {
        uint32_t i = p.first;
        uint32_t j = p.second;
        uint32_t ov = overlaps.get(i, j);
        uint32_t di = diffs.get(i, j);
        
        // zig zag over read positions to determine their individual same/diff dists
        double same = 0.0;
        double diff = 0.0;
        uint32_t k = 0;
        uint32_t l = 0;
        while (k < positions[i].size() && l < positions[j].size()) {
            if (positions[i][k] == positions[j][l]) {
                same += posToSameDist[positions[i][k]];
                diff += posToDiffDist[positions[i][k]];
                k++; l++;
            } else if (positions[i][k] < positions[j][l]) {
                k++;
            } else {
                l++;
            }
        }
        same /= ov;
        diff /= ov;
		//diff = (diff+2*same)/3;
		same = std::max(same, 0.01);
		diff = std::min(1.0, std::max(diff, same + 0.01));
        result->set(i, j, logratioSim(ov, di, same, diff));
    }
    
    // cleanup
    delete &snpPositions;
}

void ReadScoring::computeStartEndOverlapDiff(const ReadSet* readset, std::vector<uint64_t>& begins, std::vector<uint64_t>& ends, std::vector<std::vector<uint64_t>> &positions, std::vector<std::vector<uint32_t>> &alleles, TriangleSparseMatrix& overlaps, TriangleSparseMatrix& diffs, const uint32_t minOverlap) const {
    uint32_t numReads = readset->size();
    for (uint32_t i = 0; i < numReads; i++) {
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
    
    for (uint32_t i = 0; i < numReads; i++) {
        for (uint32_t j = i+1; j < numReads; j++) {
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
                overlaps.set(i, j, ov);
                diffs.set(i, j, di);
            }
        }
    }
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
