#include "ReadScoring.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>

void ReadScoring::scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    double hammingDistSame = 0;
    double hammingDistDiff = 0;
    computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, hammingDistSame, hammingDistDiff, minOverlap, ploidy);
//     std::cout<<"Global error rate = "<<hammingDistSame<<std::endl;
//     std::cout<<"Global diff rate = "<<hammingDistDiff<<std::endl;
    
    // estimte error rate for reads in same haplotype
       
//     double avgDisagr = 0.0;
//     size_t numPairs = 0;
//     size_t numBases = 0;
//     double hammingDistSame = errorrate; //2*(1.0-errorrate)*errorrate;
    
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlaps.getEntries();
    
//     for (std::pair<size_t, size_t> p : entries) {
//         size_t i = p.first;
//         size_t j = p.second;
//         numPairs++;
//         numBases += overlaps.get(i, j);
//         avgDisagr += (double)(diffs.get(i, j));
//     }
    
//     double fracSame = ploidy > 1 ? 1.0/ploidy : 0.0;
//     /*double */hammingDistDiff = 1.0;
//     
//     if (numBases == 0) {
//         avgDisagr = 0.0;
//         hammingDistDiff = (1.0+hammingDistSame)/2;
//     } else {
//         avgDisagr = avgDisagr / numBases;
//         /* 
//          * (frac_same*num_pairs)*errorrate + ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr
// 		 * => ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate
// 		 * => x = (num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate) / ((1-frac_same)*num_pairs)
//          */
// 		hammingDistDiff = (numPairs*avgDisagr - (fracSame*numPairs)*hammingDistSame) / ((1.0-fracSame)*numPairs);
// 		hammingDistDiff = std::max(hammingDistSame, std::min((1.0+hammingDistSame)/2, hammingDistDiff));
//     }
    
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

void ReadScoring::scoreReadsetLocal(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy/*, const std::vector<uint32_t> genotypes*/) const {
    uint32_t numReads = readset->size();
    
    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    double hammingDistSame = 0;
    double hammingDistDiff = 0;
    computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, hammingDistSame, hammingDistDiff, minOverlap, ploidy/*, genotypes*/);
    
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
//     std::cout<<"Classic same = "<<defaultSameDist<<std::endl;
//     std::cout<<"Classic diff = "<<defaultDiffDist<<std::endl;
//     std::cout<<"New same = "<<hammingDistSame<<std::endl;
//     std::cout<<"New diff = "<<hammingDistDiff<<std::endl;
//     double defaultSameDist = hammingDistSame;
//     double defaultDiffDist = hammingDistDiff;
    
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
			bool found = false;
            for (uint32_t k = 0; k < positions[coveredReads[idx]].size() && !found; k++) {
                if (positions[coveredReads[idx]][k] >= start) {
                    firstIdx.push_back(k);
                    found = true;
                }
            }
			if (!found)
				firstIdx.push_back(0);
			found = false;
            for (uint32_t k = 0; k < positions[coveredReads[idx]].size() && !found; k++) {
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
                uint32_t ov = 0;
                uint32_t di = 0;
                uint32_t k = firstIdx[idx1];
                uint32_t l = firstIdx[idx2];
                while (k < lastIdx[idx1] && l < lastIdx[idx2]) {
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
        same = std::max(same, 0.01);
		diff = std::min(1.0, std::max(diff, same + 0.01));
        result->set(i, j, logratioSim(ov, di, same, diff));
    }
    
    // cleanup
    delete &snpPositions;
}

void ReadScoring::scoreReadsetPatterns(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy, 
                                       const double errorrate, const uint32_t windowSize) const {    
    // compute overlap and differences for all read pairs in sparse datastrcutures
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    
    // compute length of overlap and difference for all read pairs
    double hammingDistSame = 0;
    double hammingDistDiff = 0;
    computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, hammingDistSame, hammingDistDiff, minOverlap, ploidy);
    
    // create index that maps genome positions to variant positions
    std::unordered_set<uint32_t> positionSet;
    for (std::vector<uint32_t> l : positions) {
        for (uint32_t p : l) {
            positionSet.insert(p);
        }
    }
    
    std::vector<uint32_t> allPositions(positionSet.begin(), positionSet.end());
    std::sort(allPositions.begin(), allPositions.end());
    positionSet.clear();
    std::unordered_map<uint32_t, uint32_t> posIndex;
    for (uint32_t i = 0; i < allPositions.size(); i++) {
        posIndex[allPositions[i]] = i;
    }
    
    // compute pattern counts for each window
    uint32_t numPatterns = std::pow(2, windowSize);
    uint32_t numWindows = (allPositions.size() + windowSize - 1)/windowSize;
    std::vector<std::vector<uint32_t>> patternCount(numWindows, std::vector<uint32_t>(numPatterns, 0));
    std::vector<std::vector<uint32_t>> readsInWindow(numWindows, std::vector<uint32_t>());
    std::vector<std::vector<uint32_t>> patternOfReadInWindow(numWindows, std::vector<uint32_t>());
    std::vector<std::vector<uint32_t>> presentOfReadInWindow(numWindows, std::vector<uint32_t>());
    
    for (uint32_t read = 0; read < positions.size(); read++) {
        // increment pattern count for each window, where this read covers ALL positions
//         std::cout<<"read"<<read<<": ";
        uint32_t lastWindow = 0xffffffff;
        uint32_t curPattern = 0;
        uint32_t present = 0;
        for (uint32_t rpos = 0; rpos < alleles[read].size(); rpos++) {
            uint32_t pos = positions[read][rpos];
            uint32_t varPos = posIndex[pos];
            uint32_t curWindow = varPos / windowSize;
            uint32_t offset = (windowSize - 1) - (varPos % windowSize);
            if (curWindow > lastWindow || lastWindow == 0xffffffff) {
                // new window reached, check for pattern count increase in old window
                if (present == numPatterns - 1) {
                    patternCount[lastWindow][curPattern]++;
                } else if (present > 0) {
                    readsInWindow[curWindow].push_back(read);
                    patternOfReadInWindow[curWindow].push_back(curPattern);
                    presentOfReadInWindow[curWindow].push_back(present);
                    //std::cout<<"Put Window "<<curWindow<<", Read "<<read<<" : "<<curPattern<<", "<<present<<std::endl;
                }
//                 std::cout<<"\t"<<curWindow<<":"<<present<<"/"<<curPattern;
                lastWindow = curWindow;
                present = 0;
                curPattern = 0;
            }
            curPattern |= (alleles[read][rpos] << offset);
            present |= (1U << offset);
        }
//         std::cout<<std::endl;
    }
    
    // infer local haplotype patterns for each window
    std::vector<std::vector<uint32_t>> patternsInWindow;
    std::vector<std::vector<uint32_t>> patternMultiplicityInWindow;
    for (uint32_t window = 0; window < numWindows; window++) {
        /* determine relevant patterns: find pattern with highest count, until <ploidy> many patterns
         * are chosen. To chose the same pattern k times, it must be more than k times more frequent
         * than all other patterns. */
        patternsInWindow.push_back(std::vector<uint32_t>());
        patternMultiplicityInWindow.push_back(std::vector<uint32_t>());
        std::vector<uint32_t> timesChosen(numPatterns, 0);
        for (uint32_t j = 0; j < ploidy; j++) {
            uint32_t maxCount = 0;
            uint32_t maxPattern = 0;
            for (uint32_t pattern = 0; pattern < numPatterns; pattern++) {
                if (patternCount[window][pattern] / (timesChosen[pattern]+1) > maxCount) {
                    maxCount = patternCount[window][pattern] / (timesChosen[pattern]+1);
                    maxPattern = pattern;
                }
            }
            timesChosen[maxPattern]++;
        }
        // write chosen patterns into 2D vector
//         std::cout<<"Window "<<window<<": ";
        for (uint32_t pattern = 0; pattern < numPatterns; pattern++) {
//             if (pattern % 4 == 0)
//                 std::cout<<std::endl;
//             std::cout<<"\t"<<patternCount[window][pattern]<<" ("<<timesChosen[pattern]<<")";
            if (timesChosen[pattern] > 0) {
                patternsInWindow[window].push_back(pattern);
                patternMultiplicityInWindow[window].push_back(timesChosen[pattern]);
            }
        }
//         std::cout<<std::endl;
    }
    
    // iterate over all windows
    for (uint32_t window = 0; window < numWindows; window++) {
        // write patterns as allele vectors
        std::vector<std::vector<uint32_t>> patternAlleles;
        uint32_t numLocalPatterns = patternsInWindow[window].size();
        for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
            std::vector<uint32_t> vec(windowSize, 0);
            for (uint32_t pos = 0; pos < windowSize; pos++) {
                vec[0] = (patternsInWindow[window][pattern] & (1U << (windowSize - 1 - pos))) > 0;
            }
            patternAlleles.push_back(vec);
        }
        
        std::vector<std::vector<double>> prob;
        double e = errorrate > 0.0 && errorrate < 1.0 ? errorrate : 0.05;
        // compute probabilities for each read to originate from each pattern
        for (uint32_t read = 0; read < readsInWindow[window].size(); read++) {
            std::vector<double> probRead(numLocalPatterns, 0.0);
            double sum = 0.0;
            for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
                // bit magic: inner xor detects differences between pattern and read, outer operator ensures zero bits on undefined positions for read
                uint64_t numMatch = popcount(((patternOfReadInWindow[window][read] ^ patternsInWindow[window][pattern]) ^ (numPatterns-1)) & presentOfReadInWindow[window][read]);
                uint64_t numMismatch = popcount((patternOfReadInWindow[window][read] ^ patternsInWindow[window][pattern]) & presentOfReadInWindow[window][read]);
                if (numMatch + numMismatch <= 0 || numMatch + numMismatch > windowSize) {
                    std::cout<<"Invalid match/mismatch count: Window "<<window<<", Read "<<readsInWindow[window][read]<<"("<<read<<"): "<<numMatch<<" and "<<numMismatch<<std::endl;
//                     std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" ^ "<<(numPatterns-1)<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
//                     std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
                    continue;
                }
//                 std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" ^ "<<(numPatterns-1)<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
                double factor = patternMultiplicityInWindow[window][pattern] * std::pow(e, numMismatch) * std::pow(1-e, numMatch);
                probRead[pattern] = factor;
                sum += factor;
                if (factor == 0.0) {
                    std::cout<<"Factor was zero for "<<patternOfReadInWindow[window][read]<<" ("<<presentOfReadInWindow[window][read]<<") and "<<patternsInWindow[window][pattern]<<std::endl;
                }
            }
            for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
                probRead[pattern] /= sum;
            }
            prob.push_back(probRead);
        }
        
        // compute scoring for all read pairs in current window
        for (uint32_t read1 = 0; read1 < readsInWindow[window].size(); read1++) {
            for (uint32_t read2 = read1 + 1; read2 < readsInWindow[window].size(); read2++) {
                uint32_t readId1 = readsInWindow[window][read1];
                uint32_t readId2 = readsInWindow[window][read2];
                if (overlaps.get(readId1, readId2) < minOverlap)
                    continue;
                double pSame = 0.0;
                double pDiff = 0.0;
                for (uint32_t p1 = 0; p1 < numLocalPatterns; p1++) {
                    for (uint32_t p2 = 0; p2 < numLocalPatterns; p2++) {
                        if (p1 == p2) {
                            pSame += prob[read1][p1] * prob[read2][p2];
                        } else {
                            pDiff += prob[read1][p1] * prob[read2][p2];
                        }
                    }
                }
                if (std::abs(1.0 - (pSame + pDiff)) > 0.0001) {
                    std::cout<<"pSame + pDiff = "<<(pDiff+pSame)<<" in window "<<window<<std::endl;
                } else {
                    if (pSame == 0)
                        result->set(readId1, readId2, -std::numeric_limits<float>::infinity());
                    else if (pDiff == 0) {
                        result->set(readId1, readId2, std::numeric_limits<float>::infinity());
//                         std::cout<<window<<" : "<<read1<<"-"<<read2<<" : "<<readId1<<"-"<<readId2<<" : "<<patternOfReadInWindow[window][read1]<<"-"<<patternOfReadInWindow[window][read1]<<" : "<<presentOfReadInWindow[window][read1]<<"-"<<presentOfReadInWindow[window][read1]<<std::endl;
                    }
                    else
                        result->set(readId1, readId2, result->get(readId1, readId2) + log(pSame / pDiff));
                }
            }
        }
    }
}

void ReadScoring::computeStartEndOverlapDiff (const ReadSet* readset, std::vector<uint32_t>& begins, std::vector<uint32_t>& ends, std::vector<std::vector<uint32_t>>& positions, std::vector<std::vector<uint32_t>>& alleles, TriangleSparseMatrix& overlaps, TriangleSparseMatrix& diffs, double& distSame, double& distDiff, const uint32_t minOverlap, const uint32_t ploidy/*, const std::vector<uint32_t> genotypes*/) const {
    // copy all relevant information from the readset into vectors for efficient access
    uint32_t numReads = readset->size();
    std::unordered_set<uint32_t> allPos;
    std::unordered_map<uint32_t, uint32_t> posMap;
    for (uint32_t i = 0; i < numReads; i++) {
        begins.push_back(readset->get(i)->firstPosition());
        ends.push_back(readset->get(i)->lastPosition());
        std::vector<uint32_t> pos;
        std::vector<uint32_t> all;
        for (int k = 0; k < readset->get(i)->getVariantCount(); k++) {
            pos.push_back(readset->get(i)->getPosition(k));
            all.push_back(readset->get(i)->getAllele(k));
            allPos.insert(readset->get(i)->getPosition(k));
        }
        positions.push_back(pos);
        alleles.push_back(all);
    }
    
    // create position map
    std::vector<uint32_t> posList(allPos.begin(), allPos.end());
    std::sort(posList.begin(), posList.end());
    for (uint32_t i = 0; i < posList.size(); i++) {
        posMap[posList[i]] = i;
    }
    
    // determine length of the longest read
    uint32_t longestReadLen = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        longestReadLen = std::max(longestReadLen, ends[i] - begins[i]);
    }
    
    // iterate over all read pairs (efficiently omitting those who can certainly not overlap)
    std::vector<double> relativeDiffs;
    for (uint32_t i = 0; i < numReads; i++) {
        // iterate until start position of read is behind required start
        for (uint32_t j = i+1; begins[j] <= ends[i] && j < numReads; j++) {
            if (ends[i] < begins[j] || ends[j] < begins[i])
                continue;
            uint32_t ov = 0;
            uint32_t di = 0;
            uint32_t k = 0;
            uint32_t l = 0;
            while (k < positions[i].size() && l < positions[j].size()) {
                if (positions[i][k] == positions[j][l]) {
//                     if (alleles[i][k] != alleles[j][l]) {
//                         di+=3; ov+=3;
//                     } else if (!genotypes.empty()){
//                         if (genotypes[posMap[positions[i][k]]] == 1) {
//                             if (alleles[i][k] == 0) {
//                                 ov+=2;
//                             } else {
//                                 ov+=6;
//                             }
//                         }
//                         if (genotypes[posMap[positions[i][k]]] == 2) {
//                             if (alleles[i][k] == 0) {
//                                 ov+=3;
//                             } else {
//                                 ov+=3;
//                             }
//                         }
//                         if (genotypes[posMap[positions[i][k]]] == 3) {
//                             if (alleles[i][k] == 0) {
//                                 ov+=6;
//                             } else {
//                                 ov+=2;
//                             }
//                         }
//                     } else {
//                         ov+=3;
//                     }
                    if (alleles[i][k] != alleles[j][l]) {
                        di++;
                    }
                    ov++; k++; l++; 
                } else if (positions[i][k] < positions[j][l]) {
                    k++;
                } else {
                    l++;
                }
            }
//             ov /= 3;
//             di /= 3;
            if (ov >= minOverlap) {
                overlaps.set(i, j, ov);
                diffs.set(i, j, di);
                relativeDiffs.push_back((double)di / (double)ov);
            }
        }
    }
    
    // estimate error rate of reads within equal haplotype
    computeCutoff(numReads, ploidy, relativeDiffs, distSame, distDiff);
}

void ReadScoring::computeCutoff(const uint32_t numReads, const uint32_t ploidy, std::vector<double> relDiffs, double& distSame, double& distDiff) const {
    std::sort(relDiffs.begin(), relDiffs.end());
    double sameSum = 0.0;
    int sameNum = 0;
    double diffSum = 0.0;
    int diffNum = 0;
    
    double p = (double)ploidy;
    double n = (double)numReads;
    uint32_t cutoff = std::min((uint32_t)1, (uint32_t)(relDiffs.size() - 1));
    if (ploidy < numReads) {
        cutoff = (uint32_t) std::ceil((p * (n/p) * (n/p-1) / 2) / 
                    ( (p * (n/p) * (n/p-1) / 2) + (p*(p-1)/2)*(n/p)*(n/p) ) * relDiffs.size());
        for (unsigned int i = 0; i < relDiffs.size(); i++) {
            if (i < cutoff) {
                sameSum += relDiffs[i];
                sameNum++;
            } else {
                diffSum += relDiffs[i];
                diffNum++;
            }
        }
        
        if (cutoff == 0) {
            distSame = 0.1;
        } else {
            distSame = sameSum/sameNum;
        }
        distDiff = diffSum/diffNum;
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
