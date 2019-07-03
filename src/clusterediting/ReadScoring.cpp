#include "ReadScoring.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>

void ReadScoring::scoreReadsetGlobal(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    std::vector<uint32_t> posList;
    std::unordered_map<uint32_t, uint32_t> posMap;
    uint32_t longestReadSpan = 0;
    computeStartEnd(readset, begins, ends, positions, alleles, posList, posMap, longestReadSpan);
    
    // compute length of overlap and difference for all read pairs
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    double hammingDistSame = 0;
    double hammingDistDiff = 0;
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, overlaps, diffs, hammingDistSame, hammingDistDiff, minOverlap, ploidy, longestReadSpan);
//     std::cout<<"Global error rate = "<<hammingDistSame<<std::endl;
//     std::cout<<"Global diff rate = "<<hammingDistDiff<<std::endl;
    
    // compute pair wise scores
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlaps.getEntries();
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

void ReadScoring::scoreReadsetLocal(TriangleSparseMatrix* result, ReadSet* readset, const uint32_t minOverlap, const uint32_t ploidy) const {
    uint32_t numReads = readset->size();
    
    // copy relevant information from readset for fast access
    std::vector<uint32_t> begins;
    std::vector<uint32_t> ends;
    std::vector<std::vector<uint32_t>> positions;
	std::vector<std::vector<uint32_t>> alleles;
    std::vector<uint32_t> posList;
    std::unordered_map<uint32_t, uint32_t> posMap;
    uint32_t longestReadSpan = 0;
    computeStartEnd(readset, begins, ends, positions, alleles, posList, posMap, longestReadSpan);
    
    // compute length of overlap and difference for all read pairs
    TriangleSparseMatrix overlaps;
    TriangleSparseMatrix diffs;
    double defaultSameDist = 0;
    double defaultDiffDist = 0;
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, overlaps, diffs, defaultSameDist, defaultDiffDist, minOverlap, ploidy, longestReadSpan);
    
    std::vector<std::pair<uint32_t, uint32_t>> entries = overlaps.getEntries();
    //std::cout<<"defaultSameDist="<<defaultSameDist<<" defaultDiffDist="<<defaultDiffDist<<std::endl;
    //defaultSameDist = std::max(0.05, defaultSameDist);
    //defaultDiffDist = std::min(0.5, std::max(defaultDiffDist, defaultSameDist+0.15));
    
    // compute longest read length and average read length (in base pairs) and divide by 2
    uint32_t windowSize = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        uint32_t len = readset->get(i)->lastPosition() - readset->get(i)->firstPosition();
        windowSize += len;
    }
    windowSize /= (4*numReads);
    
    //divide snp positions by window size
    std::vector<uint32_t> windowStarts;
    std::unordered_map<uint32_t, double> posToSameDist;
    std::unordered_map<uint32_t, double> posToDiffDist;
    
    uint32_t windowStartPosition = 0;
    for (uint32_t current = 0; current < posList.size(); current++) {
        if (posList[current] - windowStartPosition > windowSize || current == 0) {
            windowStarts.push_back(current);
            windowStartPosition = posList[current];
        }
    }
    windowStarts.push_back(posList.size()+1);
    
    // determine relative hamming distance for same and different haplotypes for each window
    for (uint32_t windowIdx = 0; windowIdx < windowStarts.size()-1; windowIdx++) {
        // window bounds
        uint32_t start = posList[windowStarts[windowIdx]];
        uint32_t end = posList[windowStarts[windowIdx+1]-1];

        // compute length of overlap and difference for all read pairs
        TriangleSparseMatrix overlapsLocal;
        TriangleSparseMatrix diffsLocal;
        double localSameDist = 0;
        double localDiffDist = 0;
        computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, overlapsLocal, diffsLocal, 
                           localSameDist, localDiffDist, minOverlap, ploidy, longestReadSpan, start, end);
        
        //localSameDist += 0.10;
        //localDiffDist = std::min(0.5, localSameDist+0.15);
        
        if (diffsLocal.getEntries().size() < ploidy) {
            localSameDist = defaultSameDist;
            localDiffDist = defaultDiffDist;
        }
        
        //std::cout<<"localSameDist="<<localSameDist<<" localDiffDist="<<localDiffDist<<std::endl;
        
        // store values in a map for every snp position
        for (uint32_t j = windowStarts[windowIdx]; j < windowStarts[windowIdx+1]; j++) {
            posToSameDist[posList[j]] = localSameDist;
            posToDiffDist[posList[j]] = localDiffDist;
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
}

void ReadScoring::scoreReadsetPatterns(TriangleSparseMatrix *result, ReadSet *readset, const uint32_t minOverlap, const uint32_t ploidy, 
                                       const double errorrate, const uint32_t windowSize) const {    
//     // compute overlap and differences for all read pairs in sparse datastrcutures
//     TriangleSparseMatrix overlaps;
//     TriangleSparseMatrix diffs;
//     
//     // copy relevant information from readset for fast access
//     std::vector<uint32_t> begins;
//     std::vector<uint32_t> ends;
//     std::vector<std::vector<uint32_t>> positions;
// 	std::vector<std::vector<uint32_t>> alleles;
//     
//     // compute length of overlap and difference for all read pairs
//     double hammingDistSame = 0;
//     double hammingDistDiff = 0;
//     computeStartEndOverlapDiff(readset, begins, ends, positions, alleles, overlaps, diffs, hammingDistSame, hammingDistDiff, minOverlap, ploidy);
//     
//     // create index that maps genome positions to variant positions
//     std::unordered_set<uint32_t> positionSet;
//     for (std::vector<uint32_t> l : positions) {
//         for (uint32_t p : l) {
//             positionSet.insert(p);
//         }
//     }
//     
//     std::vector<uint32_t> allPositions(positionSet.begin(), positionSet.end());
//     std::sort(allPositions.begin(), allPositions.end());
//     positionSet.clear();
//     std::unordered_map<uint32_t, uint32_t> posIndex;
//     for (uint32_t i = 0; i < allPositions.size(); i++) {
//         posIndex[allPositions[i]] = i;
//     }
//     
//     // compute pattern counts for each window
//     uint32_t numPatterns = std::pow(2, windowSize);
//     uint32_t numWindows = (allPositions.size() + windowSize - 1)/windowSize;
//     std::vector<std::vector<uint32_t>> patternCount(numWindows, std::vector<uint32_t>(numPatterns, 0));
//     std::vector<std::vector<uint32_t>> readsInWindow(numWindows, std::vector<uint32_t>());
//     std::vector<std::vector<uint32_t>> patternOfReadInWindow(numWindows, std::vector<uint32_t>());
//     std::vector<std::vector<uint32_t>> presentOfReadInWindow(numWindows, std::vector<uint32_t>());
//     
//     for (uint32_t read = 0; read < positions.size(); read++) {
//         // increment pattern count for each window, where this read covers ALL positions
// //         std::cout<<"read"<<read<<": ";
//         uint32_t lastWindow = 0xffffffff;
//         uint32_t curPattern = 0;
//         uint32_t present = 0;
//         for (uint32_t rpos = 0; rpos < alleles[read].size(); rpos++) {
//             uint32_t pos = positions[read][rpos];
//             uint32_t varPos = posIndex[pos];
//             uint32_t curWindow = varPos / windowSize;
//             uint32_t offset = (windowSize - 1) - (varPos % windowSize);
//             if (curWindow > lastWindow || lastWindow == 0xffffffff) {
//                 // new window reached, check for pattern count increase in old window
//                 if (present == numPatterns - 1) {
//                     patternCount[lastWindow][curPattern]++;
//                 } else if (present > 0) {
//                     readsInWindow[curWindow].push_back(read);
//                     patternOfReadInWindow[curWindow].push_back(curPattern);
//                     presentOfReadInWindow[curWindow].push_back(present);
//                     //std::cout<<"Put Window "<<curWindow<<", Read "<<read<<" : "<<curPattern<<", "<<present<<std::endl;
//                 }
// //                 std::cout<<"\t"<<curWindow<<":"<<present<<"/"<<curPattern;
//                 lastWindow = curWindow;
//                 present = 0;
//                 curPattern = 0;
//             }
//             curPattern |= (alleles[read][rpos] << offset);
//             present |= (1U << offset);
//         }
// //         std::cout<<std::endl;
//     }
//     
//     // infer local haplotype patterns for each window
//     std::vector<std::vector<uint32_t>> patternsInWindow;
//     std::vector<std::vector<uint32_t>> patternMultiplicityInWindow;
//     for (uint32_t window = 0; window < numWindows; window++) {
//         /* determine relevant patterns: find pattern with highest count, until <ploidy> many patterns
//          * are chosen. To chose the same pattern k times, it must be more than k times more frequent
//          * than all other patterns. */
//         patternsInWindow.push_back(std::vector<uint32_t>());
//         patternMultiplicityInWindow.push_back(std::vector<uint32_t>());
//         std::vector<uint32_t> timesChosen(numPatterns, 0);
//         for (uint32_t j = 0; j < ploidy; j++) {
//             uint32_t maxCount = 0;
//             uint32_t maxPattern = 0;
//             for (uint32_t pattern = 0; pattern < numPatterns; pattern++) {
//                 if (patternCount[window][pattern] / (timesChosen[pattern]+1) > maxCount) {
//                     maxCount = patternCount[window][pattern] / (timesChosen[pattern]+1);
//                     maxPattern = pattern;
//                 }
//             }
//             timesChosen[maxPattern]++;
//         }
//         // write chosen patterns into 2D vector
// //         std::cout<<"Window "<<window<<": ";
//         for (uint32_t pattern = 0; pattern < numPatterns; pattern++) {
// //             if (pattern % 4 == 0)
// //                 std::cout<<std::endl;
// //             std::cout<<"\t"<<patternCount[window][pattern]<<" ("<<timesChosen[pattern]<<")";
//             if (timesChosen[pattern] > 0) {
//                 patternsInWindow[window].push_back(pattern);
//                 patternMultiplicityInWindow[window].push_back(timesChosen[pattern]);
//             }
//         }
// //         std::cout<<std::endl;
//     }
//     
//     // iterate over all windows
//     for (uint32_t window = 0; window < numWindows; window++) {
//         // write patterns as allele vectors
//         std::vector<std::vector<uint32_t>> patternAlleles;
//         uint32_t numLocalPatterns = patternsInWindow[window].size();
//         for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
//             std::vector<uint32_t> vec(windowSize, 0);
//             for (uint32_t pos = 0; pos < windowSize; pos++) {
//                 vec[0] = (patternsInWindow[window][pattern] & (1U << (windowSize - 1 - pos))) > 0;
//             }
//             patternAlleles.push_back(vec);
//         }
//         
//         std::vector<std::vector<double>> prob;
//         double e = errorrate > 0.0 && errorrate < 1.0 ? errorrate : 0.05;
//         // compute probabilities for each read to originate from each pattern
//         for (uint32_t read = 0; read < readsInWindow[window].size(); read++) {
//             std::vector<double> probRead(numLocalPatterns, 0.0);
//             double sum = 0.0;
//             for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
//                 // bit magic: inner xor detects differences between pattern and read, outer operator ensures zero bits on undefined positions for read
//                 uint64_t numMatch = popcount(((patternOfReadInWindow[window][read] ^ patternsInWindow[window][pattern]) ^ (numPatterns-1)) & presentOfReadInWindow[window][read]);
//                 uint64_t numMismatch = popcount((patternOfReadInWindow[window][read] ^ patternsInWindow[window][pattern]) & presentOfReadInWindow[window][read]);
//                 if (numMatch + numMismatch <= 0 || numMatch + numMismatch > windowSize) {
//                     std::cout<<"Invalid match/mismatch count: Window "<<window<<", Read "<<readsInWindow[window][read]<<"("<<read<<"): "<<numMatch<<" and "<<numMismatch<<std::endl;
// //                     std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" ^ "<<(numPatterns-1)<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
// //                     std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
//                     continue;
//                 }
// //                 std::cout<<patternOfReadInWindow[window][read]<<" ^ "<<patternsInWindow[window][pattern]<<" ^ "<<(numPatterns-1)<<" & "<<presentOfReadInWindow[window][read]<<std::endl;
//                 double factor = patternMultiplicityInWindow[window][pattern] * std::pow(e, numMismatch) * std::pow(1-e, numMatch);
//                 probRead[pattern] = factor;
//                 sum += factor;
//                 if (factor == 0.0) {
//                     std::cout<<"Factor was zero for "<<patternOfReadInWindow[window][read]<<" ("<<presentOfReadInWindow[window][read]<<") and "<<patternsInWindow[window][pattern]<<std::endl;
//                 }
//             }
//             for (uint32_t pattern = 0; pattern < numLocalPatterns; pattern++) {
//                 probRead[pattern] /= sum;
//             }
//             prob.push_back(probRead);
//         }
//         
//         // compute scoring for all read pairs in current window
//         for (uint32_t read1 = 0; read1 < readsInWindow[window].size(); read1++) {
//             for (uint32_t read2 = read1 + 1; read2 < readsInWindow[window].size(); read2++) {
//                 uint32_t readId1 = readsInWindow[window][read1];
//                 uint32_t readId2 = readsInWindow[window][read2];
//                 if (overlaps.get(readId1, readId2) < minOverlap)
//                     continue;
//                 double pSame = 0.0;
//                 double pDiff = 0.0;
//                 for (uint32_t p1 = 0; p1 < numLocalPatterns; p1++) {
//                     for (uint32_t p2 = 0; p2 < numLocalPatterns; p2++) {
//                         if (p1 == p2) {
//                             pSame += prob[read1][p1] * prob[read2][p2];
//                         } else {
//                             pDiff += prob[read1][p1] * prob[read2][p2];
//                         }
//                     }
//                 }
//                 if (std::abs(1.0 - (pSame + pDiff)) > 0.0001) {
//                     std::cout<<"pSame + pDiff = "<<(pDiff+pSame)<<" in window "<<window<<std::endl;
//                 } else {
//                     if (pSame == 0)
//                         result->set(readId1, readId2, -std::numeric_limits<float>::infinity());
//                     else if (pDiff == 0) {
//                         result->set(readId1, readId2, std::numeric_limits<float>::infinity());
// //                         std::cout<<window<<" : "<<read1<<"-"<<read2<<" : "<<readId1<<"-"<<readId2<<" : "<<patternOfReadInWindow[window][read1]<<"-"<<patternOfReadInWindow[window][read1]<<" : "<<presentOfReadInWindow[window][read1]<<"-"<<presentOfReadInWindow[window][read1]<<std::endl;
//                     }
//                     else
//                         result->set(readId1, readId2, result->get(readId1, readId2) + log(pSame / pDiff));
//                 }
//             }
//         }
//     }
}

void ReadScoring::computeStartEnd (const ReadSet* readset,
                                   std::vector<uint32_t>& begins,
                                   std::vector<uint32_t>& ends,
                                   std::vector<std::vector<uint32_t>>& positions,
                                   std::vector<std::vector<uint32_t>>& alleles,
                                   std::vector<uint32_t>& posList,
                                   std::unordered_map<uint32_t, uint32_t>& posMap,
                                   uint32_t& longestReadSpan) const {
    // copy all relevant information from the readset into vectors for efficient access
    uint32_t numReads = readset->size();
    std::unordered_set<uint32_t> allPos;
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
    posList.clear();
    for (uint32_t pos : allPos) {
        posList.push_back(pos);
    }
    std::sort(posList.begin(), posList.end());
    for (uint32_t i = 0; i < posList.size(); i++) {
        posMap[posList[i]] = i;
    }
    
    // determine length of the longest read
    longestReadSpan = 0;
    for (uint32_t i = 0; i < numReads; i++) {
        longestReadSpan = std::max(longestReadSpan, ends[i] - begins[i]);
    }
}

void ReadScoring::computeOverlapDiff (const ReadSet* readset,
                                      const std::vector<uint32_t>& begins,
                                      const std::vector<uint32_t>& ends,
                                      const std::vector<std::vector<uint32_t>>& positions,
                                      const std::vector<std::vector<uint32_t>>& alleles,
                                      const std::vector<uint32_t>& posList,
                                      std::unordered_map<uint32_t, uint32_t>& posMap,
                                      TriangleSparseMatrix& overlaps,
                                      TriangleSparseMatrix& diffs,
                                      double& distSame,
                                      double& distDiff,
                                      const uint32_t minOverlap,
                                      const uint32_t ploidy,
                                      const uint32_t longestReadSpan,
                                      const uint32_t begin,
                                      const uint32_t end) const {
                                          
    std::vector<uint32_t> coveredReads;
    
    // determine coveredReads if begin and end position are specified (otherwise use all reads)
    if (begin == 0 && end == begins.size()) {
        for (uint32_t i = 0; i < begins.size(); i++) {
            coveredReads.push_back(i);
        }
    } else {
        // binary search to find first read, which can theoretically cover this window
        uint32_t firstIndex = std::lower_bound(begins.begin(), begins.end(), begin-longestReadSpan) - begins.begin();
        // iterate until start position of read is behind required start
        for (uint32_t j = firstIndex; begins[j] <= begin && j < begins.size(); j++) {
            if (ends[j] >= end) {
                coveredReads.push_back(j);
            }
        }
    }
    
    // iterate over all read pairs (efficiently omitting those who can certainly not overlap)
    std::vector<double> relativeDiffs;
    for (uint32_t i = 0; i < coveredReads.size(); i++) {
        // iterate until start position of read is behind required start
        uint32_t ci = coveredReads[i];
        for (uint32_t j = i+1; j < coveredReads.size() && begins[coveredReads[j]] <= ends[ci]; j++) {
            uint32_t cj = coveredReads[j];
            if (ends[ci] < begins[cj] || ends[cj] < begins[ci])
                continue;
            uint32_t ov = 0;
            uint32_t di = 0;
            uint32_t k = 0;
            uint32_t l = 0;
            while (k < positions[ci].size() && l < positions[cj].size()) {
                if (positions[ci][k] == positions[cj][l]) {
                    if (alleles[ci][k] != alleles[cj][l]) {
                        di++;
                    }
                    ov++; k++; l++; 
                } else if (positions[ci][k] < positions[cj][l]) {
                    k++;
                } else {
                    l++;
                }
            }
            if (ov >= minOverlap) {
                overlaps.set(ci, cj, ov);
                diffs.set(ci, cj, di);
                relativeDiffs.push_back((double)di / (double)ov);
            }
        }
    }
    
    // estimate error rate of reads within equal haplotype
    computeCutoff(coveredReads, ploidy, relativeDiffs, distSame, distDiff);
}

void ReadScoring::computeOverlapDiff (const ReadSet* readset,
                                      const std::vector<uint32_t>& begins,
                                      const std::vector<uint32_t>& ends,
                                      const std::vector<std::vector<uint32_t>>& positions,
                                      const std::vector<std::vector<uint32_t>>& alleles,
                                      const std::vector<uint32_t>& posList,
                                      std::unordered_map<uint32_t, uint32_t>& posMap,
                                      TriangleSparseMatrix& overlaps,
                                      TriangleSparseMatrix& diffs,
                                      double& distSame,
                                      double& distDiff,
                                      const uint32_t minOverlap,
                                      const uint32_t ploidy,
                                      const uint32_t longestReadSpan) const {
    computeOverlapDiff(readset, begins, ends, positions, alleles, posList, posMap, overlaps, diffs, distSame, distDiff, minOverlap, ploidy, longestReadSpan, 0, begins.size());
}

void ReadScoring::computeCutoff(const std::vector<uint32_t>& coveredReads, const uint32_t ploidy, std::vector<double> relDiffs, double& distSame, double& distDiff) const {
    uint32_t numReads = coveredReads.size();
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
