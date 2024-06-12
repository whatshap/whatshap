#include "correctionDBG.h"

robin_hood::unordered_map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize) {
	robin_hood::unordered_map<std::string, std::vector<unsigned>> mers;

	for (unsigned i = 0; i < sequence.length() - merSize + 1; i++) {
			mers[sequence.substr(i, merSize)].push_back(i);
	}

	return mers;
}

int getNextSrc(std::string correctedRead, unsigned beg, unsigned merSize) {
	unsigned nb = 0;
	unsigned i = beg;

	while (i < correctedRead.length() and (isUpperCase(correctedRead[i]) or nb < merSize)) {
		if (isUpperCase(correctedRead[i])) {
			nb++;
		} else {
			nb = 0;
		}
		i++;
	}

	return nb >= merSize ? i - 1 : -1;
}

int getNextDst(std::string correctedRead, unsigned beg, unsigned merSize) {
	unsigned nb = 0;
	unsigned i = beg;

	while (i < correctedRead.length() and nb < merSize) {
		if (isUpperCase(correctedRead[i])) {
			nb++;
		} else {
			nb = 0;
		}
		i++;
	}

	return nb >= merSize ? i - 1 : -1;
}


// Anchors without repeated k-mers
std::vector<std::pair<std::string, std::string>> getAnchors(robin_hood::unordered_map<kmer, unsigned>& merCounts, std::string srcZone, std::string dstZone, unsigned merSize, unsigned nb) {
	std::vector<std::pair<std::string, std::string>> res;
	unsigned i;

	robin_hood::unordered_map<std::string, std::vector<unsigned>> mersPosSrc = getKMersPos(srcZone, merSize);
	robin_hood::unordered_map<std::string, std::vector<unsigned>> mersPosDst = getKMersPos(dstZone, merSize);

	// Consider all k-mers of the src zone as potential anchors
	std::vector<std::string> candidatesSrc(srcZone.size() - merSize + 1);
	for (i = 0; i < srcZone.size() - merSize + 1; i++) {
		candidatesSrc[i] = srcZone.substr(i, merSize);
	}
	// Same with the dst zone
	std::vector<std::string> candidatesDst(dstZone.size() - merSize + 1);
	for (i = 0; i < dstZone.size() - merSize + 1; i++) {
		candidatesDst[i] = dstZone.substr(i, merSize);
	}

	// Add the anchors pairs to the result vector, without allowing repeated k-mers
	for (std::string csrc : candidatesSrc) {
		if (mersPosSrc[csrc].size() == 1) {
			for (std::string cdst : candidatesDst) {
				if (mersPosDst[cdst].size() == 1) {
					res.push_back(std::make_pair(csrc, cdst));
				}
			}
		}
	}

	// Sort the anchors vector in ascending order of the number of occurrences of (src + dst)
	std::sort(res.begin(), res.end(),
		[&merCounts](std::pair<std::string, std::string>& r1, std::pair<std::string, std::string>& r2) {
			int occ1 = merCounts[str2num(r1.first)] + merCounts[str2num(r1.second)];
			int occ2 = merCounts[str2num(r2.first)] + merCounts[str2num(r2.second)];
			return occ1 > occ2;
		}
	);

	std::vector<std::pair<std::string, std::string>> finalRes;
	for (i = 0; i < nb and i < res.size(); i++) {
		finalRes.push_back(res[i]);
	}

	return finalRes;
}

std::string polishCorrection(std::string correctedRead, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh) {
	std::set<std::string> visited;
	unsigned curBranches;
	unsigned dist;
	std::string curExt;
	std::string correctedRegion;
	unsigned maxSize;
	unsigned maxBranches = 50;
	std::vector<std::pair<std::string, std::string>> corList;
	int zone = 3;
	int srcBeg, srcEnd, dstBeg, dstEnd;
	unsigned tmpSrcBeg = 0, tmpSrcEnd = 0, tmpDstBeg = 0, tmpDstEnd = 0;
	std::string src, dst;
	std::pair<int, int> pos;
	std::vector<std::pair<std::string, std::string>> anchors;
	unsigned anchorNb;
	std::string srcZone, dstZone;
	robin_hood::unordered_map<std::string, std::vector<unsigned>> srcPos, dstPos;
	std::string oldCorrectedRead;
	int b, l;
	std::string r, c;

	// Skip uncorrected head of the read
	unsigned i = 0;
	while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
		i++;
	}

	if (i > 0 and i < correctedRead.length() and correctedRead.length() - i >= merSize) {
		int extLen = i;
		oldCorrectedRead = correctedRead;
		correctedRead = correctedRead.substr(i);
		int extSize = extendLeft(merCounts, merSize, extLen, correctedRead, solidThresh);
		if (extSize < extLen) {
			correctedRead = oldCorrectedRead.substr(0, extLen - extSize) + correctedRead;
			i = i - (extLen - extSize);
		}
	}

	// Search for poorly supported regions bordered by solid corrected regions
	while (i < correctedRead.length()) {
		srcEnd = getNextSrc(correctedRead, i, merSize + zone);
		dstEnd = getNextDst(correctedRead, srcEnd + 1, merSize + zone);
		srcBeg = srcEnd - merSize - zone + 1;
		dstBeg = dstEnd - merSize - zone + 1;

		// Polish the poorly supported region region if 2 anchors were found
		if (srcEnd != -1 and dstEnd != -1) {
			correctedRegion = "";
			srcZone = correctedRead.substr(srcBeg, merSize + zone);
			dstZone = correctedRead.substr(dstBeg, merSize + zone);
			anchors = getAnchors(merCounts, srcZone, dstZone, merSize, 5);
			srcPos = getKMersPos(srcZone, merSize);
			dstPos = getKMersPos(dstZone, merSize);

			// Attempt to link frequent anchors
			anchorNb = 0;
			while (anchorNb < anchors.size() and correctedRegion.empty()) {
				src = anchors[anchorNb].first;
				dst = anchors[anchorNb].second;
				tmpSrcBeg = srcBeg + srcPos[src][0];
				tmpSrcEnd = tmpSrcBeg + merSize - 1;
				tmpDstBeg = dstBeg + dstPos[dst][0];
				tmpDstEnd = tmpDstBeg + merSize - 1;
				
				if (src != dst) {
					curBranches = 0;
					dist = 0;
					curExt = src;
					correctedRegion = "";
					maxSize = 15.0 / 100.0 * 2.0 * (tmpDstBeg - tmpSrcEnd - 1) + (tmpDstBeg - tmpSrcEnd - 1) + merSize;
					link(merCounts, src, dst, merSize, visited, &curBranches, dist, curExt, correctedRegion, merSize, maxSize, maxBranches, solidThresh, merSize);
				}
				anchorNb++;
			}

			if (!correctedRegion.empty()) {
				// Anchor the correction to the read
				r = correctedRead.substr(tmpSrcBeg, tmpDstEnd - tmpSrcBeg + 1);
				c = correctedRegion;
				b = correctedRead.find(r);
				l = r.length();
				if ((int) b != -1) {
					correctedRead.replace(b, l, c);
					i = b;
				} else {
					i = tmpDstBeg > i ? tmpDstBeg : dstBeg;	
				}
			} else {
				i = tmpDstBeg > i ? tmpDstBeg : dstBeg;
			}
		} else {
			i = correctedRead.length();	
		}
	}

	i = correctedRead.length() - 1;
	while (i > 0 and !isUpperCase(correctedRead[i])) {
		i--;
	}

	if (i > 0 and i < correctedRead.length() - 1 and i + 1 >= merSize) {
		int extLen = correctedRead.length() - 1 - i;
		oldCorrectedRead = correctedRead;
		correctedRead = correctedRead.substr(0, i + 1);
		int extSize = extendRight(merCounts, merSize, extLen, correctedRead, solidThresh);
		if (extSize < extLen) {
			correctedRead = correctedRead + oldCorrectedRead.substr(oldCorrectedRead.length() - (extLen - extSize), extLen - extSize);
		}
	}

	return correctedRead;
}