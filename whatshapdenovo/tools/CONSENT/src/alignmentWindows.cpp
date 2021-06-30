#include "alignmentWindows.h"
#include <unistd.h>
#include <iostream>

unsigned* getCoverages(std::vector<Overlap>& alignments) {
	unsigned tplLen = alignments.begin()->qLength;
	// unsigned* coverages = new unsigned[tplLen];
	unsigned* coverages = (unsigned*) calloc(tplLen, sizeof(int));
	unsigned i;
	// for (i = 0; i < tplLen; i++) {
	// 	coverages[i] = 1;
	// }
	unsigned beg, end;

	for (Overlap al : alignments) {
		beg = al.qStart;
		end = al.qEnd;

		for (i = beg; i <= end; i++) {
			coverages[i]++;
		}
	}

	return coverages;
}

std::vector<std::pair<unsigned, unsigned>> getAlignmentWindowsPositions(unsigned tplLen, std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, int overlappingWindows) {
	unsigned* coverages = getCoverages(alignments);
	unsigned i;
	unsigned beg = 0;
	unsigned end = tplLen - 1;

	std::vector<std::pair<unsigned, unsigned>> pilesPos;

	unsigned curLen = 0;
	beg = 0;
	end = 0;
	i = 0;
	while (i < tplLen) {
		if (curLen >= windowSize) {
			pilesPos.push_back(std::make_pair(beg, beg + curLen - 1));
			if (overlappingWindows) {
				i = i - overlappingWindows;
			}
			beg = i;
			curLen = 0;
		}
		if (coverages[i] < minSupport) {
			curLen = 0;
			i++;
			beg = i;
		} else {
			curLen++;
			i++;
		}
	}

	// Special case for the last window
	int pushed = 0;
	beg = 0;
	end = tplLen - 1;
	curLen = 0;
	i = tplLen - 1;
	while (i > 0 and !pushed) {
		if (curLen >= windowSize) {
			pilesPos.push_back(std::make_pair(end - curLen + 1, end));
			pushed = 1;
			end = i;
			curLen = 0;
		}
		if (coverages[i] < minSupport) {
			curLen = 0;
			i--;
			end = i;
		} else {
			curLen++;
			i--;
		}
	}

	// delete [] coverages;
	free(coverages);

	return pilesPos;
}

std::vector<std::string> getAlignmentWindowsSequences(std::vector<Overlap>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string>& sequences, unsigned qBeg, unsigned end, unsigned merSize, unsigned maxSupport, unsigned commonKMers) {
	std::vector<std::string> curPile;
	std::vector<unsigned> curScore;
	unsigned length, shift;
	length = end - qBeg + 1;
	unsigned curPos = 0;
	unsigned tBeg, tEnd;

	if (qBeg + length - 1 >= sequences[alignments.begin()->qName].length()) {
		return curPile;
	}

	// Insert template sequence
	curPile.push_back(sequences[alignments.begin()->qName].substr(qBeg, length));

	Overlap al;
	std::string tmpSeq;
	// Insert aligned sequences
	while (curPos < alignments.size()) {// and curPile.size() < maxSupport) {
		al = alignments[curPos];
		tBeg = al.tStart;
		tEnd = al.tEnd;
		length = end - qBeg + 1;
		if (qBeg > al.qStart) {
			shift = qBeg - al.qStart;
		} else {
			shift = 0;
		}
		
		// For all alignments than span, or begin/end in the query window
		if ( ((al.qStart <= qBeg and al.qEnd > qBeg) or (end <= al.qEnd and al.qStart < end)) and al.tStart + shift <= al.tEnd) {

			if (qBeg < al.qStart and al.qEnd < end) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				length = tEnd - tBeg + 1;
			} else if (qBeg < al.qStart) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				length = std::min((int) length, std::min((int) al.tLength - 1, (int) tBeg + (int) length - 1) - (int) tBeg + 1);
			} else if (al.qEnd < end) {
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				length = std::min((int) length, (int) tEnd - std::max(0, (int) tEnd - (int) length + 1) + 1);
			}

			tmpSeq = sequences[al.tName].substr(tBeg, tEnd - tBeg + 1);
			if (al.strand) {
				tmpSeq = rev_comp::run(tmpSeq);
			}

			tmpSeq = tmpSeq.substr(shift, length);

			// Default
			if (tmpSeq.length() >= merSize) {
				curPile.push_back(tmpSeq);
			}
		}
		curPos++;
	}

	return curPile;
}

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentWindows(std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string> sequences, unsigned merSize, unsigned commonKMers) {
	unsigned tplLen = alignments.begin()->qLength;

	std::vector<std::pair<unsigned, unsigned>> windowsPos = getAlignmentWindowsPositions(tplLen, alignments, minSupport, maxSupport, windowSize, windowOverlap);

	std::vector<std::vector<std::string>> windowsSeq;

	for (std::pair<int, int> p : windowsPos) {
		windowsSeq.push_back(getAlignmentWindowsSequences(alignments, minSupport, windowSize, windowOverlap, sequences, p.first, p.second, merSize, maxSupport, commonKMers));
	}

	return std::make_pair(windowsPos, windowsSeq);
}
