#include <algorithm>
#include "correctionAlignment.h"
#include "utils.h"
#include "../BMEAN/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"

int nbSolidMers(std::string seq, robin_hood::unordered_map<kmer, unsigned> merCounts, unsigned merSize, unsigned solidThresh) {
	int nb = 0;
	for (unsigned i = 0; i < seq.length() - merSize + 1; i++) {
		if (merCounts[str2num(seq.substr(i, merSize))] >= solidThresh) {
			nb++;
		}
	}

	return nb;
}

int nbUpperCase(std::string s) {
	int nb = 0;
	for (unsigned i = 0; i < s.length(); i++) {
		if (isUpperCase(s[i])) {
			nb++;
		}
	}

	return nb;
}

std::pair<int, int> getIndels(std::string cigar){
	int ins = 0;
	int del = 0;
	int current = 0;
	for(unsigned i = 0; i < cigar.length(); i++){
		if('0' <= cigar[i] && cigar[i] <= '9'){
			current = (current * 10) + (cigar[i] - '0');
		} else {
			if (cigar[i] == 'I') {
				ins += current;
			} else if (cigar[i] == 'D') {
				del += current;
			}
			current = 0;
		}
	}
	return std::make_pair(ins, del);
}

std::string alignConsensus(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<robin_hood::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::string>& templates, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize) {
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	StripedSmithWaterman::Alignment subAlignment;	
	int32_t maskLen = 15;
	unsigned beg, end, oldEnd;
	oldEnd = 0;
	std::string outSequence;
	outSequence = sequence;
	std::transform(outSequence.begin(), outSequence.end(), outSequence.begin(), ::tolower);

	std::string corWindow;
	unsigned i = 0;
	std::string tmpSequence, consUp;
	int curPos = startPos;
	int alPos;
	int sizeAl;
	std::string curCons, oldCons;
	robin_hood::unordered_map<kmer, unsigned> oldMers;
	robin_hood::unordered_map<kmer, unsigned> curMers;
	unsigned overlap;
	std::string seq1, seq2;
	int solidMersSeq1, solidMersSeq2;
	std::pair<int, int> indels;
	unsigned ins, del;

	for (i = 0; i < consensuses.size(); i++) {
		if (consensuses[i].length() < merSize) {
			curCons = templates[i];
			std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::tolower);
		} else {
			curCons = consensuses[i];
		}
		curMers = merCounts[i];
		
		alPos = std::max(0, (int) curPos - (int) windowOverlap);
		if (alPos + windowSize + 2 * windowOverlap >= outSequence.length()) {
			sizeAl = outSequence.length() - alPos;
		} else {
			sizeAl = windowSize + 2 * windowOverlap;
		}

		aligner.Align(curCons.c_str(), outSequence.c_str() + alPos, sizeAl, filter, &alignment, maskLen);
		beg = alignment.ref_begin + alPos;
		end = alignment.ref_end + alPos;
		curCons = curCons.substr(alignment.query_begin, alignment.query_end - alignment.query_begin + 1);

		// Check if alignment positions overlap the previous window. If they do, chose the best subsequence
		if (i != 0 and oldEnd >= beg) {
 			overlap = oldEnd - beg + 1;
			if (consensuses[i].length() >= merSize and oldCons.length() >= overlap and curCons.length() >= overlap) {
				seq1 = oldCons.substr(oldCons.length() - 1 - overlap + 1, overlap);
				seq2 = curCons.substr(0, overlap);
				if (toUpperCase(seq1, 0, seq1.length() - 1) != toUpperCase(seq2, 0, seq2.length() - 1)) {
					if (overlap >= merSize) {
						solidMersSeq1 = nbSolidMers(seq1, oldMers, merSize, solidThresh);
						solidMersSeq2 = nbSolidMers(seq2, curMers, merSize, solidThresh);
					} else {
						solidMersSeq1 = nbUpperCase(seq1);
						solidMersSeq2 = nbUpperCase(seq2);
					}
					if (solidMersSeq1 > solidMersSeq2) {
						aligner.Align(seq1.c_str(), seq2.c_str(), std::min(seq1.length(), seq2.length()), filter, &subAlignment, maskLen);
						indels = getIndels(subAlignment.cigar_string);
						ins = indels.first;
						del = indels.second;
						if (overlap - ins + del < curCons.length()) {
							curCons = seq1 + curCons.substr(overlap - ins + del);
						} else {
							curCons = "";
						}
					}
				}
			}
		}

		if (curCons != "") {
			if (consensuses[i].length() >= merSize){
				consUp = curCons;
				std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::toupper);
				outSequence.replace(beg, end - beg + 1, consUp);
			}
			if (i < consensuses.size() - 1) {
				curPos = curPos + pilesPos[i+1].first - pilesPos[i].first - (end - beg + 1) + curCons.length() ;
				oldCons = curCons;
				oldMers = merCounts[i];
				oldEnd = beg + curCons.length() - 1;
			}
		}
	}

	return outSequence;
}