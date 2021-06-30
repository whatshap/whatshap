#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <set>
#include <algorithm> 
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <unistd.h>
#include <map>
#include "../BMEAN/utils.h"
#include "../CTPL/ctpl_stl.h"
#include "robin_hood.h"

struct POASeq {
	std::string seq;
	int beg;
	int end;

	bool operator<(const POASeq& s2) const {
		if (beg < s2.beg) {
			return true;
		} else if (beg == s2.beg and end < s2.end) {
			return true;
		} else {
			return false;
		}
	}

	POASeq() {
		
	}

	POASeq(std::string s, int b, int e) {
		seq = s;
		beg = b;
		end = e;
	}
};

std::string polishCorrection(std::string correctedRead, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh);

// std::vector<std::pair<std::string, std::string>> polishCorrection(std::string correctedRead, std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, robin_hood::unordered_map<std::string, unsigned>& pilesMers, unsigned merSize, int solidThresh, int minGap, int maxGap);

// std::vector<std::pair<std::string, std::string>> polishCorrection(std::string correctedRead, std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, unsigned merSize, int solidThresh, int minGap, int maxGap);

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, robin_hood::unordered_map<std::string, unsigned>& merCounts, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize);

std::string alignConsensuses(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<robin_hood::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::string>& templates, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize);

void processReads(std::vector<std::vector<std::string>>& reads, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, std::string path);

void runCorrection(std::string PAFIndex, std::string alignmentFile, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, std::string proofFile, unsigned maxMSA, std::string path);
