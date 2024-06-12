#include "correctionMSA.h"
#include "utils.h"
#include "../BMEAN/bmean.h"
#include "correctionDBG.h"

std::string weightConsensus(std::string& consensus, std::vector<std::string>& pile, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, unsigned windowSize, unsigned solidThresh) {
	std::vector<std::string> splits;
	std::string curSplit;

	std::string header = "";
	std::string sequence = "";
	std::string curFct;

	unsigned i = 0;
	while (i < consensus.length() - merSize + 1) {
		curFct = consensus.substr(i, merSize);
		curFct = toUpperCase(curFct, 0, merSize);
		if (merCounts[str2num(curFct)] >= solidThresh) {
			consensus = toUpperCase(consensus, i, i + merSize - 1);
		} else {
			consensus = toLowerCase(consensus, i, i + merSize - 1);
		}
		i++;
	}

	return consensus;
}

std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> computeConsensusReadCorrection(std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path) {
	int bmeanSup;
	bmeanSup = std::min((int) commonKMers, (int) piles.size() / 2);
	std::pair<std::vector<std::vector<std::string>>, robin_hood::unordered_map<kmer, unsigned>> rOut = MSABMAAC(piles, merSize, bmeanSup, solidThresh, minAnchors, maxMSA, path);

	if (rOut.first.size() == 0) {
		return std::make_pair(piles[0], rOut.second);
	}
	auto result = rOut.first;
	auto merCounts = rOut.second;
	std::string corTpl = result[0][0];

	// Polish the consensus
	std::vector<std::pair<std::string, std::string>> corList;
	if (corTpl.length() >= merSize) {
		corTpl = weightConsensus(corTpl, piles, merCounts, merSize, windowSize, solidThresh);
		corTpl = polishCorrection(corTpl, merCounts, merSize, solidThresh);
	}

	return std::make_pair(corTpl, merCounts);
}

std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> computeConsensusAssemblyPolishing(int id, std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path, unsigned nbThreads) {
	int bmeanSup;
	bmeanSup = std::min((int) commonKMers, (int) piles.size() / 2);
	std::pair<std::vector<std::vector<std::string>>, robin_hood::unordered_map<kmer, unsigned>> rOut = MSABMAAC(piles, merSize, bmeanSup, solidThresh, minAnchors, maxMSA, path);

	if (rOut.first.size() == 0) {
		return std::make_pair(piles[0], rOut.second);
	}
	auto result = rOut.first;
	auto merCounts = rOut.second;
	std::string corTpl = result[0][0];

	// Polish the consensus
	std::vector<std::pair<std::string, std::string>> corList;
	if (corTpl.length() >= merSize) {
		corTpl = weightConsensus(corTpl, piles, merCounts, merSize, windowSize, solidThresh);
		corTpl = polishCorrection(corTpl, merCounts, merSize, solidThresh);
	}

	return std::make_pair(corTpl, merCounts);
}