#include <string>
#include <vector>
#include <unordered_map>
#include "../BMEAN/utils.h"
#include "robin_hood.h"

std::string alignConsensus(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<robin_hood::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::string>& templates, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize);
