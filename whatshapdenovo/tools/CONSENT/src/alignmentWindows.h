#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include "Overlap.h"
#include "reverseComplement.h"
#include "robin_hood.h"

unsigned* getCoverages(std::vector<Overlap>& alignments);

std::vector<std::pair<unsigned, unsigned>> getAlignmentWindowsPositions(unsigned tplLen, std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, int overlappingWindows);

std::vector<std::string> getAlignmentWindowsSequences(std::vector<Overlap>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string>& sequences, unsigned beg, unsigned end, unsigned merSize, unsigned maxSupport, unsigned commonKMers);

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentWindows(std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string> sequences, unsigned merSize, unsigned commonKMers);
