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
#include "utils.h"
#include "Overlap.h"
#include "robin_hood.h"

std::pair<std::string, std::string> processContig(std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned maxMSA, std::string path, unsigned nbThreads);

void runCorrection(std::string PAFIndex, std::string alignmentFile, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, std::string proofFile, unsigned maxMSA, std::string path);
