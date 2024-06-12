#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
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

std::pair<std::string, std::string> processRead(int id, std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned maxMSA, std::string path);

void runCorrection(std::string PAFIndex, std::string alignmentFile, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, std::string proofFile, unsigned maxMSA, std::string path);
