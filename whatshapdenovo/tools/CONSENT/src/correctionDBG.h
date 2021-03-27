#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include "../BMEAN/utils.h"
#include "utils.h"
#include "DBG.h"
#include "robin_hood.h"

std::string polishCorrection(std::string correctedRead, robin_hood::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh);
