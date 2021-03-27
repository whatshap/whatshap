#include <vector>
#include <string>
#include <unordered_map>
#include "Overlap.h"
#include "utils.h"
#include "robin_hood.h"

std::vector<Overlap> getNextReadPile(std::ifstream& f, unsigned maxSupport);

robin_hood::unordered_map<std::string, std::string> getSequencesMap(std::vector<Overlap>& alignments, robin_hood::unordered_map<std::string, std::vector<bool>>& readIndex);
