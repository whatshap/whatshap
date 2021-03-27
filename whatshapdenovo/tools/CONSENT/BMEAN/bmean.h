#ifndef BMEAN
#define BMEAN



#include <vector>
#include <string>
#include "utils.h"
#include <cmath>
#include "robin_hood.h"


using namespace std;





std::pair<std::vector<std::vector<std::string>>, robin_hood::unordered_map<kmer, unsigned>> MSABMAAC(const vector<string>& nadine,uint32_t la,double cuisine, unsigned solidThresh, unsigned minAnchors, unsigned maxMSA, string path);



#endif
