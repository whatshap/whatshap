#include <string>
#include <vector>
#include <set>
#include <locale>
#include <unordered_map>
#include "robin_hood.h"

using namespace std;

vector<string> splitString(string s, string delimiter);

std::vector<bool> fullstr2num(const string& str);

std::string fullnum2str(vector<bool> num);

bool isUpperCase(char c);

int nbCorBases(std::string correctedRead);

bool dropRead(std::string correctedRead);

std::string toUpperCase(std::string& s, int beg, int end);

std::string toLowerCase(std::string& s, int beg, int end);

std::string trimRead(std::string correctedRead, unsigned merSize);

std::vector<std::string> splitRead(std::string name, std::string correctedRead, std::vector<std::pair<unsigned, unsigned>>& pilesPos, unsigned windowSize, unsigned windowOverlap);

void indexReads(robin_hood::unordered_map<std::string, std::vector<bool>>& index, std::string readsFile);
