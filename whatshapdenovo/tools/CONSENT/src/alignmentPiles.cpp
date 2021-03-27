#include <fstream>
#include <iostream>
#include "alignmentPiles.h"

robin_hood::unordered_map<std::string, std::string> getSequencesMap(std::vector<Overlap>& alignments, robin_hood::unordered_map<std::string, std::vector<bool>>& readIndex) {
	robin_hood::unordered_map<std::string, std::string> sequences;
	std::string header, seq;

	// Insert template sequence
	sequences[alignments.begin()->qName] = fullnum2str(readIndex[alignments.begin()->qName]);

	// Insert aligned sequences
	for (Overlap al : alignments) {
		if (sequences[al.tName] == "") {
			sequences[al.tName] = fullnum2str(readIndex[al.tName]);
		}
	}

	return sequences;
}

std::vector<Overlap> getNextReadPile(std::ifstream& f, unsigned maxSupport) {
	std::vector<Overlap> curReadAlignments;
	std::string line, curRead;
	Overlap curAl;

	getline(f, line);
	while(line.length() > 0 or !curReadAlignments.empty()) {
		if (line.length() > 0) {
			curAl = Overlap(line);
		}
		if (line.length() > 0 and (curRead == "" or curAl.qName == curRead)) {
			curRead = curAl.qName;

			// Only keep MAX best overlaps
			if (curReadAlignments.size() < maxSupport) {
				curReadAlignments.push_back(curAl);
			}

			getline(f, line);
		} else {
			if (!f.eof()) {
				f.seekg(-line.length()-1, f.cur);
			}
			return curReadAlignments;
		}
	}

	return curReadAlignments;
}