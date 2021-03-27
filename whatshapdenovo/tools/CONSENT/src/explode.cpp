#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>

std::string getId(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end;

    pos_end = s.find (delimiter, pos_start);
    return s.substr (pos_start, pos_end - pos_start);
}

void explodeFile(std::ifstream& f, std::string exploded) {
	std::string line, curRead, oldRead;
	std::set<std::string> readsList;
	int nbFiles = 1;
	std::ofstream curFile;
  	curFile.open((exploded + "_" + std::to_string(nbFiles)).c_str());
  	std::string curStr;


	getline(f, line);
	while(line.length() > 0) {
		oldRead = curRead;
		curRead = getId(line, "\t");
		if (oldRead == "" or curRead == oldRead) {
			curStr += line;
			curStr += "\n";
			getline(f, line);
		} else {
			readsList.insert(oldRead);
			if (readsList.find(curRead) != readsList.end()) {
				readsList.clear();
				curFile.close();
				nbFiles++;
				curFile.open((exploded + "_" + std::to_string(nbFiles)).c_str());
			}
			curFile << curStr;
			curStr = line;
			curStr += "\n";
			getline(f, line);
		}
	}
	if (!curStr.empty()) {
		curFile << curStr;
		curFile.close();
	}
}

int main (int argc, char* argv[]) {
	std::ifstream alignments(argv[1]);
	std::string curTpl;

	explodeFile(alignments, argv[2]);

	alignments.close();
}