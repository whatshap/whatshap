#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <string.h>

int cmpFields(const std::pair<std::string, unsigned long>& f1, const std::pair<std::string, unsigned long>& f2) {
	if (f1.first == f2.first and f1.second == f2.second) {
		return 0;
	} else if (f1.first < f2.first or (f1.first == f2.first and f1.second > f2.second)) {
		return -1;
	} else {
		return 1;
	}
}

std::pair<std::string, unsigned long> splitString(bool correction, std::string s, std::string delimiter) {
    std::pair<std::string, unsigned long> f;

    std::string token;
    std::string qName;
    std::string tName;
	std::stringstream iss(s);
	getline(iss, qName, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	getline(iss, tName, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	getline(iss, token, '\t');
	f.second = stoul(token);

	if (correction) {
		f.first = qName;
	} else {
		f.first = tName;
	}

	return f;
}

std::vector<std::string> getFilesNames(char** args) {
	std::vector<std::string> res;

	int i = 3;
	while (args[i] != NULL) {
		res.push_back(args[i]);
		i++;
	}

	return res;
}

void getNextLines(std::vector<std::ifstream*>& files, std::vector<std::string>& lines, unsigned long size) {
	std::string line;
	for (unsigned long i = 0; i < size; i++) {
		std::ifstream* f = files[i];
		getline(*f, line);
		lines[i] = std::string(line);
	}
}

std::pair<std::string, unsigned long> getFields(bool correction, std::string line) {
	std::pair<std::string, unsigned long> f;

	if (!line.empty()) {
		f = splitString(correction, line, "\t");
	} else {
		f.first = "";
		f.second = 0;
	}

	return f;
}

void getNextFields(bool correction, std::vector<std::pair<std::string, unsigned long>>& fields, std::vector<std::string>& lines, unsigned long size) {
	for (unsigned long i = 0; i < size; i++) {
		fields[i] = getFields(correction, lines[i]);
	}
}

bool areAllEmpty(std::vector<std::string>& lines, unsigned long size) {
	for (unsigned long i = 0; i < size; i++) {
		if (!lines[i].empty()) {
			return false;
		}
	}

	return true;
}

unsigned long getSmallestLine(std::vector<std::pair<std::string, unsigned long>>& fields, unsigned long size, int comp(const std::pair<std::string, unsigned long>&, const std::pair<std::string, unsigned long>&)) {
	unsigned long index = 0;
	std::pair<std::string, unsigned long> smallest = fields[0];

	while (fields[index].first.empty()) {
		index++;
	}
	smallest = fields[index];

	for (unsigned long i = index; i < size; i++) {
		if (!fields[i].first.empty() and comp(fields[i], smallest) <= 0) {
			smallest = fields[i];
			index = i;
		}
	}

	return index;
}

int main (int argc, char* argv[]) {
	std::ofstream outFile;
	outFile.open(argv[1]);

	bool correction = strcmp(argv[2], "correction") == 0 ? 1 : 0;

	std::vector<std::string> filesNames = getFilesNames(argv);
	std::vector<std::ifstream*> files(filesNames.size());
	std::vector<std::pair<std::string, unsigned long>> fields(filesNames.size());

	std::ifstream* f;
	for (unsigned long i = 0; i < filesNames.size(); i++) {
		f = new std::ifstream;
		files[i] = f;
		files[i]->open(filesNames[i]);
	}

	std::vector<std::string> lines(filesNames.size());
	getNextLines(files, lines, lines.size());
	getNextFields(correction, fields, lines, lines.size());
	bool over = areAllEmpty(lines, lines.size());
	unsigned long smallest;
	std::string line;

	while (!over) {
		smallest = getSmallestLine(fields, fields.size(), cmpFields);

		outFile << lines[smallest] << std::endl;

		getline(*files[smallest], line);
		lines[smallest] = std::string(line);
		fields[smallest] = getFields(correction, line);

		over = areAllEmpty(lines, lines.size());
	}


	for (std::ifstream* f : files) {
		f->close();
	}
}
