#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <string.h>

std::vector<std::string> splitString(std::string s, char delimiter) {
    std::vector<std::string> res;

    std::string token;
	std::stringstream iss(s);
	while (getline(iss, token, delimiter)) {
		res.push_back(token);
	}

	return res;
}

std::string processLine(std::string line) {
	std::vector<std::string> v;
	v = splitString(line, '\t');

	std::string res;
	res = v[5] + "\t" + v[6] + "\t" + v[7] + "\t" + v[8] + "\t" + v[4] + "\t" + v[0] + "\t" + v[1] + "\t" + v[2] + "\t" + v[3];
	for (unsigned i = 9 ; i < v.size(); i++) {
		res += "\t" + v[i];
	}

	return res;
}

int main(int argc, char* argv[]) {
	std::ifstream in;
	in.open(argv[1]);
	std::string line;

	std::ofstream out;
	out.open(argv[2]);

	while(getline(in, line)) {
		out << processLine(line) << std::endl;
	}

	in.close();
	out.close();
}