#include "utils.h"
#include <fstream>
#include <iostream>
#include <algorithm>

vector<string> splitString(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

std::vector<bool> fullstr2num(const string& str) {
  std::vector<bool> res;
  for(uint i(0);i<str.size();i++){
    switch (str[i]){
      case 'A':res.push_back(false);res.push_back(false);break;
      case 'C':res.push_back(false);res.push_back(true);break;
      case 'G':res.push_back(true);res.push_back(false);break;
      default:res.push_back(true);res.push_back(true);break;
    }
  }
  return res;
}

std::string fullnum2str(vector<bool> num) {
  string str(num.size()/2, 'N');
  uint j = 0;
  for(uint i(0);i<num.size();i+=2){
    if(num[i]){
      if(num[i+1]){
      	str[j] = 'T';
      }else{
        str[j] = 'G';
      }
    }else{
      if(num[i+1]){
        str[j] = 'C';
      }else{
        str[j] = 'A';
      }
    }
    j++;
  }
  return str;
}

bool isUpperCase(char c) {
	return 'A' <= c and c <= 'Z';
}

int nbCorBases(std::string correctedRead) {
	int n = 0;
	for (unsigned i = 0; i < correctedRead.length(); i++) {
		if ('A' <= correctedRead[i] && correctedRead[i] <= 'Z') {
			n++;
		}
	}

	return n;
}

bool dropRead(std::string correctedRead) {
	return (float) nbCorBases(correctedRead) / correctedRead.length() < 0.1;
}


std::string toUpperCase(std::string& s, int beg, int end) {
	std::string res = s;
	std::locale loc;
	for (int i = beg; i <= end; i++) {
		res[i] = std::toupper(res[i], loc);
	}

	return res;
}

std::string toLowerCase(std::string& s, int beg, int end) {
	std::string res = s;
	std::locale loc;
	for (int i = beg; i <= end; i++) {
		res[i] = std::tolower(res[i], loc);
	}

	return res;
}

std::string trimRead(std::string correctedRead, unsigned merSize) {
	unsigned beg, end, n;
	unsigned i;
	i = 0;
	n = 0;
	while (i < correctedRead.length() and n < merSize) {
		if (isUpperCase(correctedRead[i])) {
			n++;
		} else {
			n = 0;
		}
		i++;
	}
	beg = i - merSize;

	i = correctedRead.length() - 1;
	n = 0;
	while (i >= 0 and n < merSize) {
		if (isUpperCase(correctedRead[i])) {
			n++;
		} else {
			n = 0;
		}
		i--;
	}
	end = i + merSize;

	if (end > beg) {
		return correctedRead.substr(beg, end - beg + 1);
	} else {
		return "";
	}
}

std::vector<std::string> splitRead(std::string name, std::string correctedRead, std::vector<std::pair<unsigned, unsigned>>& pilesPos, unsigned windowSize, unsigned windowOverlap) {
	unsigned i = 0;
	unsigned prevPos = 0;
	unsigned nbUnco = 0;
	std::vector<std::string> result;

	// Skip uncorrected head
	while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
		i++;
	}
	prevPos = i;


	while (i < correctedRead.length()) {
		if (!isUpperCase(correctedRead[i])) {
			nbUnco++;
		} else {
			if (nbUnco >= windowSize) {
				// std::cerr << "nbUnco : " << nbUnco << std::endl;
				result.push_back(correctedRead.substr(prevPos, i - nbUnco - prevPos));
				prevPos = i;
			}
			nbUnco = 0;
		}
		i++;
	}

	// Remove uncorrected tail
	while (i > 0 and !isUpperCase(correctedRead[i])) {
		i--;
	}
	result.push_back(correctedRead.substr(prevPos, i - prevPos));

	return result;
}

void indexReads(robin_hood::unordered_map<std::string, std::vector<bool>>& index, std::string readsFile) {
	std::ifstream f(readsFile);
	std::string header, sequence, seq;
	int nbLines = 0;

	getline(f, header);
	while (header.length() > 0) {
		// Get header
		header.erase(0, 1);
		header = splitString(header, " ")[0];
		
		// Get sequence, watching out for multiline FASTA/FASTQ
		getline(f, seq);
		sequence = seq;
		nbLines = 1;
		getline(f, seq);
		while (seq.length() > 0 and seq[0] != '>' and seq[0] != '+') {
			sequence += seq;
			nbLines++;
			getline(f, seq);
		}

		// Index header/sequence pair
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		index[header] = fullstr2num(sequence);

		// Skip remaining lines if FASTQ
		if (seq[0] == '+') {
			getline(f, seq);
			// while (seq.length() > 0 and seq[0] != '@') {
			for (int i = 1; i < nbLines; i++) {
				getline(f, seq);
			}
			getline(f, seq);
		}

		// Next header has been read, update and loop
		header = seq;
	}
}