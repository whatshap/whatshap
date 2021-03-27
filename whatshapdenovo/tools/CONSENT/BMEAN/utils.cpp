#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <algorithm>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include "utils.h"



using namespace std;



kmer str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}



string kmer2str(kmer k,uint32_t nuc){
	string result;
	for(uint32_t i(0);i<nuc;++i){
		switch(k%4){
			case 0:  result.push_back('A');break;
			case 1:  result.push_back('C');break;
			case 2:  result.push_back('G');break;
			case 3:  result.push_back('T');break;
		}
		k>>=2;
	}
	reverse(result.begin(),result.end());
	return result;
}


kmer nuc2int(char c){
	switch(c){
		//~ case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return 0;
}



kmer nuc2intrc(char c){
	switch(c){
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		//~ case 'T': return 0;
	}
	return 0;
}




string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



void updateK(kmer& min, char nuc,kmer offsetUpdateKmer){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateKmer;
}



char randNuc(){
	switch (rand()%4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}



string rand_seq(const int n){
	string result;
	for(int32_t i(0);i<n;++i){
		result.push_back(randNuc());
	}
	return result;
}



string mutate(string& read, int n){
	for(int i(0); i<n; ++i){
		int position(rand()%read.size());
		int mutation_type(rand()%3);
		if(mutation_type==0){
		//~ if(mutation_type<4){
			read[position]=randNuc();
		}else if (mutation_type==1){
			read=read.substr(0,position)+read.substr(position+1);
		}else{
			read=read.substr(0,position)+randNuc()+read.substr(position);
		}
	}
	return read;
}





