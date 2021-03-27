#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <unordered_map>
#include "bmean.h"
#include "utils.h"
#include <time.h>
#include "Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"




using namespace std;
using namespace chrono;



void expe(duration<double>& elapsed_seconds,uint& miss,uint size ){
	vector<string> input;
	string str(rand_seq(size));
	for(uint i(0);i<50;++i){
		string str_mut(str);
		mutate(str_mut,str_mut.size()*10/100);
		input.push_back(str_mut);
	}
	auto start = std::chrono::system_clock::now();
	auto r(MSABMAAC(input,8,10,1,1,100,""));
	auto result = r.first;
	auto  end = std::chrono::system_clock::now();
	elapsed_seconds+= (end - start);

    //~ std::cout  << "elapsed time: " << elapsed_seconds.count() << "s\n";
	//~ cout<<"==========================================================="<<endl;
	//~ cout<<"================THE END===================================="<<endl;
	//~ cout<<"==========================================================="<<endl;




	//~ for(uint32_t iR(0);iR<result.size();++iR){
		//~ cout<<input[iR]<<" splitted in:			";
		//~ for(uint32_t is(0);is<result[iR].size();++is){
			//~ cout<<result[iR][is]<<"	|";
		//~ }
		//~ cout<<endl;
	//~ }
	//~ cout<<str<<endl;
	//~ cout<<result.size()<<endl<<endl;
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	int32_t maskLen=500;



	uint32_t total_size(0);
	for(uint32_t iR(0);iR<result.size();++iR){
		//~ cout<<"-----------------------------------------------------------------------------------------------------------------------------------------"<<endl;
		if(result[iR].size()==0){
			//~ cout<<"WTF de NADINE au FROMAGE"<<endl;
			continue;
		}
		total_size+=result[iR][0].size();
		for(uint32_t is(0);is<result[iR].size();++is){
			//~ cout<<result[iR][is]<<"\n";
			//~ cout<<str<<'\n'<<endl;
			aligner.Align(result[iR][is].c_str(), str.c_str(), str.size(), filter, &alignment, maskLen);
			miss+=alignment.mismatches;
			cout<<"MISSMATCHES "<<alignment.mismatches<<" "<<alignment.ref_begin<<" "<<alignment.ref_end <<" 6gare "<<alignment.cigar_string<<" "<<total_size << endl<<endl;;
			break;
		}
		//~ cout<<endl;
		//~ cout<<total_size<<endl;
		//~ if(result[iR].size()==0){
			//~ cout<<"EMPTY IM SAD"<<endl;cin.get();
		//~ }
	}
}


int main(int argc, char ** argv){
	srand (time(0));
	uint n(0),size(500),iter(200);
	duration<double> elapsed_seconds;
	uint miss(0);
	while(true){
		expe(elapsed_seconds,miss,size);
		if(++n>iter){
			break;
		}
	}
	cout<<"\ncorrected "<<iter*size/1000<<"k nucleotides"<<endl;
	cout  << "Mean elapsed time: " << elapsed_seconds.count()<< "s\n";
	cout<< "Mean missmacthes: "<<(double)miss*10000/n/size<<endl;
	return 0;
}
