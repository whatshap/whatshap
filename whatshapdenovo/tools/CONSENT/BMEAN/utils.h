#ifndef UTILS
#define UTILS



#include <vector>
#include <string>



using namespace std;



typedef uint32_t kmer;



kmer str2num(const string& str);



kmer nuc2int(char c);



kmer nuc2intrc(char c);



string intToString(uint64_t n);



void update(kmer& min, char nuc);



void updateK(kmer& min, char nuc,kmer offsetUpdateKmer);



void updateRC(kmer& min, char nuc);



void updateRCK(kmer& min, char nuc);



string kmer2str(kmer k,uint32_t nuc);



string mutate(string& read, int n);



string rand_seq(const int n);




#endif

