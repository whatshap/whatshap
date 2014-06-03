#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <cassert>
#include <string.h>

#include <boost/program_options.hpp>

#include "columnreader.h"
#include "dptable.h"

// effectively the max coverage we can handle (due to the size of a byte)
#define MAX_COVERAGE 32

using namespace std;
namespace po = boost::program_options;

// helper functions to output a haplotype (resp., pair of super-reads)
void output_haplotype(auto_ptr<DPTable::haplotype_t> h, ostream & os);
void output_super_reads(auto_ptr<DPTable::read_t> r, ostream & os, bool suppress_undecidable);
string string_rep(Entry::allele_t allele); // get string representation of allele

void usage(const char* name, const po::options_description& options_desc) {
	cerr << "Usage: " << name << " [options] <input.wif>" << endl;
	cerr << endl;
	cerr << "Reads input (in WIF format) and output two resulting \"super reads\" to stdout" << endl;
	cerr << "(also in WIF format)." << endl;
	cerr << endl;
	cerr << options_desc << endl;
	exit(1);
}

// MAIN
int main(int argc, char * const argv[]) {
	string haplotype_filename = "";
	bool all_heterozygous = false;
	bool unweighted = false;
    bool suppress_undecidable = false;     
	
	po::options_description options_desc("Allowed options");
	options_desc.add_options()
		("haplotype,h", po::value<string>(&haplotype_filename)->default_value(""), "Output resulting haplotype to given filename: the output consists of two strings over {0,1,-,X} where \"-\" means \"no data\" and \"X\" means \"data present but couldn't decide\".")
		("all_het,a", po::value<bool>(&all_heterozygous)->zero_tokens(), "Assume all positions to be heterozygous (i.e. fully trust SNP calls).")
		("unweighted,u", po::value<bool>(&unweighted)->zero_tokens(), "Solve unweighted case, i.e. set all weights to 1 regardless of weights given in input.")
        ("suppress_undecidable,s", po::value<bool>(&suppress_undecidable)->zero_tokens(), "Don't output positions for which no decision could be reached (but that were covered by reads).")
	;
	
	if (argc<2) {
		usage(argv[0], options_desc);
	}
	string inputfilename(argv[argc-1]);
	argc -= 1;

	po::variables_map options;
	try {
		po::store(po::parse_command_line(argc, argv, options_desc), options);
		po::notify(options);
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}

	unsigned int max_c = MAX_COVERAGE;

	DPTable dp_table(all_heterozygous);
	cerr << "Computing DP table ... " << endl;
	dp_table.compute_table(new ColumnReader(inputfilename, max_c, unweighted));

  //cout << "optimal score : " << dp_table.get_optimal_score() << endl;

  /*
  cout << "index path ..." << endl << endl;
  auto_ptr<vector<unsigned int> > index_path = dp_table.get_index_path();
  for(size_t i=0; i< index_path->size(); ++i) {
    cout << "index " << i << " : " << index_path->at(i) << endl;
  }
  cout << endl;
  */

  //auto_ptr<vector<bool> > partitioning = dp_table.get_optimal_partitioning();

  /*
  for (size_t i=0; i<partitioning->size(); ++i) {
    cout << "Read " << i << " --> Haplotype " << (partitioning->at(i)?"0":"1") << endl;
  }
  */

	// output super-reads
	auto_ptr<DPTable::read_t> r = dp_table.get_super_reads(new ColumnReader(inputfilename,max_c,unweighted));

	cerr << "Writing super reads ... " << endl;
	output_super_reads(r,cout, suppress_undecidable);

	// output haplotype, if requested
	if (haplotype_filename.size() > 0) {
		auto_ptr<DPTable::haplotype_t> h = dp_table.get_haplotype(new ColumnReader(inputfilename,max_c,unweighted));
		cerr << "Writing haplotype ..." << endl;
		ofstream ofs;
		try { 
			ofs.open(haplotype_filename.c_str(), ios::out);
		} catch(exception & e) { 
			cerr << "Error writing haplotype to \"" << haplotype_filename << "\": " << e.what() << endl;
			return 1;
		}
		output_haplotype(h,ofs);
		ofs.close();
	}

	return 0;
}

/****************************************************************/
// function definitions
/****************************************************************/
    
string string_rep(Entry::allele_t allele) {
  switch (allele) {
  case Entry::MAJOR_ALLELE:
    return "0";
  case Entry::MINOR_ALLELE:
    return "1";
  case Entry::BLANK:
    return "-";
  case Entry::EQUAL_SCORES:
    return "X";
  default:
    assert(false);
  }
}

void output_haplotype(auto_ptr<DPTable::haplotype_t> h, ostream & os) {

  for(size_t i=0; i< h->at(0).size(); ++i) {
    os << string_rep(h->at(0)[i]);
  }  
  os << endl;

  for(size_t i=0; i< h->at(1).size(); ++i) {
    os << string_rep(h->at(1)[i]);
  }
  os << endl;
}

void output_super_reads(auto_ptr<DPTable::read_t> r, ostream & os, bool suppress_undecidable) {
  for (size_t read_nr = 0; read_nr <= 1; ++read_nr) {
    for(size_t i=0; i< r->at(read_nr).size(); ++i) {
      Entry * e = r->at(read_nr)[i];
      if (suppress_undecidable && (e->get_allele_type() != Entry::MAJOR_ALLELE) && (e->get_allele_type() != Entry::MINOR_ALLELE)) {
        continue;
      }
      os << e->get_read_id() << " " << string_rep(e->get_allele_type()) << " " << string_rep(e->get_allele_type()) << " " << e->get_phred_score() << " : ";
    }
    os << "# 0 0 : N N" << endl; // to conform to the .wif format
  }
}
