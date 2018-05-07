#include "basic_types.h"
#include <limits>

const Cost::cost_t Cost::infinity_(std::numeric_limits<Cost::cost_t>::max());
const Cost Cost::INFTY(Cost::infinity_);


std::ostream& operator<<(std::ostream& out, const Cost& c) {
  if (c == Cost::INFTY)
    out << "INFINITY";
  else
    out << c.cost_;
  return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<bool>& v) {
  for (std::vector<bool>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << (*it? '1' : '0');
  }
  return out;
}


std::ostream& operator<<(std::ostream& out, const std::vector<char>& v) {
  for (std::vector<char>::const_iterator it= v.begin(); it != v.end(); ++it) {
    out << *it;
  }
  return out;
}



std::ostream& operator<<(std::ostream& out, const options_t& options) {
  const std::string SEP("\n");
  out
    << "Initialized? " << (options.options_initialized?"True":"False") << SEP
    << "Input filename: '" << options.input_filename << '\'' << SEP
    << "Haplotype filename: '" << options.haplotype_filename << '\'' << SEP
    << "Discard weights? " << (options.unweighted?"True":"False") << SEP
    << "Mask ambiguous positions? " << (options.no_xs?"False":"True") << SEP
    << "all-heterozygous assumption?" << (options.all_heterozygous?"True":"False") << SEP
    << "read input as unique block?" << (options.unique?"True":"False") << SEP
    << "Error rate: " << options.error_rate << SEP
    << "Alpha: " << options.alpha;
  return out;
}

#include <getopt.h>
#include <sstream>

options_t parse_arguments(int argc, char** argv) {

  options_t ret;

  // options description
  std::ostringstream oss;

  oss
    << "Program options:" << std::endl

    << "  -h [ --help ]" << std::string(4,'\t') 
    << "produce (this) help message" << std::endl

    << "  -i [ --input ] arg" << std::string(3,'\t')
    << "file containing the input reads (in WIF" << std::endl
    << std::string(5,'\t') << "format)" << std::endl

    << "  -o [ --haplotypes ] arg" << std::string(2,'\t')
    << "file where the computed haplotypes will" << std::endl
    << std::string(5,'\t') << "be written to" << std::endl

    << "  -u [ --discard-weights ]" << std::string(2,'\t')
    << "discard weights" << std::endl

    << "  -x [ --no-ambiguous ]" << std::string(3,'\t')
    << "do not mark ambiguous positions with Xs" << std::endl

    << "  -A [ --all-heterozygous ]" << std::string(2,'\t')
    << "all-heterozygous assumption" << std::endl

    << "  -U [ --unique ]" << std::string(3,'\t')
    << "input as unique block" << std::endl

    << "  -e [ --error-rate ] arg (="
    << ret.error_rate << ")" << std::string(1,'\t')
    << "read error rate" << std::endl

    << "  -a [ --alpha ] arg (="
    << ret.alpha << ")" << std::string(2,'\t')
    << "significance (smaller is better)" << std::endl

    << "  -b [ --balancing ] arg (="
    << ret.balance_cov << ")" << std::string(2,'\t')
    << "maximum coverage before balancing" << std::endl
    << std::string(5,'\t') << "is activated" << std::endl

    << "  -r [ --balance-ratio ] arg (="
    << ret.balance_ratio << ")" << std::string(1,'\t')
    << "balance ratio (larger is stricter)" << std::endl;

  std::string opts_desc = oss.str();

  // loop for obtaining the options
  int opt;
  bool sane = true;
  std::string err;

  while(1) {

    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"input", required_argument, 0, 'i'},
      {"haplotypes", required_argument, 0, 'o'},
      {"discard-weights", no_argument, 0, 'u'},
      {"no-ambiguous", no_argument, 0, 'x'},
      {"all-heterozygous", no_argument, 0, 'A'},
      {"unique", no_argument, 0, 'U'},
      {"error-rate", required_argument, 0, 'e'},
      {"alpha", required_argument, 0, 'a'},
      {"balance-cov", required_argument, 0, 'b'},
      {"balance-ratio", required_argument, 0, 'r'},
      {0, 0, 0, 0}
    };

    // get an option
    int option_index = 0;
    opt = getopt_long(argc, argv, "hi:o:uxAUe:a:b:r:", long_options, &option_index);

    if(opt == -1) // end of options
      break;

    switch(opt)
      {

      case 'h' :
	std::cout << opts_desc << std::endl;
	return ret;
      case 'i' :
	ret.input_filename = optarg;
	break;
      case 'o' :
	ret.haplotype_filename = optarg;
	break;
      case 'u' :
	ret.unweighted = true;
	break;
      case 'x' :
	ret.no_xs = true;
	break;
      case 'A' :
	ret.all_heterozygous = true;
	break;
      case 'U' :
	ret.unique = true;
	break;
      case 'e' :
	ret.error_rate = atof(optarg);
	break;
      case 'a' :
	ret.alpha = atof(optarg);
	break;
      case 'b' :
	ret.balancing = true;
	ret.balance_cov = atof(optarg);
	break;
      case 'r' :
	ret.balance_ratio = atof(optarg);
	break;
      default :
	sane = false;
	err = "unrecognized option";
      }
  }

  // sanity check on required options and ranges
  if(ret.input_filename == "") {
    sane = false;
    err = "the option '--input' is required but missing";
  }
  if(ret.haplotype_filename == "") {
    sane = false;
    err = "the option '--haplotypes' is required but missing";
  }
  if((ret.error_rate < 0.0) || (ret.error_rate > 1.0)) {
    sane = false;
    err = "error-rate must be a value between 0.0 and 1.0";
  }
  if((ret.alpha < 0.0) || (ret.alpha > 1.0)) {
    sane = false;
    err = "alpha must be a value between 0.0 and 1.0";
  }
  if((ret.balancing == true) && (ret.all_heterozygous == false)) {
    sane = false;
    err = "the option '--all-heterozygous' is required when balancing";
  }
  if((ret.balance_ratio < 0.0) || (ret.balance_ratio > 0.5)) {
      sane = false;
      err = "balance ratio must be a value between 0.0 and 0.5";
  }

  if(!sane) {
    std::cout << "ERROR while parsing the program options: ";
    std::cout << err << std::endl;
    std::cout << opts_desc << std::endl;
  }
  else
    ret.options_initialized = true;

  return ret;
}
