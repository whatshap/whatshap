#include "blockreader.h"

bool BlockReader::has_next_nounique() 
{
  if(end) {
    return false;
  } else {
    fragment_block.clear();
    read_positions.clear();
    max_position = -1;
    if(!last_fragment.empty()) {
      add_positions(last_fragment);
      fragment_block.push_back(last_fragment);
    }
    Fragment read;
    string line;

    bool end_block = false;
    while (!end_block) {
      if(!input.eof()) {
        getline(input, line, '\n');
        if(!line.empty()) {
          string_to_fragment(line, read);
          if(read[0].position <= max_position || max_position == -1) {
            add_positions(read);
            fragment_block.push_back(read);
          } else {
            end_block = true;
            already_got = false;
            last_fragment = read;
          }
        }
      } else {
        end_block = true;
        already_got = false;
        end = true;
      }
    }

    return true;
  }
}



bool BlockReader::has_next_unique() 
{
  if(end) {
    return false;
  } else {
    fragment_block.clear();
    read_positions.clear();
    max_position = -1;

    Fragment read;
    string line;

    while (!input.eof()) {
      getline(input, line, '\n');
      if(!line.empty()) {
        string_to_fragment(line, read);        
        add_positions(read);
        fragment_block.push_back(read);
      }
    }

    already_got = false;
    end = true;

    return true;
  }
}




Block BlockReader::get_block() {
  if(!already_got) {
    extract_block();
    already_got = true;
  }

  return block;
}





void BlockReader::string_to_fragment(const string &line, Fragment &read) 
{
  read.clear();
  stringstream sline(line);
  string entry;
  string token;

  int position = -1;
  bool allele = false;
  unsigned int phred = 0;

  bool flag = true;
  while(flag) {
    getline(sline, entry, ':');
    stringstream sentry(entry);

    sentry >> token;

    if(entry.empty() || token.empty() || sline.eof()) {
      cerr << "ERROR: wif input file not well formatted!" << endl;
      exit(EXIT_FAILURE);
    } else if(token.compare("#") != 0) {
      position = atoi(token.c_str());

      sentry >> token;

      sentry >> token;
      if(token.compare("0") == 0) {
        allele = false;
      } else if (token.compare("1") == 0) {
        allele = true;
      } else {
        cerr << "ERROR: found an entry in wif file with an allele not 0 or 1" << endl;
        exit(EXIT_FAILURE);
      }

      sentry >> token;
      phred = atoi(token.c_str());
      if(unweighted) {
        phred = 1;
      }

      read.push_back(EntryRead(position, allele, phred));
    } else {
      flag = false;
      if(read.empty()) {
        cerr << "ERROR: empty read are not allowed in the input wif" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
 
}



void BlockReader::add_positions(const Fragment &read) 
{
  for(Fragment::const_iterator iread = read.begin();
      iread != read.end();
      ++iread) {
    if(read_positions.find((*iread).position) ==  read_positions.end()) {
      max_position = max(max_position, (*iread).position);
      read_positions.insert((*iread).position);
    }
  }
}




//Assumption: the fragments of the input wif are sorted by starting position
void BlockReader::extract_block()
{
  fragment_pointers.clear();
  unsigned int starting_fragment = 0;
  Counter current_cov = 0;
  block.clear();

  vector<int> vread_positions(read_positions.begin(), read_positions.end());
  std::sort(vread_positions.begin(), vread_positions.end());

  for(vector<Fragment>::const_iterator ifblock = fragment_block.begin();
      ifblock != fragment_block.end();
      ++ifblock) {
    fragment_pointers.push_back((*ifblock).begin());
  }
  
  for(vector<int>::const_iterator current_position = vread_positions.begin();
      current_position != vread_positions.end();
      ++current_position) {
    block.push_back(Column());
        
    current_cov = 0;
    for(unsigned int iread = starting_fragment;
        iread < fragment_pointers.size() && ( (fragment_pointers[iread] != fragment_block[iread].begin()) || 
                                              ((*fragment_pointers[iread]).position == *current_position) );
        ++iread) {
            

      if(fragment_pointers[iread] != fragment_block[iread].end()) {
        if(++current_cov > threshold_cov) {
          cerr << "ERROR: coverage threshold excedeed:  "<< current_cov << endl;
          exit(EXIT_FAILURE);
        }
        
        if((*fragment_pointers[iread]).position == *current_position) {
          block.back().push_back(Entry(iread,
                                       ((*fragment_pointers[iread]).allele) ? Entry::MINOR_ALLELE : Entry::MAJOR_ALLELE,
                                       (*fragment_pointers[iread]).phred_score));
          ++fragment_pointers[iread];
        } else {
          if(unweighted) {
            cerr << "ERROR: HapCHAT cannot manage gaps in the unweighted version" << endl;
            exit(EXIT_FAILURE);
          } else {
            block.back().push_back(Entry(iread,
                                         Entry::BLANK,
                                         0));
          }
        }

      }
    }

    while(fragment_pointers[starting_fragment] == fragment_block[starting_fragment].end()) {
      ++starting_fragment;
    }
  }
}





