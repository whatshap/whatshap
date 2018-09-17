/*

  Copyright (C) 2017-2018 Marco Dell'Acqua

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

//include c++ libraries
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <ios>

//include Readset/HapChat libraries
#include "../readset.h"
#include "../columniterator.h"
#include "basictypes.h"
#include "../entry.h"

using namespace std;

//Struct iterator for readset

typedef vector<Read*> reads;
typedef vector<reads> blocker;

class HapChatColumnIterator {

private: 

  ColumnIterator* iterator;
  bool end;
  blocker vblock;
  unsigned int blockn;
  vector<unsigned int> min,max;
  ReadSet* readset;

public:

  //standard constructor
  HapChatColumnIterator(ReadSet* read_set){
    readset=read_set;
    iterator=new ColumnIterator(*readset);
    end=false;
    vblock=blocker();
    set_block(readset);
  }


  void set_block(ReadSet* readset){

    Read* read;
    blockn=-1;
    bool overflag;
    unsigned int readn;
    readn=iterator->get_read_count();
    unsigned int minn,maxx;

    for(unsigned int i=0;i<readn;i++){
      read=readset->get(i);
      minn=read->firstPosition();
      maxx=read->lastPosition();

      if(min.size()==0){
	min.push_back(minn);
	max.push_back(maxx);
	vblock.push_back(reads());
	vblock.at(0).push_back(read);
      }
      else{

	for(unsigned int j=0;j<min.size();j++){
	  overflag=minn<min.at(j)&&maxx>max.at(j);
	  
	  if((minn>=min.at(j)&&minn<=max.at(j))||(maxx>=min.at(j)&&maxx<=max.at(j))||overflag){
	    minn=std::min(min.at(j),minn);
	    maxx=std::max(max.at(j),maxx);
	    min.at(j)=minn;
	    max.at(j)=maxx;
	    vblock.at(j).push_back(read);
	    minn=0;
	    break;
	  }
	}

	if(minn!=0){
	  min.push_back(minn);
	  max.push_back(maxx);
	  vblock.push_back(reads());
	  vblock.at(vblock.size()-1).push_back(read);
	}
      }
    }
  }


  bool has_next_block() {

    blockn++;
    
    if(blockn>=vblock.size()) return false;
    ReadSet* readset=new ReadSet();

    for(unsigned int i=0;i<vblock.at(blockn).size();i++){
      readset->add(vblock.at(blockn).at(i));		
    }

    readset->reassignReadIds();
    readset->sort();
    this->iterator=new ColumnIterator(*readset);
    return true;
  }


  //take the current column 
  Column get_column(){

    if(has_next()){
      unique_ptr<vector<const Entry *> > next;
      //if(iterator->it.has_next()) { cout <<"has next" << endl; } // sanity check
			
      next= iterator->get_next();
      vector<const Entry*>* p=next.release();
      Column column;
      
      for(unsigned int i=0;i<p->size();i++){
	column.push_back(*p->at(i));
      }

      return column;
    }

    end=true;
    return Column(0, Entry(-1, Entry::BLANK, 0));
  }


  //return true if there is other column
  bool has_next(){

    return iterator->has_next();
  }


  //set the pointer to the first column
  void reset(){

    iterator->jump_to_column(0);
    end=false;
  }

		
  //print the readid,allele,quality of each entry of the column(for test)
  void print(Column column){
  
    cout << "column: ";

    for(unsigned int i=0; i<column.size();i++){
      cout << column[i].get_read_id()<<","<<column[i].get_allele_type()<<","<<column[i].get_phred_score()<<";";
    }
    cout << endl;
  }

	
  const vector <unsigned int> get_positions(){

    iterator=new ColumnIterator(*readset);
    return *iterator->get_positions();
  }


  unsigned int column_count(){

    return iterator->get_column_count();    
  }


  bool is_ended(){

    return end;
  }

};
