//include c++ libraries
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <ios>

//include Readset/HapChat libraries
#include "readset.h"
#include "columniterator.h"
#include "basic_types.h"
#include "entry.h"




using namespace std;

//Struct iterator for readset
typedef struct Iterator{
	ColumnIterator it;
	Iterator(ReadSet* readset):it(*readset,nullptr){};
	Iterator():it(*new ReadSet(),nullptr){};

}Iterator;
typedef vector<Read*> reads;
typedef vector<reads> blocker;
class HapCHATcore {

	private: 
		Iterator* iterator;
		bool end,unique;
		blocker vblock;
		unsigned int blockn;
		vector<unsigned int> min,max;
	public:
		//standard constructor
		HapCHATcore(ReadSet* read_set){
		read_set->reassignReadIds();
		this->iterator=new Iterator(read_set);
		end=false;
		vblock=blocker();
		};
		HapCHATcore(string filename,bool que){
		 unique=que;
		 ifstream input;
		 int position,allele;
		 unsigned int phred,nread=0;
		 ReadSet* readset =new ReadSet();
		 Read* read;
		 bool flag;
		 stringstream sline;
  		string entry;
  		string token;
		 try {
     			 input.open(filename, ios::in);
    		 } catch(exception & e) { 
      			cerr << "ERROR: failing opening the input file: " << filename << "\": " << e.what() << endl;
      			exit(EXIT_FAILURE);
    		 }

    if(!input.is_open()) {
    		  cerr << "ERROR: failing opening the input file: " << filename << endl;
    	  exit(EXIT_FAILURE);
    	}
    	string line;
    	while(!input.eof()){
    	getline(input, line, '\n');
    	read=new Read(to_string(nread),0,0,0);
    	nread++;
    	if(!line.empty()){
    			sline=stringstream(line);
    			flag=true;
    			while(flag){
    			getline(sline,entry,':');
    			stringstream sentry(entry);
    			sentry >> token;
    			if(!token.compare("#") == 0){
    			position=atoi(token.c_str());
    			sentry >> token;
    			sentry >> token;
					allele = atoi(token.c_str());    			
    			sentry >> token;
    			phred = atoi(token.c_str());
    			read->addVariant(position,allele,phred);
    			}else{flag=false;}
    			}
    		readset->add(read);    		
    	}   	
    	}
    	
    	end=false;
    	readset->reassignReadIds();
    		readset->sort();
    	this->iterator=new Iterator(readset);
    		
    	if(unique){
			
		
			}else{
			blockn=-1;
			bool overflag;
			unsigned int readn;
			readn=iterator->it.get_read_count();
			unsigned int minn,maxx;
			for(unsigned int i=0;i<readn;i++){
			 read=readset->get(i);
			 minn=read->firstPosition();
			 maxx=read->lastPosition();
			 if(min.size()==0){
			 	min.push_back(minn);
			 	max.push_back(maxx);
			 	vblock.push_back(reads());
			 	}
			 else{
			 	for(unsigned int j=0;j<min.size();j++){
			 	overflag=minn<min.at(j)&&maxx>max.at(j);
			 		if((minn>=min.at(j)&&minn<=max.at(j))||(maxx>=min.at(j)&&maxx<max.at(j))||overflag){
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
			 
			 };
			 
			 
			}
			}
					};
		bool hasNextBlock(){
		if(unique){
		unique=false;
		return true;
		}
		blockn++;
		if(blockn>=vblock.size()) return false;
		ReadSet* readset=new ReadSet();
		for(unsigned int i=0;i<vblock.at(blockn).size();i++){
			readset->add(vblock.at(blockn).at(i));		
		}
		readset->reassignReadIds();
    readset->sort();
    this->iterator=new Iterator(readset);
		return true;
		}
		//take the current column 
		Column getColumn(){
			if(hasNext()){
				unique_ptr<vector<const Entry *> > next;
				//if(iterator->it.has_next()) { cout <<"has next" << endl; } // sanity check
			
				next= iterator->it.get_next();
				vector<const Entry*>* p=next.release();
				Column column;
		
				for(unsigned int i=0;i<p->size();i++){
					column.push_back(*p->at(i));
				}
				return column;
			}
				end=true;
    		return Column(0, Entry(-1, Entry::BLANK, 0));
		};
	
		//return true if there is other column
		bool hasNext(){
			return iterator->it.has_next();
		}
		//set the pointer to the first column
		void reset(){
			iterator->it.jump_to_column(0);
			end=false;
		};
		
		//print the readid,allele,quality of each entry of the column(for test)
		void print(Column column){
			cout << "column: ";
			for(unsigned int i=0; i<column.size();i++){
			cout << column[i].get_read_id()<<","<<column[i].get_allele_type()<<","<<column[i].get_phred_score()<<";";
			};
			cout << endl;
		};
		
		unsigned int columnCount(){
			return iterator->it.get_column_count();

		};
		bool isEnded(){
			return end;
		}
};
