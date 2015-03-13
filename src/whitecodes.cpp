#include <bitset>
#include <limits>
#include <cassert>
#include <vector>
#include <math.h>
#include <boost/iterator/iterator_concepts.hpp>

#include "whitecodes.h"

using namespace std;

WhiteCodes::WhiteCodes (int length, int decimalindx)
{
    this->length=length;
    this->k=20;
    this->changed_bit=(int*)malloc(sizeof(int)*30);
    this->c= decimalindx;
    this->iter=1;
    
    for(int i=0;i<30;i++)
     this->changed_bit[i]=-1;
    

    std::string binaryinput="";
    for (int i=0; i < length; ++i)  // assuming a 32 bit int
     binaryinput.append(std::to_string(decimalindx & (1 << i) ? 1 : 0));
    
    int start=0;
    for(int i=0;i<=k;i++)
     ham(start, i, binaryinput, length);

 }
 
WhiteCodes::~WhiteCodes()
{  
  delete changed_bit; 
  
}

void WhiteCodes::ham(int start, int k, std::string bit, int length)
{
  
  if(k==0){     
   all_enum.push_back(bit);     
  return;
  }
  
  if(k>1 && start==length-1)
   return;
  
  if(start>length-1)
   return;
  
  ham(start+1,k,bit,length); 

  if(bit.at(start) =='0'){        
   bit.replace(start, 1, "1");
   ham(start+1, k-1, bit, length);
  }
  else if(bit.at(start)=='1'){   
   bit.replace(start,1,"0");
   ham(start+1,k-1,bit,length);
  }
}

bool WhiteCodes::has_next()
{
  int total_enum=all_enum.size();
  return iter<= total_enum;  
}

WhiteCodes::int_t WhiteCodes::get_next (int** changed_bit)
{
  int result=c;
  
  for(int i=0; i<30; i++)
   (*changed_bit)[i]=this->changed_bit[i];
  
  int total_enum=all_enum.size();
  if(iter<total_enum) {
  for(int i=0;i<30;i++)
    this->changed_bit[i]=-1;
  
   
  int decimal=0;

  for (int j=0;j<length;j++)  
     decimal= (pow(2,j)*(all_enum.at(iter)[j]-'0')) + decimal;
  c=decimal;
 
  
  int j=0;
  for(int i=0;i<length;i++)
  {
   if(all_enum.at(iter)[i]!=all_enum.at(iter-1)[i]){
     this->changed_bit[j]=i;
     j++;
   }
  }
  iter++;
  }else{
    iter=total_enum+1;
  }

  return result;
}
