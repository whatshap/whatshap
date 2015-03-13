#ifndef WHITECODES_H
#define WHITECODES_H

#include <iostream>
#include <vector>

class WhiteCodes{
public:
   typedef unsigned int int_t;
  
   WhiteCodes(int length, int decimalindx);
   ~WhiteCodes();

   bool has_next();
   
   int_t get_next(int** changed_bit);
   
   void ham(int start, int k, std::string bit, int length);
   
private:
  int length;
  int *changed_bit;
  int c;
  std::vector<std::string> all_enum;
  int k;
  int iter;
   
};
#endif