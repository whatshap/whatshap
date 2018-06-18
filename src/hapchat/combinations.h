/*

  Copyright (C) 2015-2018 Yuri Pirola, Simone Zaccaria

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <stdlib.h>
#include <vector>
#include <bitset>
#include <iostream>
#include <string.h>

#include "binomialcoefficient.h"
#include "basictypes.h"

using namespace std;

class Combinations {

public:

  Combinations()
    : n(0), t(0), count(0), end(false), combination(2, -1), 
      result(), ending_mask(), j(0), x(0), cumulative(false), t_cumulative(0), first(false)
  {
    BinomialCoefficient::binomial_coefficient(32, 32);
  }


  void initialize(const unsigned int n, const unsigned int t)
  {
    this->n = n;
    this->t = t;
    count = -1;
    end = false;
    combination.resize(t + 2);
    fill(combination.begin(), combination.end(), -1);
    result.reset();
    ending_mask.reset();
    j = 0;
    x = 0;
    cumulative = false;
    first = true;
  }


  void initialize_cumulative(const unsigned int n, const unsigned int t)
  {
    this->n = n;
    this->t = 0;
    count = -1;
    end = false;
    combination.resize(t + 2);
    fill(combination.begin(), combination.end(), -1);
    result.reset();
    ending_mask.reset();
    j = 0;
    x = 0;
    cumulative = true;
    t_cumulative = t;
    first = true;
  }


  bool has_next()
  {
    return !end;
  }


  void next()
  {
    if(!end)
      {
        if(first)
          {
            //cout << endl << "QUA 1" << endl;
            first = false;
            initializing_step();
          } else {
          //cout << endl << "QUA 2" << endl;
          basic_step();
        }
        if(cumulative && end && !(t == t_cumulative))
          {
            //cout << endl << "QUA 3" << endl;
            end = false;
            t++;
            combination.resize(t + 2);
            fill(combination.begin(), combination.end(), -1);
            ending_mask.reset();
            j = 0;
            x = 0;
            first = true;
          }
      }
  }


  void get_combination(bitset<MAX_COVERAGE> &result)
  {
    result.reset();
    result |= this->result;
  }


  int get_index()
  {
    return count;
  }


  unsigned int indexof(bitset<MAX_COVERAGE> comb)
  {
    int k = 0;
    int c_k = 0;
    int temp = 0;
    int result = 0;
    //binomial_coefficient(32, 32);
    while(comb.any())
      {
        temp = ffsl(comb.to_ulong());
        c_k += temp;
        k++;
        result += BinomialCoefficient::binomial_coefficient(c_k - 1, k);
        comb>>=(temp);
        //std::cout << comb.to_string() << std::endl;
        //std::cout << "   (" << c_k-1 << ", " << k << ") = " << binomial_coefficient(c_k - 1, k) << "     "; 
      }
    return result;
  }


  unsigned int cumulative_indexof(bitset<MAX_COVERAGE> comb, const unsigned int n_elements)
  {
    unsigned int k = comb.count();
    unsigned int result = indexof(comb);
    for(unsigned int i = 0; i < k; i++)
      {
        result += BinomialCoefficient::binomial_coefficient(n_elements, i);
      }
    return result;
  }


  void combinationof(const unsigned int index, const unsigned int n_elements,
                     const unsigned int k, bitset<MAX_COVERAGE> &result)
  {
    unsigned int iterator = n_elements - 1;
    unsigned int counter = k;
    unsigned int position = 0;
          
    result.reset();

    while(counter > 0)
      {
        if(index >= BinomialCoefficient::binomial_coefficient(iterator, counter) + position)
          {
            position += BinomialCoefficient::binomial_coefficient(iterator, counter);
            result.set(iterator, 1);
            counter--;
          }
        iterator--;
      }
  }


  void cumulative_combinationof(const unsigned int index, const unsigned int n_elements,
                                const unsigned int max_k, bitset<MAX_COVERAGE> &result)
  {
    unsigned int position = BinomialCoefficient::binomial_coefficient(n_elements, 0);
    int counter = max_k - 1;

    for(int i = counter; i > 0; i--)
      {
        position += BinomialCoefficient::binomial_coefficient(n_elements, i);
      }
    while(index < position)
      {
        position -= BinomialCoefficient::binomial_coefficient(n_elements, counter);
        counter--;
      }
    counter++;
    combinationof(index - position, n_elements, counter, result);
  }


  void start_from(bitset<MAX_COVERAGE> comb, const unsigned int n_elements)
  {
    int counter = 0;
    int sum = 0;
    int temp;
    int max_k = comb.count();

    initialize(n_elements, max_k);
    count = indexof(comb);
          
    while(comb.any())
      {
        temp = ffsl(comb.to_ulong());
        sum += temp;
        comb>>=(temp);
        this->combination[counter] = sum - 1;
        result.set(sum - 1, 1);
        counter++;
      }
    combination[counter] = n_elements;
    combination[counter + 1] = 0;


    j = 0;
    while(combination[j + 1] <= j)
      {
        j++;
      }
          
    first = false;

    for(int y = 1; y < max_k + 1; y++)
      {
        ending_mask.set(n - y, 1);
      }
          
    check_end();
  }


  void cumulative_start_from(bitset<MAX_COVERAGE> comb, const unsigned int n_elements, 
                             const unsigned int max_k)
  {
    start_from(comb, n_elements);
    cumulative = true;
    count = cumulative_indexof(comb, n_elements);
    t_cumulative = max_k;

    if(end && !(t == t_cumulative))
      {
        end = false;
        t++;
        combination.resize(t + 2);
        fill(combination.begin(), combination.end(), -1);
        ending_mask.reset();
        j = 0;
        x = 0;
        first = true;
      }
  }


private:
  unsigned int n; //Length of the sequence, the starting set of elements
  unsigned int t; //Length of the combinations
  unsigned int count;
  bool end; 
  vector<int> combination;
  bitset<MAX_COVERAGE> result;
  bitset<MAX_COVERAGE> ending_mask;
  int j;
  int x;
  bool cumulative;
  unsigned int t_cumulative;
  bool first;


  void initializing_step()
  {
    j = 0;
    x = 0;
    for(j = 1; j < (int)t + 1; j++)
      {
        combination[j - 1] = j - 1;
        ending_mask.set(n - j, 1);
      }
    combination[t] = n;
    combination[t + 1] = 0;
    j = t;
    make_combination();
    check_end();
  }


  void basic_step()
  {
    if(j > 0)
      {
        x = j;
        combination[j - 1] = x;
        j = j - 1;
      } else {
      if((combination[0] + 1) < combination[1])
        {
          combination[0] = combination[0] + 1;
        } else {
        j = 2;
        bool flag = false;
        do {
          combination[j - 2] = j - 2;
          x = combination[j - 1] + 1;
          flag = false;
          if(x == combination[j])
            {
              j++;
              flag = true;
            }
        } while(flag);
        //XXX: Check (int) carefully
        if (j > (int)t)
          {
            end = true;
          } else {
          combination[j - 1] = x;
          j--;
        }				
      }
    }
    make_combination();
    check_end();
  }


  void make_combination()
  {
    result.reset();
    for(unsigned int k = 0; k < t; k++)
      {
        result.set(combination[k], 1);		
      }
    count++;
  }


  void check_end()
  {
    end = (result == ending_mask);
  }

};

#endif
