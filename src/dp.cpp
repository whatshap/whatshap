#include <iostream>

#include "read.h"
#include "readset.h"
#include "dptable.h"

using namespace std;

void output_vector(const auto_ptr<vector<bool> > v) {
  for(int j=v->size()-1; j>= 0; --j) {
    cout << v->at(j) << " ";
  }
}

int main(int argc, char * const argv[]) {

  // F from pWhatsHap paper

  Read * f1 = new Read("f_1",50);
  f1->addVariant(0,'C',0,5);

  Read * f2 = new Read("f_2",50);
  f2->addVariant(0,'G',1,3);
  f2->addVariant(1,'C',0,2);

  Read * f3 = new Read("f_3",50);
  f3->addVariant(0,'G',1,6);
  f3->addVariant(1,'G',1,1);

  Read * f4 = new Read("f_4",50);
  f4->addVariant(1,'C',0,2);

  ReadSet * input = new ReadSet();
  input->add(f1);
  input->add(f2);
  input->add(f3);
  input->add(f4);
  input->finalize(); // reads cannot be modified now

  cout << "input : " << endl << input->toString() << endl;

  // initialize dp table and phase

  DPTable dp_table(input);

  cout << "score : " << dp_table.get_optimal_score() << endl;

  // output the parititioning

  auto_ptr<vector<bool> > partitioning = dp_table.get_optimal_partitioning();
  cout << endl << "optimal partitioning :" << endl;
  output_vector(partitioning);
  cout << endl;

  return 0;
}
