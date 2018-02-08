#ifndef GENOTYPE_DP_TABLE
#define GENOTYPE_DP_TABLE

#include <array>
#include <vector>
#include <memory>

#include "columnindexingscheme.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"
#include "pedigree.h"
#include "pedigreepartitions.h"
#include "vector2d.h"
#include "backwardcolumniterator.h"
#include "transitionprobabilitycomputer.h"

class GenotypeDPTable
{
private:

  // stores genotype likelihoods for a given individual and a given position
  struct genotype_likelihood_t {
    // stores likelihoods in this order: 0/0, 1/0, 1/1
    std::vector<long double> likelihoods;
    size_t ind_id;
    size_t position;
    genotype_likelihood_t():likelihoods(3,0.0L),ind_id(0),position(0){}
    genotype_likelihood_t(long double abs, long double het, long double hom, size_t ind_id, size_t pos) :likelihoods(3,0.0L),ind_id(ind_id),position(pos)
    {
      likelihoods[0] = abs;
      likelihoods[1] = het;
      likelihoods[2] = hom;
    }
    // divide likelihoods by given value
    void divide_likelihoods_by(long double& val){
      std::transform(likelihoods.begin(), likelihoods.end(), likelihoods.begin(), std::bind2nd(std::divides<long double>(), val));
    }

    friend inline std::ostream& operator<<(std::ostream& out, const genotype_likelihood_t& g){
      out << g.likelihoods[0] << " " << g.likelihoods[1] << " " << g.likelihoods[2];
      return out;
    }
  };

  // the input sequencing reads
  ReadSet* read_set;
  // stores sample index for each read
  std::vector<unsigned int> read_sources;
  // the recombination cost vector
  const std::vector<unsigned int>& recombcost;
  // the pedigree containing all the individuals
  const Pedigree* pedigree;
  std::vector<PedigreePartitions*> pedigree_partitions;
  // indexing schemes
  std::vector<ColumnIndexingScheme*> indexers;
  // projection_column_table[c] contains the projection column between columns c and c+1
  std::vector<Vector2D<long double>* > forward_projection_column_table;
  std::vector<Vector2D<long double>* > backward_projection_column_table;
  // genotype likelihoods for each individual at each position
  Vector2D<genotype_likelihood_t> genotype_likelihood_table;
  //iterator used to iterate the columns of the input matrix (forward)
  ColumnIterator input_column_iterator;
  // iterator used to iterate the columns of the input matrix (backward)
  BackwardColumnIterator backward_input_column_iterator;
  // stores the transmission probability computers for each column
  std::vector<TransitionProbabilityComputer*> transition_probability_table;
  // scaling parameters
  std::vector<long double> scaling_parameters;

  // helper to pull read ids out of read column
  std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);
  // initializes all members associated with the DP table
  void clear_forward_table();
  void clear_backward_table();
  // forward pass: computes the forward probabilities
  void compute_forward_prob();
  // backward pass: computes the backward probabilities
  void compute_backward_prob();
  // computes the index for each column
  void compute_index();

  // computes column of forward probabilities of given index, assuming previous column was already computed (from left to right)
  void compute_forward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

  // computes column of backward probabilities of given index, assuming previous column was already computed (from right to left)
  void compute_backward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

  // returns the number of bits set
  static size_t popcount(size_t x);

  // given two transmission vectors and their length, compute probability of changing from t1 to t2
  long double compute_transition_prob(size_t t1, size_t t2, size_t length, unsigned int r);

  // used to initialize/clear tables
  template<class T>
  void init(std::vector<T*>& v, size_t size)
  {
    for(size_t i=0; i<v.size(); ++i) {
        if (v[i] != nullptr) {
            delete v[i];
        }
    }
    v.assign(size,nullptr);
  }

public:
  /** Constructor
   * @param read_set   DP table is constructed for the given reads. Ownership is retained by caller.
   *			Pointer must remain valid during the lifetime of this GenotypeDPTable.
   * @param recombcost phred scaled recombination probabilities
   * @param pedigree the pedigree giving individuals and their relationships
   * @param positions positions to work on. If 0, all positions given in the read_set are used.
   * 		      caller retains ownership.
   */
  GenotypeDPTable(ReadSet* read_set, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, const std::vector<unsigned int>* positions = nullptr);
  ~GenotypeDPTable();

  // returns the computed genotype likelihoods for a given individual and a given SNP position
  std::vector<long double> get_genotype_likelihoods(unsigned int individual, unsigned int position);

};
#endif
