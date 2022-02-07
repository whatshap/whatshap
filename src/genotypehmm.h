#ifndef GENOTYPE_HMM
#define GENOTYPE_HMM

#include <array>
#include <vector>
#include <memory>

#include "column.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"
#include "pedigree.h"
#include "pedigreepartitions.h"
#include "vector2d.h"
#include "backwardcolumniterator.h"
#include "transitionprobabilitycomputer.h"
#include "emissionprobabilitycomputer.h"

class GenotypeHMM
{
private:

  // stores genotype likelihoods for a given individual and a given position
  struct genotype_likelihood_t {
    // stores likelihoods in this order: 0/0, 1/0, 1/1
    std::vector<long double> likelihoods;
    size_t ind_id;
    size_t position;
    genotype_likelihood_t() : likelihoods() {}
    genotype_likelihood_t(unsigned int n_genotypes):likelihoods(n_genotypes ,0.0L),ind_id(0),position(0){}
   
    // divide likelihoods by given value
    void divide_likelihoods_by(long double& val){
      std::transform(likelihoods.begin(), likelihoods.end(), likelihoods.begin(), std::bind2nd(std::divides<long double>(), val));
    }

    friend inline std::ostream& operator<<(std::ostream& out, const genotype_likelihood_t& g){
      out << g.likelihoods[0] << " " << g.likelihoods[1] << " " << g.likelihoods[2];
      return out;
    }
  };

  // number of reference samples
  unsigned int n_references;
  
  // what allele each reference sample in the vcf file has {{<alleles in position 1>},{<alleles in position 2>},{<alleles in position 3>},....}
  const std::vector<std::vector<int> >* allele_references;
  
  // position of the variants
  const std::vector<unsigned int>* variant_positions;
  
  // number of alleles per variant position
  const std::vector<unsigned int>* variant_n_allele_positions;
  
  // the input sequencing reads
  ReadSet* read_set;
  
  // stores sample index for each read
  std::vector<unsigned int> read_sources;
  
  // the recombination cost vector
  const std::vector<float>& recombcost;
  
  // the pedigree containing all the individuals
  const Pedigree* pedigree;
  // indexing schemes
  std::vector<Column*> hmm_columns;
  
  // projection_column_table[c] contains the projection column between columns c and c+1
  std::vector<std::vector<long double>* > forward_pass_column_table;
  std::vector<std::vector<long double>* > backward_pass_column_table;

  std::vector<std::vector<unsigned int>* > active_reads;
  
  // genotype likelihoods for each individual at each position
  Vector2D<genotype_likelihood_t> genotype_likelihood_table;
  
  //iterator used to iterate the columns of the input matrix (forward)
  ColumnIterator input_column_iterator;
  
  // iterator used to iterate the columns of the input matrix (backward)
  BackwardColumnIterator backward_input_column_iterator;
  
  // stores the transmission probability computers for each column. object at index i contains the probability computer between index i and i+1.
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

  // updates the emission probabilities
  void update_emission_probability(Vector2D<long double>* em_prob, const int bit_changed, const ColumnIndexingIterator& iterator, std::vector<const Entry *>& entries);

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
  GenotypeHMM(ReadSet* read_set, const std::vector<float>& recombcost, const Pedigree* pedigree, const unsigned int& n_references, const std::vector<unsigned int>* positions = nullptr, const std::vector<unsigned int>* n_allele_positions = nullptr, const std::vector<std::vector<int> >* allele_references = nullptr);
  ~GenotypeHMM();

  // returns the computed genotype likelihoods for a given individual and a given SNP position
  std::vector<long double> get_genotype_likelihoods(unsigned int individual, unsigned int position);

};
#endif
