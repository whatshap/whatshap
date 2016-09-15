#ifndef PEDIGREE_DP_TABLE_H
#define PEDIGREE_DP_TABLE_H

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

typedef struct index_and_inheritance_t {
	unsigned int index;
	unsigned int inheritance_value;
	index_and_inheritance_t()  : index(0), inheritance_value(0) {};
} index_and_inheritance_t;

class PedigreeDPTable {
private:
	ReadSet* read_set;
	// stores the sample index for each read
	std::vector<unsigned int> read_sources;
	const std::vector<unsigned int>& recombcost;
	const Pedigree* pedigree;
	bool distrust_genotypes;
	std::vector<PedigreePartitions*> pedigree_partitions;
	// vector of indexingschemes
	std::vector<ColumnIndexingScheme*> indexers;
	// optimal score and its index in the rightmost DP table column
	unsigned int optimal_score;
	unsigned int optimal_score_index;
	unsigned int optimal_transmission_value;
	// transmission value preceeding the optimal one (in the column before / in the backtrace)
	unsigned int previous_transmission_value;
	// projection_column_table[c] contains the projection column "between" columns c and c+1
	std::vector<Vector2D<unsigned int>* > projection_column_table;
	// index_backtrace_table[c][i][t] indicates the index (=bipartition) in column c from which the
	// i-th entry in the FORWARD projection of column c comes from, assuming a transmission value of t
	std::vector<Vector2D<unsigned int>* > index_backtrace_table;
	// let x := index_backtrace_table[c][i][t] and dp[x][t] the corresponding DP entry
	// and j be the BACKWARD projection of x.
	// Then t' = transmission_backtrace_table[c][i][t] is the transmission index (from {0,1,2,3})
	// that gave rise to dp[x][t].
	std::vector<Vector2D<unsigned int>* > transmission_backtrace_table;
	ColumnIterator input_column_iterator;
	// optimal path obtained from backtrace
	std::vector<index_and_inheritance_t> index_path;

	// helper function to pull read ids out of read column
	std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);

	/** Initializes/clears all member variables associated with the DP table, i.e. indexers, index_backtrace_table,
	 *  transmission_backtrace_table, optimal_score, optimal_score_index, optimal_transmission_value, and previous_transmission_value. */
	void clear_table();
	void compute_table();
	/** Computes the DP column at the given index, assuming that the previous column
	 *  has already been computed. */
	void compute_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

	/** Returns the number of set bits. */
	static size_t popcount(size_t x);

	template <class T>
	void init(std::vector<T*>& v, size_t size) {
		for(size_t i=0; i<v.size(); ++i) {
			if (v[i] != nullptr) {
				delete v[i];
			}
		}
		v.assign(size, nullptr);
	}

public:
	/** Constructor.
	 *  @param read_set DP table is constructed for the contained reads. Ownership is retained
	 *                  by caller. Pointer must remain valid during the lifetime of this PedigreeDPTable.
	 *  @param distrust_genotypes If true, then the genotypes may be changed at costs given as genotype likelihoods
	 *                            (in the given pedigree object).
	 *  @param positions Positions to work on. If 0, then all positions given in read_set will be used. Caller retains
	 *                   ownership.
	 */
	PedigreeDPTable(ReadSet* read_set, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const std::vector<unsigned int>* positions = nullptr);
 
	~PedigreeDPTable();

	unsigned int get_optimal_score();

	/** Computes optimal haplotypes and adds them (in the form of "super reads") to 
	 *  the given read_set.
	 *
	 *   @param output_read_set Must have as many entries as there are individuals. The haplotypes for individual
	 *                          with index i in the pedigree (given at construction time) are added to output_read_set->at(i).
	 */
	void get_super_reads(std::vector<ReadSet*>* output_read_set, std::vector<unsigned int>* transmission_vector);

	/** Performs a backtrace through the DP table and returns optimal partitioning of the reads.
	 *  Pointer ownership is transferred to caller. */
	std::vector<bool>* get_optimal_partitioning();
};

#endif
