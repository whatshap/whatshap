#ifndef PEDMECHEURISTIC_H
#define PEDMECHEURISTIC_H

#include "polyphase/allelematrix.h"
#include "mecheader.h"
#include <vector>
#include <unordered_map>

struct PedSolution {
    Bipartition bp;
    Bipartition bpNew;
    Transmission trans;
    MecScore score;
    MecScore mutationScore;
    RowIndex btRow;
    std::vector<Balance> balances;
    
    PedSolution() :
    bp(0),
    bpNew(0),
    trans(0),
    score(0),
    mutationScore(0),
    btRow(0),
    balances(2, Balance(1, 0)) {}
    
    PedSolution(Bipartition bip, Transmission trans, MecScore score, uint32_t numSamples) :
    bp(bip),
    bpNew(0),
    trans(trans),
    score(score),
    mutationScore(0),
    btRow(0),
    balances(numSamples * 2, Balance(1, 0)) {}
    
    PedSolution(const Bipartition& bip, const Bipartition& bpNew, Transmission trans, MecScore score, MecScore mutationScore, RowIndex btRow, std::vector<Balance>& balances) :
    bp(bip),
    bpNew(bpNew),
    trans(trans),
    score(score),
    mutationScore(mutationScore),
    btRow(btRow),
    balances(balances){
        bp.insert(bp.end(), bpNew.begin(), bpNew.end());
    }
    
    bool operator<(const PedSolution& other) const {
        return score + mutationScore < other.score + mutationScore;
    }
    void finalize() {
        bp.insert(bp.end(), bpNew.begin(), bpNew.end());
    }
};

/**
 * Heuristic colver class for PedMEC problem.
 */
class PedMecHeuristic {

public:

    /**
     * Constructs a solver object for the given input. The read set must not be changed before the
     * solve()-method has been called.
     * @param rs The read set to phase. Sample ids must be zero-indexed and consecutive.
     * @param recombCost Vector of recombination cost per position. Length must match the number of
     *                   positions in read set.
     * @param pedigree Set of all trio relationships and sample genotypes. Must not contain circles
     *                 or children with more than one pair of parents.
     * @param distrust_genotypes If true, phasing may deviate from target genotypes if appropriate
     * @param positions Sorted vector of positions to phase. If undefined, use all positions defined
     *                  in readset. Note that the pedigree stores genotypes by variant index, not by
     *                  genomic variant position.
     * @param rowLimit Number of intermediate solutions to store. Higher means more accurarcy but 
     *                 also more time. Maximum is 65535, default is 1024.
     * @param allowMutations Whether de-novo mutations are allowed in the result. If set to false,
     *                       the input must not contain Mendelian conflicts. Default is true
     * @param verbosity 0 (default) for no output, 1 for output per column, 2 for output per read
     */
    PedMecHeuristic(ReadSet* rs, const std::vector<unsigned int>& recombcost, const Pedigree* pedigree, bool distrust_genotypes, const std::vector<unsigned int>* positions = nullptr, uint32_t rowLimit = 1024, bool allowMutations = true, uint32_t verbosity = 0);
    ~PedMecHeuristic();
    /**
     * Solves the given input instance and stores the results.
     */
    void solve();
    
    /** Returns the internally used scoring scheme for penalizing introduced mutations. */
    std::vector<MecScore> getMutationCost() const;
    
    /** Returns the score of the optimal solution. */
    MecScore getOptScore() const;
    
    /** Returns the computed bipartition of all reads. */
    Bipartition* getOptBipartition() const;
    
    /** Returns the transmission of the optimal solution. Caller inherits ownership of all pointers.*/
    std::vector<Transmission>* getOptTransmission() const;
    
    /** Returns the computed haplotypes. Access order: sample id, haplotype id, position id. */
    std::vector<std::vector<std::vector<Allele>>> getOptHaplotypes() const;
    
    /** Returns the computed haplotypes as a vector of read sets. Caller inherits ownership of all pointers. */
    void getSuperReads(std::vector<ReadSet*>* superReads) const;
    
    /** Returns a boolean representation of the haplotypes of all samples, indicating whether an allele was mutated. Caller inherits ownership of all pointers. */
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>>* getMutations() const;

private:
    // input
    ReadSet* rs;
    std::vector<MecScore> recombCost;
    const Pedigree* pedigree;
    bool distrustGenotypes;
    std::vector<uint32_t>* positions;
    uint32_t rowLimit;
    bool allowMutations;
    uint32_t verbosity;
    
    // computation
    bool solved;
    size_t tmBits;
    uint32_t numSamples;
    std::vector<uint32_t> globalSampleIds;
    std::unordered_map<uint32_t, uint32_t> sampleMap;
    std::vector<Pedigree::triple_entry_t> trios;
    std::vector<std::vector<int>> genotypes; // per sample per position: 0=ref/ref, 1=ref/alt, 2=alt/alt
    std::vector<std::vector<std::pair<Allele, Allele>>> alleles; // per sample per position: ref and alt allele
    std::unordered_map<uint32_t, Position> posMap;
    std::vector<MecScore> mutationCost;
    
    // results
    MecScore optScore;
    Bipartition optBipart;
    std::vector<std::vector<std::vector<Allele>>> optHaps;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> mutations;
    std::vector<Transmission> optTrans;
    
    /**
     * Returns whether two bipartition vectors are equal.
     * @param a first bipartition
     * @param b second bipartition
     */
    bool bpEqual(const Bipartition a, const Bipartition b) const;
    
    /**
     * Copies the score of an existing solution from the previous column to a new solution from
     * the current column, if the score would be improved. Sets backtracking pointer to said old
     * solution if it was an improvement.
     * @param newSol reference to new solution object
     * @param oldSol reference to old solution object
     * @param oldIdx index of old solution in its column
     */
    void updateSolution(PedSolution& newSol, const PedSolution& oldSol, RowIndex oldIdx) const;
    
    /**
     * Returns the recombination cost between two transmissions.
     * @param t1 first transmission
     * @param t2 second transmission
     * @param p position/column (for recombCost vector)
     */
    MecScore getRecombinationCost(const Transmission t1, const Transmission t2, Position p) const;
    
    /**
     * Returns the mutation cost for a set of balance vectors and a given transmission.
     * @param balances The balance vectors representing allele votes per sample
     * @param t transmission vector to use
     * @param p position/column (for mutationCost vector)
     * @param allowFlips if set, uses the allele flip costs to "repair" a mutation if that
     *                   would be cheaper than full mutation cost
     * @param ahead number of following positions/columns to include as well.
     *              Default is 0 (= just column p)
     */
    MecScore getMutationCost(const std::vector<Balance>& balances, const Transmission& t, Position p, bool allowFlips = false, size_t ahead = 0) const;
    
    MecScore getOptPhasing(const std::vector<MecScore>& balances, const Transmission& t, Position p, std::vector<Allele>* optPhasing = nullptr, std::vector<bool>* mutated = nullptr) const;
    
    /**
     * Adds a given allele balance vector to an existing one. Along the process, collects
     * the additionally generated MEC score and returns this as result. All input vectors
     * must have the same length.
     * @param basis balance vector to which the new one is added
     * @param coBasis balance vector for the other partition of the same sample (relevant
     *                if genotype is constrained)
     * @param add added balance
     * @param target target genotypes for the corresponding columns:
     *               0=0/0, 1=0/1, 2=1/1, -1 no constraint)
     */
    MecScore addBalance(Balance& basis, const Balance& coBasis, const Balance& add, const std::vector<int>& target) const;
    
    /**
     * Given a solution from a vector of solutions, generate alternatives with same biparition,
     * but different transmissions, if this would result in better overall solutions (because
     * of mutation).
     * @param sols solution vector to which newly generated solutions are added
     * @param toExt index of solution to extend
     * @param p position/column (for mutation/transmission costs)
     */
    void extendSolutions(std::vector<PedSolution>& sols, uint32_t toExt, Position p) const;
    
    /**
     * Removes solutions from a given vector, such that at most rowLimit many solutions or only
     * optimal solutions are left.
     * @param sols vector of solutions to filter
     */
    void filterSolutions(std::vector<PedSolution>& sols) const;
    
    /** Debug */
    void printColumnInfo(Position p, std::vector<ReadId>& startIndex, std::vector<PedSolution>& col) const;
};

#endif
