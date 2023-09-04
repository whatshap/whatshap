#ifndef ALLELEMATRIX_H
#define ALLELEMATRIX_H

#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdint>
#include "../readset.h"
#include "../read.h"

/**
 * Stores the allele information of a ReadSet in a more accessible and memory-efficient
 * way. Offers additional functionality, such as allele depths and split functions.
 * Once created, the matrix cannot be modified anymore.
 */
class AlleleMatrix {
    
public:

    typedef int8_t Allele;
    typedef uint32_t Position;
    typedef std::unordered_map<Position, Allele> AlleleRow;
    typedef std::pair<Position, Allele> AlleleItem;

    /**
     * Creates a new allele matrix from a list of reads and list of positions.
     *
     * @param readList the list of reads with allele information.
     * @param posList sorted (!) list of global positions.
     * @param idList list of global read ids.
     */
    AlleleMatrix(const std::vector<AlleleRow>& readList, const std::vector<uint32_t>& posList, const std::vector<uint32_t>& idList);

    /**
     * Creates a new allele matrix from a ReadSet.
     *
     * @param rs the ReadSet to encode.
     */
	AlleleMatrix(ReadSet* rs);
    
    /**
     * Returns the number of reads (rows) inside the matrix.
     */
    uint64_t size() const;

    /**
     * Returns the number of positions (columns) of the matrix.
     */
    uint64_t getNumPositions() const;

    /**
     * Returns the genome positions of the variants as a list.
     */
    std::vector<Position> getPositions() const;

    /**
     * Returns the highest number of different alleles for any variant.
     */
    Allele getMaxNumAllele() const;
    
    /**
     * Returns the allele of the specified read for the given variant. Returns -1 if
     * not defined.
     *
     * @param readId index of read
     * @param position index of variant
     */
    Allele getAllele(const uint32_t readId, const Position position) const;

    /**
     * Returns the allele of the specified read for the given genome position.
     * Returns -1 if not defined.
     *
     * @param readId index of read
     * @param genPosition genome position
     */
    Allele getAlleleGlobal(const uint32_t readId, const Position genPosition) const;

    /**
     * Returns the entire read as dictionary from variant id to allele.
     *
     * @param readId index of read
     */
    std::vector<AlleleItem> getRead(const uint32_t readId) const;

    /**
     * Returns the first defined variant of a read.
     *
     * @param readId index of read
     */
    Position getFirstPos(const uint32_t readId) const;

    /**
     * Returns the last defined variant of a read.
     *
     * @param readId index of read
     */
    Position getLastPos(const uint32_t readId) const;
    
    /**
     * Returns the global id of a read.
     *
     * @param readId index of read
     */
    Position getGlobalId(const uint32_t readId) const;

    /**
     * Converts a genome position to its variant id. Returns -1U if position refers
     * to no variant.
     *
     * @param genPosition genome position
     */
    Position globalToLocal(const Position genPosition) const;

    /**
     * Converts a variant id to its genome position.
     *
     * @param position index of variant
     */
    Position localToGlobal(const Position position) const;
    
    /**
     * Returns the count of each allele for a variant as dictionary from allele to count.
     *
     * @param position index of variant
     */
    std::vector<uint32_t> getAlleleDepths(const Position position) const;

    /**
     * Constructs a submatrix containing only positions of the specified interval. Interval is
     * left-including and right-excluding. Reads without variants can be removed with an
     * additional flag. Global read ids and global positions will be kept, but local indices
     * for reads and positions will change according to the remaining matrix.
     * 
     * @param start first position in the matrix (including)
     * @param end first position NOT in the matrix
     * @param removeEmpty set to true to remove empty reads
     */
    AlleleMatrix* extractInterval(Position start, Position end, bool removeEmpty) const;

    /**
     * Constructs a submatrix containing only the specified positions. Reads without variants
     * can be removed with an additional flag. Global read ids and global positions will be kept,
     * but local indices for reads and positions will change according to the remaining matrix.
     * 
     * @param positions positions for the new matrix
     * @param readIds local read ids to consider
     * @param removeEmpty set to true to remove empty reads
     */
    AlleleMatrix* extractSubMatrix(const std::vector<Position>& positions,
                                   const std::vector<uint32_t>& readIds,
                                   bool removeEmpty) const;

private:
    // data per read
	std::vector<AlleleRow> m;
    std::vector<Position> starts;
    std::vector<Position> ends;
    std::vector<uint32_t> globalReadIds;
    
    // data per position
    std::vector<std::vector<uint32_t>> depths;
    std::vector<Position> genPos;
    std::unordered_map<Position, Position> posIdx;
    Allele maxAllele;
};

#endif

                                                                                                                                                                                      
