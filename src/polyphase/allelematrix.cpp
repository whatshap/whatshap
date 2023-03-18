#include "allelematrix.h"
#include <algorithm>
#include <set>

using Allele = AlleleMatrix::Allele;
using Position = AlleleMatrix::Position;
using AlleleRow = AlleleMatrix::AlleleRow;
using AlleleItem = AlleleMatrix::AlleleItem;

AlleleMatrix::AlleleMatrix(const std::vector<AlleleRow>& readList, const std::vector<uint32_t>& posList, const std::vector<uint32_t>& idList) :
    m(readList.size()),
    starts(readList.size()),
    ends(readList.size()),
    globalReadIds(readList.size()),
    depths(posList.size()),
    genPos(posList.size()),
    maxAllele(0)
    {
        // create position map
        for (uint32_t i = 0; i < posList.size(); i++) {
            genPos[i] = posList[i];
        }
        std::sort(genPos.begin(), genPos.end());
        for (uint32_t i = 0; i < genPos.size(); i++) {
            posIdx[genPos[i]] = i;
        }
        
        // copy read information
        for (uint32_t i = 0; i < readList.size(); i++) {
            AlleleRow ar;
            std::vector<Position> definedPos;
            for (AlleleItem v: readList[i]) {
                definedPos.push_back(v.first);
                ar[v.first] = v.second;
                if (v.second >= maxAllele) {
                    maxAllele = v.second + 1;
                    for (auto& d: depths)
                        d.resize(maxAllele);
                }
                depths[v.first][v.second] += 1;
            }
            m[i] = ar;
            std::sort(definedPos.begin(), definedPos.end());
            if (definedPos.empty()) {
                starts[i] = -1;
                ends[i] = 0;
            } else {
                starts[i] = definedPos[0];
                ends[i] = definedPos[definedPos.size() - 1];
            }
            globalReadIds[i] = idList[i];
        }
    }

AlleleMatrix::AlleleMatrix(ReadSet* rs) :
    m(rs->size()),
    starts(rs->size()),
    ends(rs->size()),
    globalReadIds(rs->size()),
    depths(rs->get_positions()->size()),
    maxAllele(0)
    {
        // create position map
        for (uint32_t p : *rs->get_positions())
            genPos.push_back(p);
        std::sort(genPos.begin(), genPos.end());
        genPos.shrink_to_fit();
        for (uint32_t i = 0; i < genPos.size(); i++) {
            posIdx[genPos[i]] = i;
        }

        // copy read information
        for (uint32_t i = 0; i < rs->size(); i++) {
            starts[i] = posIdx[rs->get(i)->firstPosition()];
            ends[i] = posIdx[rs->get(i)->lastPosition()];
            globalReadIds[i] = i;
            for (int k = 0; k < rs->get(i)->getVariantCount(); k++) {
                Allele a = (Allele)(rs->get(i)->getAllele(k));
                Position p = this->globalToLocal(rs->get(i)->getPosition(k));
                m[i][p] = a;
                if (a >= maxAllele) {
                    maxAllele = a + 1;
                    for (auto& d: depths)
                        d.resize(maxAllele);
                }
                depths[p][a] += 1;
            }
        }
    }
    
uint64_t AlleleMatrix::size() const {
    return m.size();
}

uint64_t AlleleMatrix::getNumPositions() const {
    return genPos.size();
}

std::vector<Position> AlleleMatrix::getPositions() const {
    return genPos;
}

Allele AlleleMatrix::getMaxNumAllele() const {
    return maxAllele;
}

Allele AlleleMatrix::getAllele(const uint32_t readId, const Position position) const {
    if (m[readId].find(position) == m[readId].end())
        return -1;
    return m[readId].at(position);
}

Allele AlleleMatrix::getAlleleGlobal(const uint32_t readId, const Position genPosition) const {
    if (posIdx.find(genPosition) == posIdx.end())
        return -1;
    return getAllele(readId, posIdx.at(genPosition));
}

std::vector<AlleleItem> AlleleMatrix::getRead(const uint32_t readId) const {
    std::vector<AlleleItem> read(m[readId].begin(), m[readId].end());
    std::sort(read.begin(), read.end(), [] (const AlleleItem& a, const AlleleItem& b) { return a < b; });
    return read;
}

Position AlleleMatrix::getFirstPos(const uint32_t readId) const {
    return starts[readId];
}

Position AlleleMatrix::getLastPos(const uint32_t readId) const {
    return ends[readId];
}

Position AlleleMatrix::getGlobalId(const uint32_t readId) const {
    return globalReadIds[readId];
}

Position AlleleMatrix::globalToLocal(const Position genPosition) const {
    if (posIdx.find(genPosition) == posIdx.end())
        return -1;
    return posIdx.at(genPosition);
}

Position AlleleMatrix::localToGlobal(const Position position) const {
    return genPos[position];
}

std::vector<uint32_t> AlleleMatrix::getAlleleDepths(const Position position) const {
    return depths[position];
}

AlleleMatrix* AlleleMatrix::extractInterval(Position start, Position end, bool removeEmpty) const {
    std::vector<AlleleRow> newReads;
    std::unordered_set<uint32_t> defPos;
    std::vector<uint32_t> idList;
    for (uint32_t i = 0; i < m.size(); i++) {
        if (removeEmpty && (starts[i] >= end || ends[i] < start))
            continue;
        AlleleRow newRead;
        for (auto& entry : m[i]) {
            if (entry.first >= start && entry.first < end) {
                newRead[entry.first - start] = entry.second;
                defPos.insert(localToGlobal(entry.first));
            }
        }
        idList.push_back(globalReadIds[i]);
        newReads.push_back(newRead);
    }
    std::vector<uint32_t> posList(defPos.begin(), defPos.end());
    std::sort(posList.begin(), posList.end());
    return new AlleleMatrix(newReads, posList, idList);
}

AlleleMatrix* AlleleMatrix::extractSubMatrix(const std::vector<Position>& positions,
                                             const std::vector<uint32_t>& readIds,
                                             bool removeEmpty) const {
    std::vector<AlleleRow> newReads;
    std::vector<uint32_t> posList;
    std::vector<uint32_t> idList;
    std::unordered_map<Position, Position> projPos;
    for (uint32_t i = 0; i < positions.size() && positions[i] < this->getNumPositions(); i++) {
        projPos[positions[i]] = i;
        posList.push_back(this->localToGlobal(positions[i]));
    }
    Position start = -1;
    Position end = 0;
    if (positions.size() > 0) {
        start = *std::min_element(positions.begin(), positions.end());
        end = *std::max_element(positions.begin(), positions.end());
    }
    for (uint32_t i : readIds) {
        if (i >= this->size())
            continue;
        if (removeEmpty && (starts[i] >= end || ends[i] < start))
            continue;
        AlleleRow newRead;
        for (auto& entry : m[i])
            if (projPos.find(entry.first) != projPos.end())
                newRead[projPos[entry.first]] = entry.second;
        if (removeEmpty && newRead.empty())
            continue;
        idList.push_back(globalReadIds[i]);
        newReads.push_back(newRead);
    }
    std::sort(posList.begin(), posList.end());
    return new AlleleMatrix(newReads, posList, idList);
}
