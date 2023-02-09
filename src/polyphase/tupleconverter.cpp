/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2021  <copyright holder> <email>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tupleconverter.h"

TupleConverter::TupleConverter(const std::vector<GlobalClusterId>& oldClusters, const std::vector<GlobalClusterId>& newClusters, uint32_t ploidy) {
    std::unordered_map<GlobalClusterId, LocalClusterId> idMap;
    for (LocalClusterId c = 0; c < newClusters.size(); c++)
        idMap[newClusters[c]] = c;
    for (LocalClusterId c = 0; c < oldClusters.size(); c++) {
        GlobalClusterId g = oldClusters[c];
        if (idMap.find(g) != idMap.end())
            oldToNew[c] = idMap[g];
    }
    idMap.clear();
    for (LocalClusterId c = 0; c < oldClusters.size(); c++)
        idMap[oldClusters[c]] = c;
    for (LocalClusterId c = 0; c < newClusters.size(); c++) {
        GlobalClusterId g = newClusters[c];
        if (idMap.find(g) != idMap.end())
            newToOld[c] = idMap[g];
    }
    this->ploidy = ploidy;
}

TupleConverter::TupleConverter(const TupleConverter& other) :
    oldToNew(other.oldToNew),
    newToOld(other.newToOld) {}

ClusterTuple TupleConverter::convertNewToOld(const ClusterTuple newTuple) const {
    std::vector<LocalClusterId> v;
    for (uint32_t i = 0; i < ploidy; i++) {
        LocalClusterId c = newTuple.get(i);
        if (newToOld.find(c) == newToOld.end())
            return ClusterTuple::INVALID_TUPLE;
        else
            v.push_back(newToOld.at(c));
    }
    return ClusterTuple(v);
}

ClusterTuple TupleConverter::convertOldToNew(const ClusterTuple oldTuple) const {
    std::vector<LocalClusterId> v;
    for (uint32_t i = 0; i < ploidy; i++) {
        LocalClusterId c = oldTuple.get(i);
        if (oldToNew.find(c) == oldToNew.end())
            return ClusterTuple::INVALID_TUPLE;
        else
            v.push_back(oldToNew.at(c));
    }
    return ClusterTuple(v);
}

ClusterTuple TupleConverter::permuteAgainstOld(const ClusterTuple newTuple, const ClusterTuple oldTuple) const {
    // copy new cluster ids into v
    std::vector<int32_t> v;
    for (uint32_t i = 0; i < ploidy; i++) {
        v.push_back(newTuple.get(i));
    }
    
    // allocate ploidy-sized u for construction of result
    std::vector<LocalClusterId> u(ploidy, 0);
    std::vector<uint32_t> resOld;
    for (uint32_t i = 0; i < ploidy; i++) {
        LocalClusterId c = oldTuple.get(i);
        // if old id exists on new position and it occurs in v -> write into u, else -> add to residual indices
        if (oldToNew.find(c) != oldToNew.end()) {
            int32_t d = oldToNew.at(c);
            for (uint32_t j = 0; j < ploidy; j++) {
                // entries in v are invalidated to account for duplicate ids
                if (v[j] == d) {
                    u[i] = d;
                    v[j] = -1;
                    d = -1;
                    break;
                }
            }
            if (d >= 0) {
                resOld.push_back(i);
            }
        } else {
            resOld.push_back(i);
        }
    }
    
    // iterate over non-taken entries of v and write them into residual indices in u
    uint32_t resIdx = 0;
    for (uint32_t i = 0; i < ploidy; i++) {
        if (v[i] >= 0) {
            u[resOld[resIdx]] = v[i];
            resIdx++;
        }
    }
    
    return ClusterTuple(u);
}
