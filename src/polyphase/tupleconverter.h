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

#ifndef TUPLECONVERTER_H
#define TUPLECONVERTER_H

#include "tuple.h"
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>

/**
 * A converter for the haplo-threader algorithm. Takes the global cluster lists of two positions as
 * input and transforms tuples with (position-wise) local cluster ids for one position into a tuples
 * with local cluster ids for the other position, if possible.
 */
class TupleConverter {
    
public:
    /**
     * Constructs the converter based on two global cluster id lists and a ploidy.
     * 
     * @param oldClusters list of old clusters
     * @param newClusters list of new clusters
     * @param ploidy number of elements in a tuple
     */
    TupleConverter(const std::vector<GlobalClusterId>& oldClusters, const std::vector<GlobalClusterId>& newClusters, const uint32_t ploidy);

    /**
     * Copy constructor
     *
     * @param other converter to copy
     */
    TupleConverter(const TupleConverter& other);
    
    /**
     * For a tuple of the new position, returns a tuple for the old position, where the local cluster indices
     * are adjusted accordingly. If any of the used local cluster ids in the input tuple has no correspondence
     * for the old position, an invalid tuple is returned.
     * 
     * @param newTuple Tuple to transform
     */
    ClusterTuple convertNewToOld(const ClusterTuple newTuple) const;
    
    /**
     * For a tuple of the old position, returns a tuple for the new position, where the local cluster indices
     * are adjusted accordingly. If any of the used local cluster ids in the input tuple has no correspondence
     * for the new position, an invalid tuple is returned.
     * 
     * @param oldTuple Tuple to transform
     */
    ClusterTuple convertOldToNew(const ClusterTuple oldTuple) const;
    
    /**
     * For a tuple of the new position, its local ids are permuted in such a way that the similarity to the
     * provided old tuple is maximized. This means that any cluster from the old tuple which also exists for
     * the new position and occurs in the new tuple will keep its position.
     * 
     * @param newTuple Tuple to permute
     * @param oldTuple Tuple against which to align
     */
    ClusterTuple permuteAgainstOld(const ClusterTuple newTuple, const ClusterTuple oldTuple) const;
    
private:
    std::unordered_map<LocalClusterId, LocalClusterId> oldToNew;
    std::unordered_map<LocalClusterId, LocalClusterId> newToOld;
    uint32_t ploidy;
};

#endif // TUPLECONVERTER_H
