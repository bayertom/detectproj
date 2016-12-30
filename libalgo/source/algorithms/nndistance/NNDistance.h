// Description: Various NN-distances

// Copyright (c) 2010 - 2013
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


#ifndef NNDistance_H
#define NNDistance_H

#include "libalgo/source/structures/list/Container.h"

//Forward declaration
template <class T>
class KDTree2D;

//Compute several variants of NN-distance
class NNDistance
{
        public:
                template <typename Point1, typename Point2>
                static typename Point1::Type getCrossNearestNeighbourDistance ( const Container <Point1 *> &list1, const Container <Point2 *> &list2 );

                template <typename Point1, typename Point2>
                static typename Point1::Type getAverageNearestNeighbourDistance ( const KDTree2D <Point1> &tree, const Point2 * point, unsigned int k = 10 );

                template <typename Point1, typename Point2>
                static typename Point1::Type compare2DatasetsUsingAverageNearestNeighbourDistance ( const Container <Point1 *> &list1, const Container <Point2 *> &list2, unsigned int k = 10 );
};


#include "NNDistance.hpp"

#endif
