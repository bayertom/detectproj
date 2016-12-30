// Description: 2D convex hull using Q-hull algorithm

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


#ifndef ConvexHull_H
#define ConvexHull_H

#include "libalgo/source/structures/list/Container.h"

//Forward declaration
class TIndexList;


//Construct convex hull using Q-HULL algorithm
class ConvexHull
{
        public:
                template <typename Point>
                static void getConvexHull ( const Container <Point> &nl, Container <Point, NonDestructable > &ch, const bool print_exception = true, std::ostream * output = &std::cout );

        private:
                template <typename Point>
                static void QHULL ( const unsigned int p1, const unsigned int p2, const Container <Point> &nl, const TIndexList & il , Container <Point, NonDestructable > &hull, const unsigned int n );

};

#include "ConvexHull.hpp"

#endif
