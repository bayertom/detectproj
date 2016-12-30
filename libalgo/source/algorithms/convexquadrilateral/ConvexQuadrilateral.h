// Description: Tests, if  quadrilateral is strictly convex

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

#ifndef ConvexQuadrilateral_H
#define ConvexQuadrilateral_H


#include "libalgo/source/algorithms/pointlineposition/PointLinePosition.h"


//Tests, if  quadrilateral is strictly convex
class ConvexQuadrilateral
{
        public:
                template <typename Point>
                static unsigned short isStrictlyConvex ( const Point * p1, const Point * p2, const Point * p3, const Point * p4 );
};

#include "ConvexQuadrilateral.hpp"

#endif
