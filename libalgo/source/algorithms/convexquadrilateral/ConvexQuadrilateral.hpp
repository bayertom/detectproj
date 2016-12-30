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


#ifndef ConvexQuadrilateral_HPP
#define ConvexQuadrilateral_HPP


#include "libalgo/source/exceptions/BadDataException.h"


template <typename Point>
unsigned short ConvexQuadrilateral::isStrictlyConvex ( const Point * p1, const Point * p2, const Point * p3, const Point * p4 )
{
        //Test, if quadrilateral given by 4 points is strictly convex, points are sorted counterclockwise
        //Results: 0   non-convex quadrilateral
        //         1   convex quadrilateral
        //         2   strictly convex quadrilateral

        if ( p1 == NULL || p2 == NULL || p3 == NULL || p4 == NULL )
        {
                throw BadDataException ( "BadDataException: can not test the convexity of the quadrilateral ", " stop computing." );
        }

        //First edge test
        const unsigned short t1 = PointLinePosition::getPointLinePosition2D ( p1, p3, p4 );

        //Strictly convex quadrilateral
        if ( t1 == 2 )
        {
                return 2;
        }

        //Second edge test
        const unsigned short t2 = PointLinePosition::getPointLinePosition2D ( p2, p4, p1 );

        //Strictly convex quadrilateral
        if ( t2 == 2 )
        {
                return 2;
        }

        //Third edge test
        const unsigned short t3 = PointLinePosition::getPointLinePosition2D ( p3, p1, p2 );

        //Strictly convex quadrilateral
        if ( t3 == 2 )
        {
                return 2;
        }

        //Fourth edge test
        const unsigned short t4 = PointLinePosition::getPointLinePosition2D ( p4, p2, p3 );

        //Strictly convex quadrilateral
        if ( t4 == 2 )
        {
                return 2;
        }

        //There are no 3 colinear points, this quadrilateral is convex, but not strictly
        if ( ( t1 + t2 + t3 + t4 == 4 ) || ( t1 + t2 + t3 + t4 == 0 ) )
        {
                return 1;
        }

        //Non-convex quadrilateral
        return 0;
}

#endif

