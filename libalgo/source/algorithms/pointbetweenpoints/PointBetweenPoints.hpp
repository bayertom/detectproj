// Description: Test if a point is betwwen 2 points

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


#ifndef PointBetweenPoints_HPP
#define PointBetweenPoints_HPP

#include "../pointlineposition/PointLinePosition.h"

template <typename T>
bool PointBetweenPoints::isPointBetweenPoints ( const T xp, const T yp, const T x1, const T y1, const T x2, const T y2 )
{
        //Test if a point [xp, yp] is between 2 points (belonging to segment [(x1, y1), (x2, y2)] )

        //Segment [(x1, y1), (x2, y2)] is not vertical, compare x coordinates
        if ( x1 != x2 )
        {
                return ( ( x1 <= xp ) && ( xp <= x2 ) ) ||
                       ( ( x1 >= xp ) && ( xp >= x2 ) ) ;
        }

        //Otherwise compare y coordinates
        else
        {
                return ( ( y1 <= yp ) && ( yp <= y2 ) ) ||
                       ( ( y1 >= yp ) && ( yp >= y2 ) ) ;
        }
}

#endif
