// Description: Compute Euclidian distance

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


#ifndef EuclDistance_H
#define EuclDistance_H

#include <cmath>

#include "libalgo/source/structures/point/Point3DCartesian.h"


//Euclidian distance of 2 points
class EuclDistance
{
        public:
                template <typename T>
                static T getEuclDistance ( const T x1, const T y1, const T z1, const T x2, const T y2, const T z2 ) { return sqrt ( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) + ( z2 - z1 ) * ( z2 - z1 ) );}

                template <typename T>
                static T getEuclDistance ( const Point3DCartesian <T> *p1, const Point3DCartesian <T> *p2 ) { return getEuclDistance ( p1->getX(), p1->getY(), p1->getZ(), p2->getX(), p2->getY(), p2->getZ() );}

                template <typename T>
                static T getEuclDistance2D ( const Point3DCartesian <T> *p1, const Point3DCartesian <T> *p2 ) { return getEuclDistance ( p1->getX(), p1->getY(), ( T ) 0,  p2->getX(), p2->getY(), ( T ) 0 ); }

                template <typename T>
                static T getEuclDistance2D ( const T x1, const T y1, const T x2, const T y2 ) { return getEuclDistance ( x1, y1, ( T ) 0, x2, y2, ( T ) 0 ); }
};

#endif

