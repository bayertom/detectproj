// Description: Position of the point and line

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


#ifndef PointLinePosition_HPP
#define PointLinePosition_HPP

#include "libalgo/source/structures/point/Point3DCartesian.h"


template <typename T>
unsigned short PointLinePosition::getPointLinePosition2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 )
{
        /*
         Tests if point P = [xp,yp] is in the left half plane. Return
        		1 = point lies in left half plane,
        		0 = point lies in the right half plane,
         		2 = point lies on the line
         */

        return VectorVectorOrientation::getVectorVectorOrientation2D ( p2->getX() - p1->getX(), p2->getY() - p1->getY(),  p->getX() - p1->getX(), p->getY() - p1->getY() );
}

#endif
