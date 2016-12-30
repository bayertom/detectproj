// Description: Compute angle given by 3 points

// Copyright (c) 2010 - 2015
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


#ifndef Angle3Points_HPP
#define Angle3Points_HPP


#include <math.h>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
T Angle3Points::getAngle3Points ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3 )
{
        //Compute angle between ordered set of three point p_left, p_mid, p_right, angle in interval <0, 360>
        if ( ( ( fabs ( p1->getX() - p2->getX() ) < MIN_FLOAT ) && ( fabs ( p1->getY() - p2->getY() ) < MIN_FLOAT ) ) ||
                        ( ( fabs ( p1->getX() - p3->getX() ) < MIN_FLOAT ) && ( fabs ( p1->getY() - p3->getY() ) < MIN_FLOAT ) ) ||
                        ( ( fabs ( p2->getX() - p3->getX() ) < MIN_FLOAT ) && ( fabs ( p2->getY() - p3->getY() ) < MIN_FLOAT ) ) )
        {
                throw BadDataException ( "BadDataException: can not compute angle between 3 points, ", " at least two points are indentical." );
        }

        //Angle difference
        const T angle = ( atan2 ( p2->getY() - p3->getY(), p2->getX() - p3->getX() ) - atan2 ( p2->getY() - p1->getY(), p2->getX() - p1->getX() ) ) * 180 / M_PI;

        return angle < 0 ? angle + 360 : angle;
}


#endif
