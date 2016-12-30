// Description: Distance point-loine

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


#ifndef PointLineDistance_HPP
#define PointLineDistance_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
T PointLineDistance::getPointLineDistance2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 )
{
        //Compute distance point and line
        return getPointLineDistance2D ( p->getX(), p->getY(), p1->getX(), p1->getY(), p2->getX(), p2->getY() );
}


template <typename T>
T PointLineDistance::getPointLineDistance2D ( const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2 )
{
        //Compute unsigned distance point and line
        return fabs ( getPointLineDistance2DSigned ( xa, ya, x1, y1, x2, y2 ) );
}



template <typename T>
T PointLineDistance::getPointLineDistance2DSigned ( const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2 )
{
        //Compute signed distance point and line
        //         If distance:
        //		> 0 = point lies in left half plane,
        //		< 0 = point lies in the right half plane,
        //		= 0 = point lies on the line,

        T denominator = sqrt ( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) );

        //Throw exception
        if ( fabs ( denominator ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute distance point - line, ", "end points of the line are identical.", denominator );
        }

        return ( xa * ( y1 - y2 ) + x1 * ( y2 - ya ) + x2 * ( ya - y1 ) ) / denominator ;
}


template <typename T>
T PointLineDistance::getPointLineSegmentDistance2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 )
{
        //Compute distance point and line segment
        const T nx = p1->getY() - p2->getY(), ny = p2->getX() - p1->getX();

        //Point p3 is on normal n1: given by p1 and perpendicular to (p1, p2)
        const Point3DCartesian <T> p3 ( p1->getX() + nx,  p1->getY() + ny );

        //Point p4 is on normal n2: given by p2 and perpendicular to (p1, p2)
        const Point3DCartesian <T> p4 ( p2->getX() + nx,  p2->getY() + ny );

        //Position of the point according to the both normals using signed distance
        const T dist1 = getPointLineDistance2DSigned ( p->getX(), p->getY(), p1->getX(), p1->getY(), p3.getX(), p3.getY() );
        const T dist2 = getPointLineDistance2DSigned ( p->getX(), p->getY(), p2->getX(), p2->getY(), p4.getX(), p4.getY() );

        //Point is between both normals n1 and n2
        if ( ( dist1 < 0 ) && ( dist2 > 0 ) || ( dist1 > 0 ) && ( dist2 < 0 ) )
        {
                return getPointLineDistance2D ( p->getX(), p->getY(), p1->getX(), p1->getY(), p2->getX(), p2->getY() );
        }

        //Point is left to the n1
        if ( dist1 >= 0 )
        {
                return EuclDistance::getEuclDistance2D ( p, p1 );
        }

        //Point is right to n2
        return EuclDistance::getEuclDistance2D ( p, p2 );
}

#endif
