// Description: Angle given by two lines

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


#ifndef LineLineAngle_HPP
#define LineLineAngle_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/exceptions/MathInvalidArgumentException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
T LineLineAngle::getLineLineAngle ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, const Point3DCartesian <T> * p4 )
{
        //2 lines angle (0, 180) deg
        return  getLineLineAngle ( p1->getX(), p1->getY(), p2->getX(), p2->getY(), p3->getX(), p3->getY(), p4->getX(), p4->getY() );
}


template <typename T>
T LineLineAngle::getLineLineAngle ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4 )
{
        //2 lines angle (0, 180 deg)
        const T ux = x2 - x1, uy = y2 - y1;
        const T vx = x4 - x3, vy = y4 - y3;

        //Compute denominator
        const T denominator = sqrt ( ( ux * ux + uy * uy ) * ( vx * vx + vy * vy ) );

        //Throw exception, zero devision
        if ( fabs ( denominator ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException ", "can not compute line x line angle:  Identical points, denominator = 0.", denominator );
        }

        //Compute angle
        const T acos_angle = ( ux * vx + uy * vy ) / denominator;

        //Throw exception
        if ( acos_angle > 1 + ARGUMENT_ROUND_ERROR  || acos_angle < -1 - ARGUMENT_ROUND_ERROR )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException, ", "can not compute line x line angle: fabs(acos(argument) > 1 + ARGUMENT_ROUND_ERROR.", acos_angle );
        }

        //Correct longitude
        if ( acos_angle > 1 )
        {
                return ( T ) 0.0;
        }

        //Correct longitude
        if ( acos_angle < -1 )
        {
                return ( T ) - 180.0;
        }

        return ( acos ( acos_angle ) * 180 / M_PI );
}


#endif
