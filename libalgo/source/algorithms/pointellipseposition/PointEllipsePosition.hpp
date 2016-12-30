// Description: Position and ellipse rotated by the angle (counterclockwise direction)

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


#ifndef PointEllipsePosition_HPP
#define PointEllipsePosition_HPP

#include <math.h>

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"

template <typename T>
unsigned short PointEllipsePosition::getPointEllipsePosition ( const T x, const T y, const T xc, const T yc, const T a, const T b, const T angle )
{
        //Return position of the point and ellipse rotated by the angle -phi
        //		1 = Point inside ellipse,
        //		0 = Point outside ellipse,
        // 		2 = Point on ellipse.

        if ( a == 0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException, bad ellipse parameter.", " Can not compute point and ellipse position", a );
        }

        if ( b == 0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException, bad ellipse parameter.", " Can not compute point and ellipse position", b );
        }

        //Compute test
        const T a_cnd = ( ( x - xc ) * cos ( angle * M_PI / 180 ) - ( y - yc ) * sin ( angle * M_PI / 180 ) );
        const T a_and = ( ( x - xc ) * sin ( angle * M_PI / 180 ) + ( y - yc ) * cos ( angle * M_PI / 180 ) );
        const T t = a_cnd * a_cnd / ( a * a ) + a_and * a_and / ( b * b ) ;

        //Round result
        if ( fabs ( t - 1 ) < POSITION_ROUND_ERROR )
        {
                return 2;
        }

        //Point inside ellipse
        if ( t < 1 )
        {
                return 1;
        }

        //Point outside ellipse
        return 0;
}

#endif
