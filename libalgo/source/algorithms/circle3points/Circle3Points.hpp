// Description: Parameters of the circle given by three points

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


#ifndef Circle3Points_HPP
#define Circle3Points_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"
#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
void Circle3Points::getCentreAndDiameterCircle ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, T & xc, T & yc, T & r )
{
        //Get inscribed circle center and center diameter
        const T xt = ( p1->getX() + p2->getX() + p3->getX() ) / 3.0;
        const T yt = ( p1->getY() + p2->getY() + p3->getY() ) / 3.0;

        //Reduce coordinates to keep the accuracy for triangles of the inapproriate shapes
        const T x1r = p1->getX() - xt, y1r = p1->getY() - yt;
        const T x2r = p2->getX() - xt, y2r = p2->getY() - yt;
        const T x3r = p3->getX() - xt, y3r = p3->getY() - yt;

        //Compute k1 - k12
        const T k1 = x1r * x1r + y1r * y1r;
        const T k2 = x2r * x2r + y2r * y2r;
        const T k3 = x3r * x3r + y3r * y3r;
        const T k4 = y1r - y2r;
        const T k5 = y1r - y3r;
        const T k6 = y2r - y3r;
        const T k7 = x1r - x2r;
        const T k8 = x1r - x3r;
        const T k9 = x2r - x3r;
        const T k10 = x1r * x1r;
        const T k11 = x2r * x2r;
        const T k12 = x3r * x3r;

        //Compute denominator
        const T denom1 = x3r * ( -k4 ) + x2r * k5 + x1r * ( -k6 );

        //Throw exception
        if ( fabs ( denom1 ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute mid point of the circle. ", "Circle given by 3 points, 2 points are identical.", denom1 );
        }

        const T denom2 = y1r * ( -k9 ) + y2r * k8 + y3r * ( -k7 );

        //Throw exception
        if ( fabs ( denom2 ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute mid point of the circle.", "Circle given by 3 points, 2 points are identical.", denom2 );
        }

        //Center of the circle
        xc = 0.5 * ( k12 * ( -k4 ) + k11 * k5 - ( k10 + k4 * k5 ) * k6 ) / denom1 + xt;
        yc = 0.5 * ( k1 * ( -k9 ) + k2 * k8 + k3 * ( -k7 ) ) / denom2 + yt;

        //X coordinates are too big, throw exception
        if ( fabs ( xc ) > MAX_POINT_COORDINATE )
        {
                throw BadDataException ( "BadDataException: too big coordinates of the circle mid point, xc > MAX_POINT_COORDINATE", " (MAX_POINT_COORDINATE = 1e+10)." ) ;
        }

        //Y coordinates are too big, throw exception
        if ( fabs ( yc ) > MAX_POINT_COORDINATE )
        {
                throw BadDataException ( "BadDataException: too big coordinates of the circle mid point, yc > MAX_POINT_COORDINATE", " (MAX_POINT_COORDINATE = 1e+10)." ) ;
        }

        //Radius of the circle
        r = sqrt ( ( x1r - ( xc - xt ) ) * ( x1r - ( xc - xt ) ) + ( x1r - ( yc - yt ) ) * ( x1r - ( yc - yt ) ) );
}

#endif
