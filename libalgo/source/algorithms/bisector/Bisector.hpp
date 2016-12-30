// Description: Compute normalized bisector and left-oriented normalized bisector

// Copyright (c) 2010 - 2013
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//libalgo/source/
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


#ifndef Bisector_HPP
#define Bisector_HPP

#include <cmath>

#include "libalgo/source/algorithms/vectorvectororientation/VectorVectorOrientation.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
void Bisector::getNormalizedLeftBisector2D ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, T & p_bis_x, T & p_bis_y,
                T & n_vect_bis_x, T & n_vect_bis_y )
{
        //Compute normalized bisector vector (nbx, nby) oriented left to the edge (p1, p2)
        const T norm = sqrt ( ( p2->getX() - p1->getX() ) * ( p2->getX() - p1->getX() ) +
                              ( p2->getY() - p1->getY() ) * ( p2->getY() - p1->getY() ) );

        //Throw exception
        if ( norm < MIN_FLOAT )
        {
                throw  MathZeroDevisionException <T> ( "MathZeroDevisionException, ", "can not compute normalized bisector: 1/norm, norm = 0." );
        }

        //Normalized vector representing left bisector from mid point of the line
        p_bis_x = 0.5 * ( p1->getX() + p2->getX() ) ;
        p_bis_y = 0.5 * ( p1->getY() + p2->getY() ) ;

        n_vect_bis_x = ( p1->getY() - p2->getY() ) / norm;
        n_vect_bis_y = ( p2->getX() - p1->getX() ) / norm;
}


template <typename T>
void Bisector::getNormalizedBisector2D ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3,
                T & n_vect_bis_x, T & n_vect_bis_y )
{
        //Compute normalized left bisector of two lines given by points (p1, p2) and (p3, p4)
        const T nu = sqrt ( ( p1->getX() - p2->getX() ) * ( p1->getX() - p2->getX() ) +
                            ( p1->getY() - p2->getY() ) * ( p1->getY() - p2->getY() ) );
        const T nv = sqrt ( ( p3->getX() - p2->getX() ) * ( p3->getX() - p2->getX() ) +
                            ( p3->getY() - p2->getY() ) * ( p3->getY() - p2->getY() ) );

        //Throw exception
        if ( nu < MIN_FLOAT )
        {
                throw  MathZeroDevisionException <T> ( "MathZeroDevisionException", "can not compute normalized bisector: 1/norm (u), norm(u) = 0.", nu );
        }

        if ( nv < MIN_FLOAT )
        {
                throw  MathZeroDevisionException <T> ( "MathZeroDevisionException, ", "can not compute normalized bisector: 1/norm (v), norm(v) = 0.", nv );
        }

        //Compute bisector
        T vect_bis_x = ( p1->getX() - p2->getX() ) / nu + ( p3->getX() - p2->getX() ) / nv;
        T vect_bis_y = ( p1->getY() - p2->getY() ) / nu + ( p3->getY() - p2->getY() ) / nv;

        //Vectors are collinear and negatively oriented
        if ( ( fabs ( vect_bis_x ) < POSITION_ROUND_ERROR ) && ( fabs ( vect_bis_x ) < POSITION_ROUND_ERROR ) )
        {
                vect_bis_x =  - ( p1->getY() - p2->getY() ) / nu;
                vect_bis_y = ( p1->getX() - p2->getX() ) / nu;
        }

        //Switch orientation, if neccessary: result = left bisector
        if ( VectorVectorOrientation::getVectorVectorOrientation2D ( p1->getX() - p2->getX(), p1->getY() - p2->getY(), vect_bis_x, vect_bis_y ) == 1 )
        {
                vect_bis_x *= -1;
                vect_bis_y *= -1;
        }

        //Compute norm
        const T vect_bis_n = sqrt ( vect_bis_x * vect_bis_x + vect_bis_y * vect_bis_y );

        //Normalize bisector
        n_vect_bis_x = vect_bis_x / vect_bis_n;
        n_vect_bis_y = vect_bis_y / vect_bis_n;
}

#endif
