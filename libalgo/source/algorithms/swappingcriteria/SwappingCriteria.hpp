// Description: VArious swapping criteria for data depending triangulations

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


#ifndef SwappingCriteria_HPP
#define SwappingCriteria_HPP

#include <cmath>

#include "libalgo/source/structures/point/Node3DCartesian.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"
#include "libalgo/source/exceptions/MathInvalidArgumentException.h"
#include "libalgo/source/exceptions/BadDataException.h"

using namespace std;

template <typename T>
T SwappingCriteria ::getAbn ( const Point3DCartesian <T> *n1, const Point3DCartesian <T> *n2, const Point3DCartesian <T> *n3, const Point3DCartesian <T> *n4, const Point3DCartesian <T> *n5, const Point3DCartesian <T> *n6 )
{
        //Angle between normals criterion

        // Vectors 1
        const T	u11 = n2->getX() - n1->getX();
        const T	u12 = n2->getY() - n1->getY();
        const T	u13 = n2->getZ() - n1->getZ();
        const T	v11 = n3->getX() - n1->getX();
        const T	v12 = n3->getY() - n1->getY();
        const T	v13 = n3->getZ() - n1->getZ();

        // Normal 1
        const T	a_cnd = u12 * v13 - v12 * u13;
        const T	b1 = u13 * v11 - v13 * u11;
        const T	c1 = u11 * v12 - v11 * u12;

        // Vectors  2
        const T	u21 = n5->getX() - n4->getX();
        const T	u22 = n5->getY() - n4->getY();
        const T	u23 = n5->getZ() - n4->getZ();
        const T	v21 = n6->getX() - n4->getX();
        const T	v22 = n6->getY() - n4->getY();
        const T	v23 = n6->getZ() - n4->getZ();

        // Normal 2
        const T	a_and = u22 * v23 - v22 * u23;
        const T	b2 = u23 * v21 - v23 * u21;
        const T	c2 = u21 * v22 - v21 * u22;

        // Angle
        const T dot = a_cnd * a_and + b1 * b2 + c1 * c2;
        const T norm1 = sqrt ( a_cnd * a_cnd + b1 * b1 + c1 * c1 );
        const T norm2 = sqrt ( a_and * a_and + b2 * b2 + c2 * c2 );

        //Correct normals
        if ( norm1 == 0 || norm2 == 0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: ", "Can not compute ABN swap criterion.", norm1 );
        }

        //Compute acos (ABN)
        T abn_acos = dot / ( norm1 * norm2 );

        //Throw exception: fabs (abn_acos) - 1 > MIN_FLOAT
        if ( ( abn_acos > 1 + ARGUMENT_ROUND_ERROR ) || ( abn_acos < - 1 - ARGUMENT_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException, ", "ABN criterion: acos(arg), fabs(arg) > 1.", abn_acos );
        }

        //Correct round errors
        else if ( abn_acos > 1 )
        {
                abn_acos = 1;
        }

        else if ( abn_acos < - 1 )
        {
                abn_acos = -1;
        }

        //ABN criterion
        return  fabs ( acos ( abn_acos ) * 180 / M_PI );
}


template <typename T>
T SwappingCriteria::getSco ( const Point3DCartesian <T> *n1, const Point3DCartesian <T> *n2, const Point3DCartesian <T> *n3, const Point3DCartesian <T> *n4, const Point3DCartesian <T> *n5, const Point3DCartesian <T> *n6 )
{
        //Smoothness of contours
        // Vectors 1
        const T	u11 = n2->getX() - n1->getX();
        const T	u12 = n2->getY() - n1->getY();
        const T	u13 = n2->getZ() - n1->getZ();
        const T	v11 = n3->getX() - n1->getX();
        const T	v12 = n3->getY() - n1->getY();
        const T	v13 = n3->getZ() - n1->getZ();

        // Normal 1
        const T	a_cnd = u12 * v13 - v12 * u13;
        const T	b1 = u13 * v11 - v13 * u11;

        // Vectors  2
        const T	u21 = n5->getX() - n4->getX();
        const T	u22 = n5->getY() - n4->getY();
        const T	u23 = n5->getZ() - n4->getZ();
        const T	v21 = n6->getX() - n4->getX();
        const T	v22 = n6->getY() - n4->getY();
        const T	v23 = n6->getZ() - n4->getZ();

        // Normal 2
        const T	a_and = u22 * v23 - v22 * u23;
        const T	b2 = u23 * v21 - v23 * u21;

        // Angle
        const T dot = a_cnd * a_and + b1 * b2 ;
        const T norm1 = sqrt ( a_cnd * a_cnd + b1 * b1 );
        const T norm2 = sqrt ( a_and * a_and + b2 * b2 );

        //Throw exception
        if ( norm1 == 0 || norm2 == 0 )
        {
                return 0;
        }

        T sco_acos = dot / ( norm1 * norm2 ) ;


        //Throw exception: fabs (abn_acos) - 1 > MIN_FLOAT
        if ( ( sco_acos > 1 + ARGUMENT_ROUND_ERROR ) || ( sco_acos < - 1 - ARGUMENT_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException, ", "ABN criterion: acos(arg), fabs(arg) > 1.", sco_acos );
        }

        //Correct round errors
        else if ( sco_acos > 1 )
        {
                sco_acos = 1;
        }

        else if ( sco_acos < - 1 )
        {
                sco_acos = -1;
        }


        //SCO criterion
        return  acos ( sco_acos ) * 180 / M_PI ;
}


template <typename T>
bool SwappingCriteria::getClineRenka ( const Node3DCartesian <T> *p, const Node3DCartesian <T> *p1, const Node3DCartesian <T> *p2, const Node3DCartesian <T> *p3 )
{
        //Cline-Renka legality test: true, if quadrilateral illegal (it needs to be swapped)
        const T	x4 = p->getX();		// Tested point
        const T	y4 = p->getY();		// Tested point

        //Points coiunter clock wise oriented
        const T	x1 = p1->getX();	// First vertex in triangle adjacent to legalized
        const T	y1 = p1->getY();	// First vertex in triangle adjacent to legalized
        const T	x2 = p2->getX();	// Second vertex in triangle adjacent to legalized
        const T	y2 = p2->getY();	// Second vertex in triangle adjacent to legalized
        const T	x3 = p3->getX();	// Third vertex in triangle adjacent to legalized
        const T	y3 = p3->getY();	// Third vertex in triangle adjacent to legalized

        //Cline-Renka test
        T cos_a = ( x1 - x3 ) * ( x2 - x3 ) + ( y1 - y3 ) * ( y2 - y3 ) ;
        T cos_b = ( x2 - x4 ) * ( x1 - x4 ) + ( y2 - y4 ) * ( y1 - y4 ) ;

        //Triangle is legal: swap = false
        if ( ( cos_a >= 0 ) && ( cos_b >= 0 ) )
        {
                return false;
        }

        //Triangle is illegal: swap = true
        if ( ( cos_a < 0 ) && ( cos_b < 0 ) )
        {
                return true;
        }

        T sin_ab = ( ( x1 - x3 ) * ( y2 - y3 ) - ( x2 - x3 ) * ( y1 - y3 ) ) * cos_b +
                   ( ( x2 - x4 ) * ( y1 - y4 ) - ( x1 - x4 ) * ( y2 - y4 ) ) * cos_a;

        //Triangle is illegal: swap = true
        if ( sin_ab < 0 )
        {
                return true;
        }

        //Triangle is legal: swap = false
        return false;
}


template <typename T>
T SwappingCriteria::getEmptyCircleTest ( const Node3DCartesian <T> *p, const Node3DCartesian <T> *p1, const Node3DCartesian <T> *p2, const Node3DCartesian <T> *p3 )
{
        //Modified empty circle test by Okabe, 2000
        //Test legality of Delaunay triangle and point, also used in incremental construction of Voronoi diagram
        const T	x4 = p->getX();		// Test point
        const T	y4 = p->getY();		// Test point

        //Points counter clock wise oriented
        const T	x1 = p1->getX();	// First vertex in triangle adjacent to legalized
        const T	y1 = p1->getY();	// First vertex in triangle adjacent to legalized
        const T	x2 = p2->getX();	// Second vertex in triangle adjacent to legalized
        const T	y2 = p2->getY();	// Second vertex in triangle adjacent to legalized
        const T	x3 = p3->getX();	// Third vertex in triangle adjacent to legalized
        const T	y3 = p3->getY();	// Third vertex in triangle adjacent to legalized

        //Are not 2 points identical?
        if ( *p1 == *p2 || *p1 == *p3 || *p2 == *p3 )
        {
                throw BadDataException ( "BadDataException, can not compute empty circle test", "2 points are identical" );
        }

        //First determinant
        const T Jijk2 = ( y1 - y3 ) * ( ( x2 - x3 ) * ( x2 - x3 ) + ( y2 - y3 ) * ( y2 - y3 ) ) - ( ( x1 - x3 ) * ( x1 - x3 ) + ( y1 - y3 ) * ( y1 - y3 ) ) * ( y2 - y3 );

        //Second determinant
        const T Jijk3 = ( x1 - x3 ) * ( ( x2 - x3 ) * ( x2 - x3 ) + ( y2 - y3 ) * ( y2 - y3 ) ) - ( ( x1 - x3 ) * ( x1 - x3 ) + ( y1 - y3 ) * ( y1 - y3 ) ) * ( x2 - x3 );

        //Third determinant
        const T Jijk4 = ( x1 - x3 ) * ( y2 - y3 ) - ( y1 - y3 ) * ( x2 - x3 );

        //Criterion
        return  Jijk2 * ( x4 - x3 ) - Jijk3 * ( y4 - y3 ) + Jijk4 * ( ( x4 - x3 ) * ( x4 - x3 ) + ( y4 - y3 ) * ( y4 - y3 ) );
}

#endif
