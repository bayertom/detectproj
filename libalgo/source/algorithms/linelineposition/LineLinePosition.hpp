// Description: Position of two lines

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


#ifndef LineLinePosition_HPP
#define LineLinePosition_HPP

#include "libalgo/source/algorithms/pointbetweenpoints/PointBetweenPoints.h"


template <typename T>
unsigned short LineLinePosition::get2LineSegmentsPosition ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, const Point3DCartesian <T> * p4 , T & x_int, T & y_int )
{
        //  Tests, if a line segment(p1,p2) is intersected by the line segment (p3,p4) and compute itersection.
        return get2LineSegmentsPosition ( p1->getX(), p1->getY(), p2->getX(), p2->getY(), p3->getX(), p3->getY(), p4->getX(), p4->getY(), x_int, y_int );
}


template <typename T>
unsigned short LineLinePosition::get2LinesPosition ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, const Point3DCartesian <T> * p4 , T & x_int, T & y_int )
{
        //  Tests, if a line (p1,p2) is intersected by the line (p3,p4) and compute itersection.
        return get2LinesPosition ( p1->getX(), p1->getY(), p2->getX(), p2->getY(), p3->getX(), p3->getY(), p4->getX(), p4->getY(), x_int, y_int );
}


template <typename T>
unsigned short LineLinePosition::get2LineSegmentsPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int )
{
        /*  Tests, if a line segment(p1,p2) is intersected by the line segment (p3,p4) and compute itersection. Results

        		7: Parallel,
        		6: Coliner (2 identical points),	x-------x
        		5: Colinear (1 identical point),	x-------x-------x
        		4: Colinear (No identical point),	x-------x       x-------x
        		3: Intersect in two end points,
        		2: Intersect in one end point,
        		1: Intersect (but not in end points),
        		0: Do not intersect ( not collinear nor parallel),
        */
        unsigned short intersection_code = 2;

        //Initialize intersections
        x_int = 0, y_int = 0;

        //Compute denominator
        const T denom =	( y4 - y3 ) * ( x2 - x1 ) - ( x4 - x3 ) * ( y2 - y1 ) ;

        //Segments are parallel
        if ( fabs ( denom ) < POSITION_ROUND_ERROR )
        {
                return get2ParallelsPosition ( x1, y1, x2, y2, x3, y3, x4, y4, x_int, y_int );
        }

        //Compute numerators
        const T numer1 = ( ( x4 - x3 ) * ( y1 - y3 ) - ( y4 - y3 ) * ( x1 - x3 ) );
        const T numer2 = ( ( x2 - x1 ) * ( y1 - y3 ) - ( y2 - y1 ) * ( x1 - x3 ) );

        //Compute parameters s,t
        const T s = numer1 / denom;
        const T t = numer2 / denom;

        //Both segments intersect in 2 end points
        if ( ( fabs ( s ) < POSITION_ROUND_ERROR )  && ( fabs ( t ) < POSITION_ROUND_ERROR ) ||
                        ( fabs ( s ) < POSITION_ROUND_ERROR )  && ( fabs ( t - 1.0 ) < POSITION_ROUND_ERROR ) ||
                        ( fabs ( s - 1.0 ) < POSITION_ROUND_ERROR )  && ( fabs ( t ) < POSITION_ROUND_ERROR ) ||
                        ( fabs ( s - 1.0 ) < POSITION_ROUND_ERROR ) && ( fabs ( t - 1.0 ) < POSITION_ROUND_ERROR ) )
        {
                intersection_code =  3;
        }

        //Both segments intersect in one end point
        else if ( ( fabs ( s ) < POSITION_ROUND_ERROR )  && ( fabs ( t ) < POSITION_ROUND_ERROR ) && ( fabs ( t - 1.0 ) < POSITION_ROUND_ERROR ) ||
                        ( fabs ( t ) < POSITION_ROUND_ERROR )  && ( fabs ( s ) < POSITION_ROUND_ERROR ) && ( fabs ( s - 1.0 ) < POSITION_ROUND_ERROR ) )
        {
                intersection_code =  2;
        }

        //Segments do not intersect: do not compute any intersection
        else if ( ( s < 0.0 ) || ( s > 1 ) ||
                        ( t < 0.0 ) || ( t > 1 ) )
        {
                return  0;
        }

        //Segments intersect, but not in end points
        else if ( ( s > 0.0 ) && ( s < 1.0 ) && ( t > 0.0 ) && ( t < 1.0 ) )
        {
                intersection_code =  1;
        }

        //Compute intersection
        x_int = x1 + s * ( x2 - x1 );
        y_int = y1 + s * ( y2 - y1 );

        //Segments intersect in one end point
        return intersection_code;
}


template <typename T>
unsigned short LineLinePosition::get2ParallelsPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int )
{
        //Get position of 2 parallells

        //Lines are parallel, otherwise are colinear
        if ( PointLineDistance::getPointLineDistance2D ( x3, y3, x1, y1, x2, y2 ) > POSITION_ROUND_ERROR )
        {
                return 7;
        }

        //Test 1, if end point are identical (p1 = p3, p2 = p4)
        const bool t1 = ( ( x1 - x3 ) * ( x1 - x3 ) + ( y1 - y3 ) * ( y1 - y3 ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR );
        const bool t2 = ( ( x2 - x4 ) * ( x2 - x4 ) + ( y2 - y4 ) * ( y2 - y4 ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR );

        //Both points are identical: there is no valid intersection
        if ( t1 && t2 ) return 6;

        //1 point is identical: there is an intersection
        if ( ( t1 && ( true || ( x_int = x1 ) ) && ( true || ( y_int = y1 ) ) ) ||
                        ( t2 && ( true || ( x_int = x2 ) ) && ( true || ( y_int = y2 ) ) ) )
        {
                return 5;
        }

        //Test 2, if end point are identical (p1 = p4, p2 = p3)
        const bool t3 = ( ( x1 - x4 ) * ( x1 - x4 ) + ( y1 - y4 ) * ( y1 - y4 ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR );
        const bool t4 = ( ( x2 - x3 ) * ( x2 - x3 ) + ( y2 - y3 ) * ( y2 - y3 ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR );

        //Both points are identical: there is no valid intersection
        if ( t3 && t4 ) return 6;

        //1 point is identical: there is an intersection
        if ( ( t3 && ( true || ( x_int = x1 ) ) && ( true || ( y_int = y1 ) ) ) ||
                        ( t4 && ( true || ( x_int = x2 ) ) && ( true || ( y_int = y2 ) ) ) )
        {
                return 5;
        }

        //No point is identical
        return 4;
}



template <typename T>
bool LineLinePosition::get2LinesPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int )
{
        /* Tests, if a line (p1,p2) is intersected by the line  (p3,p4) and compute itersection. Return:
        		1: Intersection exist
        		0: Intersecion does not exist
        */

        //Initialize intersections
        x_int = 0, y_int = 0;

        //Vectors
        const T ux = x2 - x1, uy = y2 - y1;
        const T vx = x4 - x3, vy = y4 - y3;
        const T wx = x1 - x3, wy = y1 - y3;

        //Compute k1, k2
        const T k1 = vx * wy - vy * wx;
        const T k2 = vy * ux - vx * uy;

        //Intersection exists, compute coordinates of the intersection point
        if ( fabs ( k2 ) > POSITION_ROUND_ERROR )
        {
                //Compute s
                const T s = k1 / k2;

                //Intersections
                x_int = x1 + s * ux;
                y_int = y1 + s * uy;

                //Intersection code could be set using i
                return true;
        }

        //Intersection does not exist
        return false;
}



#endif
