// Description: Position of the point and face using ray algorithm (O�Rourke implementation)

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


#ifndef PointFacePosition_HPP
#define PointFacePosition_HPP

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/pointlinedistance/PointLineDistance.h"


template <typename T>
unsigned short PointFacePosition:: getPointFacePosition ( const Point3DCartesian <T> * q, const Face <T> *face )
{
        /*Returns position of the point to the Face
        0: q is strictly interior
        1: q is strictly exterior
        2: q is on the edge but not an endpoint
        3: q is a vertex
        */

        //There is a valid Face
        if ( face != NULL )
        {
                unsigned short l_inters = 0, r_inters = 0;

                HalfEdge <T> *e = face->getHalfEdge();
                const HalfEdge <T> *e_start = e;

                //Get point P[i-1]
                const Node3DCartesian <T> *pi1 = e->getPreviousEdge()->getPoint();

                //Process all edges
                do
                {
                        //Get point P[i]
                        const Node3DCartesian <T> *pi = e->getPoint();

                        //Test point is identical with end points: d(q, pi) < POSITION_ROUND_ERROR
                        if ( ( q->getX() - pi->getX() ) * ( q->getX() - pi->getX() ) + ( q->getY() - pi->getY() ) * ( q->getY() - pi->getY() ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR )
                        {
                                return 3;
                        }

                        //Perform  tests: different position of segment end points to ray
                        bool t1 = ( pi->getY() > q->getY() + POSITION_ROUND_ERROR ) != ( pi1->getY() > q->getY() + POSITION_ROUND_ERROR );
                        bool t2 = ( pi->getY() < q->getY() - POSITION_ROUND_ERROR ) != ( pi1->getY() < q->getY() - POSITION_ROUND_ERROR );

                        //Segment intersects left or ray
                        if ( t1 || t2 )
                        {
                                //Compute intersection
                                T x = ( ( pi->getX() - q->getX() ) * ( pi1->getY() - q->getY() ) - ( pi1->getX() - q->getX() ) * ( pi->getY() - q->getY() ) )
                                      / double ( ( pi1->getY() - pi->getY() ) );

                                //Aditional test to avoid roundness error: point too close to edge
                                if ( PointLineDistance::getPointLineSegmentDistance2D ( q, pi, pi1 ) < POSITION_ROUND_ERROR )
                                {
                                        return 2;
                                }

                                //Increment number of right intersections
                                if ( t1 && x > 0 )
                                {
                                        r_inters ++;
                                }

                                //Increment  number of left intersections
                                if ( t2 && x < 0 )
                                {
                                        l_inters ++;
                                }
                        }

                        //Assign point
                        pi1 = pi;

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != e_start );

                //Point lies on the edge, but not at endpoints
                if ( ( l_inters % 2 ) != ( r_inters % 2 ) )
                {
                        return 2;
                }

                //Points is strictly interior
                if ( l_inters % 2 == 1 )
                {
                        return 0;
                }

                //Point is strictly exterior
                return 1;
        }

        //Throw exception
        else
        {
                throw BadDataException ( "BadDataException: no face incident to edge, ", "can not test point vs. face position" );
        }
}



template <typename Point, TDestructable destructable>
unsigned short PointFacePosition:: getPointFacePosition ( const Point * q, const Container <Point*, destructable> &face )
{
        /*Returns position of the point to the Face
        0: q is strictly interior
        1: q is strictly exterior
        2: q is on the edge but not an endpoint
        3: q is a vertex
        */
        unsigned short l_inters = 0, r_inters = 0;
        const unsigned int n = face.size();

        //Process all edges
        for ( unsigned int i = 0; i < n; i++ )
        {
                //Get point P[i]
                const Point * pi = face [i];

                //Test point is identical with end points: d(q, pi) < POSITION_ROUND_ERROR
                if ( ( q->getX() - pi->getX() ) * ( q->getX() - pi->getX() ) + ( q->getY() - pi->getY() ) * ( q->getY() - pi->getY() ) < POSITION_ROUND_ERROR * POSITION_ROUND_ERROR )
                {
                        return 3;
                }

                //Get previous index
                int i1 = ( i + n - 1 ) % n;

                //Get point P[i]
                const Point * pi1 = face [i1];

                //Perform  tests: different position of segment end points to ray
                bool t1 = ( pi->getY() > q->getY() + POSITION_ROUND_ERROR ) != ( pi1->getY() > q->getY() + POSITION_ROUND_ERROR );
                bool t2 = ( pi->getY() < q->getY() - POSITION_ROUND_ERROR ) != ( pi1->getY() < q->getY() - POSITION_ROUND_ERROR );

                //Segment intersects left or ray
                if ( t1 || t2 )
                {
                        //Compute intersection
                        typename Point::Type x = ( ( pi->getX() - q->getX() ) * ( pi1->getY() - q->getY() ) - ( pi1->getX() - q->getX() ) * ( pi->getY() - q->getY() ) )
                                                 / double ( ( pi1->getY() - pi->getY() ) );

                        //Aditional test to avoid roundness error: point too close to edge
                        if ( PointLineDistance::getPointLineSegmentDistance2D ( q, pi, pi1 ) < POSITION_ROUND_ERROR )
                        {
                                return 2;
                        }

                        //Increment number of right intersections
                        if ( ( t1 ) && ( x > 0 ) )
                        {
                                r_inters ++;
                        }

                        //Increment  number of left intersections
                        if ( ( t2 ) && ( x < 0 ) )
                        {
                                l_inters ++;
                        }
                }
        }

        //Point lies on the edge, but not at endpoints
        if ( ( l_inters % 2 ) != ( r_inters % 2 ) )
        {
                return 2;
        }

        //Points is strictly interior
        if ( l_inters % 2 == 1 )
        {
                return 0;
        }

        //Point is strictly exterior
        return 1;
}


#endif
