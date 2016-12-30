// Description: Lawson oriented walk, find face containing point p using vector orientation test

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


#ifndef LawsonOrientedWalk_HPP
#define LawsonOrientedWalk_HPP

#include <stdio.h>

#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"


template <typename T>
HalfEdge <T> *LawsonOrientedWalk::findTriangleWalk ( const Node3DCartesian <T> *p, short * status, HalfEdge <T> *e )
{
        /*
         * Find triangle, inside that lies the point p = [xp,yp] using Lawson oriented search ;
         * status=1: point lies inside the triangle, status=2: point lies on the edge of
         * the triangle ;
         * Used heuristic: finding starts from the last triangle defined by edge e
         */

        unsigned short  t1 = 0, t2 = 0, t3 = 0;

        // Check first edge
        t1 = PointLinePosition::getPointLinePosition2D ( p, e->getPoint(), ( e->getNextEdge() )->getPoint() );

        if ( t1 > 0 )
        {
                // Check next edge, move pointer
                e = e->getNextEdge();
        }

        // Change orientation of the half plane
        else
        {
                // Set new orientation
                t1 = 1;

                // Change orientation of the edge, increment edge
                e = e->getTwinEdge();

                // Triangle lies on the boundary
                if ( e == NULL )
                {
                        return NULL;
                }

                // Set next edge
                e = e->getNextEdge();
        }

        //Run until any suitable triangle is found
        for ( ;; )
        {
                // Lies the point in the intersection of two half planes?
                t2 = PointLinePosition::getPointLinePosition2D ( p, e->getPoint(), e->getNextEdge()->getPoint() );

                if ( t2 > 0 )
                {
                        // Check the third edge, move pointer
                        e = e->getNextEdge();

                        // Points lies inside a triangle (intersection of three half planes)
                        t3 = PointLinePosition::getPointLinePosition2D ( p, e->getPoint(), ( e->getNextEdge() )->getPoint() );

                        if ( t3 > 0 )
                        {
                                // Point inside triangle
                                if ( t1 * t2 * t3 == 1 )
                                {
                                        *status = 1;
                                }

                                // Point lies on the edge of the triangle
                                else
                                {
                                        *status = 2;

                                        // Return edge point lies on
                                        if ( t3 == 2 )
                                        {
                                                // Actual edge
                                                return e;
                                        }

                                        if ( t1 == 2 )
                                        {
                                                // First edge
                                                return e->getNextEdge();
                                        }

                                        if ( t2 == 2 )
                                        {
                                                // Second edge
                                                return e->getNextEdge()->getNextEdge();
                                        }
                                }

                                // Triangle was found, return half edge
                                return e;
                        }

                        // Change orientation of the half plane
                        else
                        {
                                // This edge becomes first edge, move pointer
                                e = e->getTwinEdge();

                                // No triangle found
                                if ( e == NULL )
                                {
                                        return NULL;
                                }

                                // But we test the next edge, move pointer
                                e = e->getNextEdge();

                                // Set orientation
                                t1 = 1;
                        }
                }

                // Change orientation of the half plane
                else
                {
                        // This edge becomes first edge
                        e = e->getTwinEdge();

                        // Triangle lies on the boundary
                        if ( e == NULL )
                        {
                                return NULL;
                        }

                        // But we test the next edge
                        e = e->getNextEdge();

                        // Set orientation for this edge
                        t1 = 1;
                }
        }

        //No triangle was found;
        return NULL;
}

template <typename T>
Face <T> *LawsonOrientedWalk::findFaceWalk ( const Node3DCartesian <T> *p, short * status, Face <T> *pol )
{
        /*
         * Find CONVEX Face, inside that lies the point p = [xp,yp] using Lawson oriented search ;
         * status = 0: point lies outside of Face,
         * status = 1: point lies inside the Face,
         * status = 2: point lies on the edge of the Face ;
         * Used heuristic: searching begins from the last Face defined by edge e
         */
        unsigned short  t = 0;

        //Start edge of the Face
        HalfEdge <T> *e = pol->getHalfEdge();
        const HalfEdge <T> *e_start = e;

        //Process all halfedges of all Faces
        for ( ;; )
        {
                //Compute half plane test
                t = PointLinePosition::getPointLinePosition2D ( p, e->getPoint(), ( e->getNextEdge() )->getPoint() ); //Was off

                //Point colinear with some edge
                if ( t == 2 )
                {
                        *status = 2;

                        //Increment edge
                        e = e->getNextEdge();
                }

                //Point lefts to the edge
                else if ( t == 1 )
                {
                        //Increment edge
                        e = e->getNextEdge();
                }

                //Point rights to the edge
                else if ( t == 0 )
                {
                        //Change Face
                        e = e->getTwinEdge();

                        //We reach boundary edge, return NULL
                        if ( e == NULL )
                        {
                                *status = 0;

                                //No Face has been found
                                return NULL;
                        }

                        //Remeber this edge, new start edge of the Face
                        e_start = e;

                        //Increment edge
                        e = e->getNextEdge();
                }

                //We reach start edge of the Face, return e;
                if ( e == e_start )
                {
                        //Point does not lie on the edge
                        if ( *status != 2 )
                        {
                                *status = 1;
                        }

                        //Get Face
                        pol = e->getFace();

                        //Return Face
                        return e->getFace() ;
                }
        }

        //Return null
        return NULL;
}


template <typename T>
HalfEdge <T> *LawsonOrientedWalk::findFaceWalk2 ( const Node3DCartesian <T> *p, short * status, HalfEdge <T> *e, const T max_steps )
{
        /*
         * Find CONVEX Face, inside that lies the point p = [xp,yp] using Lawson oriented search ;
         * status = 0: point lies outside of Face,
         * status = 1: point lies inside the Face,
         * status = 2: point lies on the edge of the Face ;
         * Used heuristic: searching begins from the last Face defined by edge e
         */
        *status = 0;

        //Start edge of the Face
        HalfEdge <T> *e_start = e;

        //Process all halfedges of all Faces
        for ( unsigned int i = 0; i < max_steps; i++)
        {
                //Compute half plane test
                *status = PointLinePosition::getPointLinePosition2D ( p, e->getPoint(), ( e->getNextEdge() )->getPoint() );

                //Point colinear with some edge
                if ( *status == 2 )
                {
                        //Point lies on the edge, perform aditional test:
                        //a) to avoid the test of next edges of the triangle
                        //b) to avoid the numerical inacuracy

                        //Distance (p, e->start)
                        T dist_p_start = EuclDistance::getEuclDistance2D ( p->getX(), p->getY(), e->getPoint()->getX(), e->getPoint()->getY() );

                        //Distance (p, e->end)
                        T dist_p_end = EuclDistance::getEuclDistance2D ( p->getX(), p->getY(), e->getNextEdge()->getPoint()->getX(), e->getNextEdge()->getPoint()->getY() );

                        //Distance (e->start, e->end)
                        T dist_start_end = EuclDistance::getEuclDistance2D ( e->getPoint()->getX(), e->getPoint()->getY(),
                                           e->getNextEdge()->getPoint()->getX(), e->getNextEdge()->getPoint()->getY() );

                        //Does a point lie on the line (e->start, e->end)
                        if ( ( dist_p_start <= dist_start_end ) && ( dist_p_end  <= dist_start_end ) )
                        {
                                return e;
                        }

                        //Increment edge
                        e = e->getNextEdge();

                        //Change Face
                        e = e->getTwinEdge();

                        //We reach boundary edge, return NULL
                        if ( e == NULL )
                        {
                                *status = 0;

                                //No Face has been found
                                return NULL;
                        }

                        //Set start edge
                        e_start = e;

                        continue;
                }

                //Point lefts to the edge
                else if ( *status == 1 )
                {
                        //Increment edge
                        e = e->getNextEdge();
                }

                //Point rights to the edge
                else if ( *status == 0 )
                {
                        //Change Face
                        e = e->getTwinEdge();

                        //We reach boundary edge, return NULL
                        if ( e == NULL )
                        {
                                *status = 0;

                                //No Face has been found
                                return NULL;
                        }

                        //Remeber this edge, new start edge of the Face
                        e_start = e;

                        //Increment edge
                        e = e->getNextEdge();
                }

                //We reach start edge of the Face, return e;
                if ( e == e_start )
                {
                        //Return half edge
                        return e ;
                }
        }

        //Return null
        return NULL;
}

#endif
