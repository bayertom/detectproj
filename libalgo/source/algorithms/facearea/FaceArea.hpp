// Description: Compute area of the face

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


#ifndef  FaceArea_HPP
#define  FaceArea_HPP

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
T FaceArea::getFaceArea ( const Face <T> *face )
{
        //Get unsigned area of the Face stored in circular list
        return fabs ( getFaceAreaSigned ( face ) );
}


template <typename T>
T FaceArea::getFaceAreaSigned ( const Face <T> *face )
{
        //Get signed area of the Face stored in circular list
        if ( face != NULL )
        {
                T area = 0;
                const HalfEdge <T> *e_start = face->getHalfEdge();
                HalfEdge <T> *e = const_cast <HalfEdge <T> *> ( e_start );

                //There is a valid edge
                if ( e_start != NULL )
                {
                        //Get points
                        const Node3DCartesian <T> *pi = e->getPreviousEdge()->getPoint();
                        const Node3DCartesian <T> *pii = e->getPoint();

                        //Proces all edges of the Face
                        do
                        {
                                //Get point
                                const Node3DCartesian <T> *piii = e->getNextEdge()->getPoint();

                                //Compute area
                                area += pii->getX() * ( piii->getY() - pi->getY() );

                                //Assign points
                                pi = pii;
                                pii = piii;

                                //Increment edge
                                e = e->getNextEdge();

                        }
                        while ( e != e_start );
                }

                //Get area
                return 0.5 * area;
        }

        //Throw exception
        else
        {
                throw BadDataException ( "BadDataException: can not compute face area", " line segment = NULL." );
        }
}


template <typename Point, TDestructable destructable>
typename Point::Type FaceArea::getFaceArea ( const Container <Point *, destructable> *points )
{
        //Get unsigned area of the Face stored in list of vertices
        return fabs ( getFaceAreaSigned ( points ) );
}


template <typename Point, TDestructable destructable>
typename Point::Type FaceArea::getFaceAreaSigned ( const Container <Point *, destructable> *points )
{
        //Get signed area of the Face stored in list of vertices
        typename Point::Type area = 0;
        const unsigned int n = points->size();

        //Process all point
        for ( unsigned int i = 1; i < n - 1; i++ )
        {
                area += ( *points ) [i]->getX() * ( ( *points ) [i + 1]->getY() - ( *points ) [i - 1]->getY() );
        }

        //Compute x[0] * ( y[1] - y[n - 1]) and x[n-1] * ( y[0] - y[n - 2])
        area += ( *points ) [0]->getX() * ( ( *points ) [1]->getY() - ( *points ) [n - 1]->getY() );
        area += ( *points ) [n - 1]->getX() * ( ( *points ) [0]->getY() - ( *points ) [n - 2]->getY() );

        return 0.5 * area;
}

#endif
