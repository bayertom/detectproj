// Description: Compute perimeter of the face

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


#ifndef FacePerimeter_HPP
#define FacePerimeter_HPP

#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
T FacePerimeter::getFacePerimeter ( const Face <T> *f )
{
        //Get Face perimeter
        if ( f != NULL )
        {
                return getFacePerimeter ( f->getHalfEdge() );
        }

        //Throw exception
        else
        {
                throw BadDataException ( "BadDataException: can not compute face perimeter", " face = NULL." );
        }
}


template <typename T>
T FacePerimeter:: getFacePerimeter ( const HalfEdge <T> *e_start )
{
        //Get perimeter of the Face given by start half edge
        T perimeter = 0;
        HalfEdge <T> *e = const_cast <HalfEdge <T> *> ( e_start );

        //There is a valid edge
        if ( e_start != NULL )
        {
                //Get point
                Node3DCartesian <T> *pi = e->getPoint();

                //Proces all edges of the Face
                do
                {
                        //Get point
                        Node3DCartesian <T> *pii = e->getNextEdge()->getPoint();

                        //Compute area
                        perimeter += EuclDistance::getEuclDistance2D ( pi->getX(), pi->getY(), pii->getX(), pii->getY() );

                        //Assign point
                        pi = pii;

                        //Assign edge
                        e = e->getNextEdge();

                }
                while ( e != e_start );
        }

        //Get area
        return perimeter;
}

#endif
