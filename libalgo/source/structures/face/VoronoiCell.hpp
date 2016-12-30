// Description: Structure representing Voronoi cell

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

#ifndef VoronoiCell_HPP
#define VoronoiCell_HPP

#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"


template <typename T>
VoronoiCell <T> ::~VoronoiCell()
{
        this->generator->setFace ( NULL );
        this->generator = NULL;
}


template <typename T>
void VoronoiCell <T>::removeAdjacency()
{
        //Set twin edges of Faces adjacent to actual Face to NULL
        if ( this->edge != NULL )
        {
                //Get actual edge in cell
                HalfEdge <T> *e = this->edge;

                //Proces all edges of the Face
                do
                {
                        //Get twin edge
                        HalfEdge <T> *e_twin = e->getTwinEdge();

                        //Get dual edge
                        HalfEdge <T> *e_dual = e->getDualEdge();

                        //Adjacent cell exists
                        if ( e_twin != NULL )
                        {
                                //Remove adjacency: Set pointer to NULL
                                e_twin -> setTwinEdge ( NULL );
                        }

                        //Dual edge exists
                        if ( e_dual != NULL )
                        {
                                //Remove duality: Set pointer to NULL
                                e_dual->setDualEdge ( NULL );
                        }

                        //Increment edge (including zero length edges)
                        e = e->getNextEdge();

                }
                while ( e != this->edge );
        }
}


template <typename T>
void VoronoiCell<T>::print ( std::ostream * output ) const
{
        //Print Voronoi cell vertices
        *output << std::fixed << std::setprecision ( 3 );

        *output << std::endl << "Generator: " << generator -> getX() << "  " << generator->getY() << std::endl;
        *output << "Bounded: " << bounded << std::endl << std::endl;
        Face <T>::print ( output );
}


#endif
