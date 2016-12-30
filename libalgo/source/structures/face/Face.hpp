// Description: Class storing face

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


#ifndef Face_HPP
#define Face_HPP

#include <iostream>

#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"

template <typename T>
template <TDestructable destructable>
Face <T> :: Face ( const Container <Node3DCartesian <T> *, destructable> &nl, Container <HalfEdge <T> *> &hl )
{
        //Create Face from the list of nodes
        this->vertices_count = nl.size();
        this->edge = NULL;

        //Nodes list contains some points
        if ( this->vertices_count > 2 )
        {
                //Create first half edge
                HalfEdge <T> *h = new HalfEdge <T> ( nl [0], NULL, NULL, NULL ); //Create first edge
                HalfEdge <T> *h_start = h;
                hl.push_back ( h );

                //Create other half edges
                HalfEdge <T> *h_old = h;

                //Process all vertices
                for ( unsigned int i = 1; i < this->vertices_count; i++, h_old = h )
                {
                        //Create new edge and set previous edge
                        h = new HalfEdge <T> ( nl [i], h_old, NULL, NULL );

                        //Set next edge
                        h_old->setNextEdge ( h );

                        //Set face for each edge
                        h->setFace ( this );

                        //Add half edge to the list
                        hl.push_back ( h );
                }

                //Set previous edge of the first edge
                h_start->setPreviousEdge ( h );

                //Set next edge for the last edge
                h->setNextEdge ( h_start );

                //Set face for h start
                h_start->setFace ( this );

                //Set start edge of the Face
                this->edge = h_start;
        }
}


template <typename T>
Face <T> ::~Face()
{
        edge = NULL;
}


template <typename T>
unsigned int Face <T> ::countVertices() const
{
        //Count vertices of the Face
        unsigned int vertices = 0;

        //There is a valid edge
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                do
                {
                        //Increment vertices
                        vertices ++;

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != edge );
        }

        return vertices;
}


template <typename T>
void Face <T>::removeAdjacency()
{
        //Set twin edges of Faces adjacent to actual Face to NULL
        if ( edge != NULL )
        {
                //Get actual edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                do
                {
                        //Get twin edge
                        HalfEdge <T> *e_twin = e->getTwinEdge();

                        //Adjacent cell exists
                        if ( e_twin != NULL )
                        {
                                //Set pointer to NULL
                                e_twin -> setTwinEdge ( NULL );
                        }

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != edge );
        }
}



template <typename T>
void Face <T> ::print ( std::ostream * output ) const
{
        //Print Face
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                do
                {
                        //Print edge
                        e->print ( output );

                        //Add new line
                        *output << '\n';

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != edge );
        }

        std::cout << std::endl;
}


template <typename T>
void Face <T> ::toPointsList ( Container <Point3DCartesian <T> >&pl ) const
{
        //Convert Face to nodes list
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                do
                {
                        //Add point to the list
                        pl.push_back ( * ( e->getPoint() ) );

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != edge );
        }
}


template <typename T>
void Face <T> ::toNodesList ( Container <Node3DCartesian <T> *, NonDestructable > &nl ) const
{
        //Convert Face to nodes list
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                do
                {
                        //Add point to the list
                        nl.push_back ( e->getPoint() );

                        //Increment edge
                        e = e->getNextEdge();

                }
                while ( e != edge );
        }
}

#endif
