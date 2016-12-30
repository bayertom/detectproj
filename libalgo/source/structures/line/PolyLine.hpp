// Description: Polyline definition

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


#ifndef PolyLine_HPP
#define PolyLine_HPP

#include "libalgo/source/structures/line/HalfEdge.h"


template <typename T>
PolyLine <T> :: PolyLine ( const Container <Node3DCartesian <T> *> *nl )
{
        //Create polyline from the list of nodes
        this->vertices_count = nl->size();

        //Create half edges and set previous edges
        HalfEdge <T> *h = NULL, * h_old = NULL, *h_start = NULL;

        for ( unsigned int i = 0; i < this->vertices_count; i++ )
        {
                if ( h != NULL )
                {
                        h = new HalfEdge <T> ( ( *nl ) [i], h_old, NULL, NULL );
                }

                //Create first edge
                else
                {
                        h = new HalfEdge <T> ( ( *nl ) [i], NULL, NULL, NULL );
                        h_start = h;
                }

                //Assign edge
                h_old = h;
        }

        //Set next edge for other half edges
        while ( h_old != h_start )
        {
                //Get previous half edge
                HalfEdge <T> *h_prev = h_old->getPreviousEdge();

                //Set next edge for previous edge
                h_prev->setNextEdge ( h_old );

                //Decrement edge
                h_old = h_old->getPreviousEdge();
        }

        //Set start edge of the Face
        this->edge = h_start;

        //Set Face for all half edges
        h = h_start;

        for ( unsigned int i = 0; i < vertices_count; i++, h = h->getNextEdge() )
        {
                //Set Face
                h->setFace ( NULL );
        }
}


template <typename T>
PolyLine <T> ::~PolyLine()
{
        //Get half edge
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;
                HalfEdge <T> *e_prev = NULL;

                //Proces all edges of the Face
                for ( unsigned int i = 0; i < vertices_count ; i++ )
                {
                        //Get edge of the polyline
                        e_prev = e;

                        //First increment edge
                        e = e->getNextEdge();

                        //Then delete edge
                        delete e_prev;

                        //Set edge to NULL
                        e_prev = NULL;
                }
        }
}


template <typename T>
void PolyLine <T> ::printPolyLine ( std::ofstream & file ) const
{
        //Print polyline
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                for ( unsigned int i = 0; i < vertices_count; i++, e = e->getNextEdge() )
                {
                        //Print half edge
                        e->printHalfEdge ( file );
                }
        }
}


template <typename T>
void PolyLine <T> ::toPoints3DCartesianList ( Container <Point3DCartesian<T> >&pl ) const
{
        //Convert polyline to list
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                //Proces all edges of the Face
                for ( unsigned int i = 0; i < vertices_count; i++, e = e->getNextEdge() )
                {
                        //Add point to the list
                        pl.push_back ( * ( e->getPoint() ) );
                }
        }
}


template <typename T>
void PolyLine <T> ::toNodes3DCartesianList ( Container <Node3DCartesian <T> *> *nl ) const
{
        //Convert polyline to nodes list
        if ( edge != NULL )
        {
                //Get next edge in cell
                HalfEdge <T> *e = edge;

                for ( unsigned int i = 0; i < vertices_count; i++, e = e->getNextEdge() )
                {
                        //Add node to the list
                        nl->push_back ( new Node3DCartesian <T> ( e->getPoint() ) );
                }
        }
}

#endif
