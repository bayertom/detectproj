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


#ifndef Face_H
#define Face_H

#include <stdio.h>
#include <iostream>
#include <ostream>

#include "libalgo/source/structures/list/Container.h"


//Forward declaration
template <typename T>
class Point3DCartesian;

template <typename T>
class Node3DCartesian;

template <typename T>
class HalfEdge;


//Structure storing 2D Face using circular list
template <typename T>
class Face
{
        protected:
                HalfEdge <T> *edge;					//Any edge  of the Face
                mutable unsigned int vertices_count;                    //Total vertices of the Face

        public:
                template <TDestructable destructable>
                Face ( const Container <Node3DCartesian <T> *, destructable > &nl, Container <HalfEdge <T> *> &hl ) ;
                Face ( HalfEdge <T> *edge_ ) {edge = edge_; vertices_count = countVertices(); }
                Face ( const  Face <T> *f ) : edge ( f->edge ), vertices_count ( f->vertices_count ) {}
                virtual ~Face();

        protected:
                unsigned int countVertices() const;

        public:
                //Other methods
                HalfEdge <T> *getHalfEdge() const {return edge;}
                HalfEdge <T> *getPreviousHalfEdge () {return edge->getPreviousEdge();}
                HalfEdge <T> *getNextHalfEdge () {return edge -> getNextEdge();}
                unsigned int getVerticesCount() const {return  vertices_count == 0 ? vertices_count = countVertices() : vertices_count;}

                void setVerticesCount ( const unsigned int vertices_count_ ) {vertices_count = vertices_count_;}
                void setHalfEdge ( HalfEdge <T> *edge_ ) {edge = edge_;}

                void toPointsList ( Container <Point3DCartesian <T> > &pl ) const;
                void toNodesList ( Container <Node3DCartesian <T> *, NonDestructable> &nl ) const;

        public:
                virtual void removeAdjacency();
                virtual void print ( std::ostream * output = &std::cout ) const;
};

#include "Face.hpp"


template <typename T>
void operator << ( std::ostream & output, const Face <T> &f )
{
        //Print face: overloaded operator, common for all derived class
        f.print ( &output ) ;
}


#endif
