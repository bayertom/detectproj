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


#ifndef PolyLine_H
#define PolyLine_H

#include <stdio.h>

#include "libalgo/source/structures/list/Container.h"

//Forward declaration
template <typename T>
class HalfEdge;

template <typename T>
class PolyLine
{
        protected:
                HalfEdge <T> *edge;                     //First edge  of the polyline
                unsigned int vertices_count;             //Total vertices of the polyline

        public:
                PolyLine () {edge = NULL; vertices_count = 0; }
                PolyLine ( const Container <Node3DCartesian <T> *> *nl ) ;
                PolyLine ( HalfEdge <T> *edge_ ) : edge ( edge_ ) { }
                virtual ~PolyLine();

        public:
                //Other methods
                HalfEdge <T> *getHalfEdge() const {return edge;}
                HalfEdge <T> *getPreviousHalfEdge () {return edge -> getPreviousEdge();}
                HalfEdge <T> *getNextHalfEdge() {return edge -> getNextEdge();}
                void setHalfEdge ( HalfEdge <T> *edge_ ) {edge = edge_;}
                unsigned int getVerticesCount() const {return vertices_count;}

                void toPoints3DCartesianList ( Container <Point3DCartesian<T> >&pl ) const;
                void toNodes3DCartesianList ( Container <Node3DCartesian <T> *> *nl ) const;
                void printPolyLine ( std::ofstream & file ) const;
};

#include "PolyLine.hpp"

#endif
