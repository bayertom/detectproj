// Description: 2D half edge class

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


#ifndef HalfEdge_H
#define HalfEdge_H

#include <iostream>
#include <ostream>

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class Face;

//Half edge data structure
template <typename T>
class HalfEdge
{
        private:
                Node3DCartesian <T> *start; 		//Start point of the half edge
                HalfEdge <T> *previous;			//Pointer to previoust half edge in triangle (not used in DT)
                HalfEdge <T> *next;			//Pointer to next half edge in triangle
                HalfEdge <T> *twin;			//Pointer to half edge with the different orientation
                HalfEdge <T> *dual;			//Pointer to the dual structure (Delaunay / Voronoi) edge
                Face <T> *face;				//Pointer to Face contains this half edge
                bool simplex;				//Is this half edge member of the simplex triangle? (sometimes as different indicator)


        public:
                HalfEdge() : start ( NULL ), previous ( NULL ), next ( NULL ), twin ( NULL ), dual ( NULL ), face ( NULL ), simplex ( false ) {};
                HalfEdge ( Node3DCartesian <T> *start_, HalfEdge <T> *next_, HalfEdge <T> *twin_ ) :
                        start ( start_ ), previous ( NULL ), next ( next_ ), twin ( twin_ ), dual ( NULL ), face ( NULL ), simplex ( false ) {}
                HalfEdge ( Node3DCartesian <T> *start_, HalfEdge <T> *previous_, HalfEdge <T> *next_, HalfEdge <T> *twin_ ) :
                        start ( start_ ), previous ( previous_ ), next ( next_ ), twin ( twin_ ), dual ( NULL ), face ( NULL ), simplex ( false )  {}
                HalfEdge ( const HalfEdge <T> *e ) : start ( e->start ), previous ( e->previous ), next ( e->next ), twin ( e->twin ),
                        dual ( e->dual ), face ( e->face ), simplex ( e->simplex ) {}

                ~HalfEdge();

        public:
                HalfEdge <T> *operator ++ () {return next;}
                HalfEdge <T> *operator -- () {return previous;}

        public:
                Node3DCartesian <T> *getPoint () const {return start;}

                HalfEdge <T> *getPreviousEdge() const {return previous; }
                HalfEdge <T> *getNextEdge() const {return next ;}

                HalfEdge <T> *getTwinEdge() const {return twin;}
                HalfEdge <T> *getDualEdge() const {return dual;}
                Face <T> *getFace() const {return face;}

                bool isBoundaryEdge() const {bool result; twin == NULL ? result = true : result = false; return result;}
                bool isSimplexEdge() const {return simplex;}

                void setPoint ( Node3DCartesian <T> *p ) {start = p;}
                void setPreviousEdge ( HalfEdge <T> *e ) {previous = e;}
                void setNextEdge ( HalfEdge <T> *e ) {next = e;}
                void setTwinEdge ( HalfEdge <T> *e ) {twin = e;}
                void setDualEdge ( HalfEdge <T> *e ) {dual = e;}
                void setFace ( Face <T> *face_ ) {face = face_;}
                void setEdgeAsSimplex ( const bool simplex_ ) {simplex = simplex_;}

        public:
                virtual void print ( std::ostream * output = &std::cout ) const;
                friend void operator << ( std::ostream & output, const HalfEdge <T> &e ) { e.print ( &output ) ; }
};

#include "HalfEdge.hpp"

#endif
