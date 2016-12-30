// Description: 3D cartesian node storing pointers to face / half edge used for topological operations
// Used in DCEL model

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

#ifndef Node3DCartesian_H
#define Node3DCartesian_H

#include <stdio.h>
#include <iostream>
#include <ostream>

#include "Point3DCartesian.h"

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class Face;


//2D point storing topological links (incident half edge and face)
//Used in DCEL model
template <typename T>
class Node3DCartesian: virtual public Point3DCartesian <T>
{
        protected:

                //Pointer to incident half edge
                HalfEdge <T> *half_edge;

                //Pointer to incident Face
                Face <T> *face;

        public:
                Node3DCartesian() : Point3DCartesian <T> () , half_edge ( NULL ), face ( NULL ) { }
                Node3DCartesian ( const T x, const T y, const T z = 0 ) : Point3DCartesian <T> ( x, y, z ), half_edge ( NULL ), face ( NULL ) { }
                Node3DCartesian ( const std::string & point_label, const T x, const T y, const T z = 0 ) : Point3DCartesian <T> ( point_label, x, y, z ), half_edge ( NULL ), face ( NULL ) { }
                Node3DCartesian ( const std::string & point_label, const T x, const T y, const T z, HalfEdge <T> * half_edge_, Face <T> *face_ ) : Point3DCartesian <T> ( point_label, x, y, z ),
                        half_edge ( half_edge_ ), face ( face_ ) { }
                Node3DCartesian ( const Node3DCartesian <T> *p ) : Point3DCartesian <T> ( p ), half_edge ( p->half_edge ), face ( p->face ) {};
                virtual ~Node3DCartesian() { half_edge = NULL; face = NULL;}

        public:
                //Operators
                bool operator == ( const Node3DCartesian <T> &p ) const;
                bool operator != ( const Node3DCartesian <T> &p ) const { return ! ( *this == p );}

        public:

                //Functions
                HalfEdge <T> *getHalfEdge() const {return half_edge;}
                Face <T> *getFace() const {return face;}
                void setHalfEdge ( HalfEdge <T> *half_edge_ ) {half_edge = half_edge_;}
                void setFace ( Face <T> *face_ ) {face = face_;}

        public:
                //Other functions
                virtual void print ( std::ostream * file = &std::cout ) const;
                virtual Node3DCartesian <T> *clone() {return new Node3DCartesian <T> ( this->point_label, this->x, this->y, this->z, this->half_edge, this->face );}
};

#include "Node3DCartesian.hpp"

#endif
