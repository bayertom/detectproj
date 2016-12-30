// Description: 2D Delaunay triangle by incremental insertion + Lawson Oriented Walk, topological model, de Berg et al, 2000

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


#ifndef DT2D_H
#define DT2D_H

#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class HalfEdge;

//2D Delaunay triangulation
class DT2D
{
        public:
                template <typename T>
                static void DT ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges_dt, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static void swapDiagonal ( HalfEdge <T> *e_twin, HalfEdge <T> *e12, HalfEdge <T> *e13, HalfEdge <T> *e21, HalfEdge <T> *e22, HalfEdge <T> *e23 );

        private:

                template <typename T>
                static void createSimplexTriangle ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges_dt );

                template <typename T>
                static void DTInsertPoint ( Node3DCartesian <T> *p, const Node3DCartesian <T> *s1, const Node3DCartesian <T> *s2, const Node3DCartesian <T> *s3, HalfEdge <T> **e1, Container <HalfEdge <T> *> &half_edges_dt );

                template <typename T>
                static void legalizeTriangle ( Node3DCartesian <T> *p, HalfEdge <T> *e11, const Node3DCartesian <T> *s1, const Node3DCartesian <T> *s2, const Node3DCartesian <T> *s3 );

                template <typename T>
                static void removeSimplexTriangles ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges_dt );
};

#include "DT2D.hpp"

#endif
