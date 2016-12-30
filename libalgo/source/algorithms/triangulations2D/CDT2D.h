// Description: Constrained 2D Delaunay triangulation (Sloan, 1992)

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


#ifndef CDT2D_H
#define CDT2D_H

#include <list>

#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class Edge;


//New user types
template <typename T>
struct THalfEdgesCDTList
{
        typedef std::list < HalfEdge <T> *> Type;
};

//2D CDT class
class CDT2D
{
        public:
                template <typename T>
                static void CDT ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges, Container <Edge <Node3DCartesian <T> > *> &constrained_edges, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

        private:
                template <typename T>
                static void insertConstrainedEdge ( Edge <Node3DCartesian <T> > *e_constrained, HalfEdge <T> *he, Container <Edge <Node3DCartesian <T> > *> &constrained_edges );

                template <typename T>
                static void findIntersectedEdges ( Edge <Node3DCartesian <T> > *e_constrained, HalfEdge <T> *he, typename THalfEdgesCDTList <T>::Type & half_edges_removed, Container <Edge <Node3DCartesian <T> > *> &constrained_edges );

                template <typename T>
                static void removeIntersectedEdges ( const Edge <Node3DCartesian <T> > *e_constrained, typename THalfEdgesCDTList <T>::Type & half_edges_removed, typename THalfEdgesCDTList <T>::Type & half_edges_created );

                template <typename T>
                static void restoreDT ( const Edge <Node3DCartesian <T> > *e_constrained, typename THalfEdgesCDTList <T>::Type & half_edges_removed, typename THalfEdgesCDTList <T>::Type & half_edges_created );

                template <typename T>
                static bool legalizeConstrainedTriangle ( HalfEdge <T> *e_twin );

};

#include "CDT2D.hpp"

#endif
