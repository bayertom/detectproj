// Description: Face overlay operations: union, difference, intersection, xor
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


#ifndef FaceOverlay_H
#define FaceOverlay_H

#include <vector>
#include <map>

#include "libalgo/source/structures/list/Container.h"


//Forward declaration
template <typename T>
class Face;

template <typename T>
class Node3DCartesian;


//Type of overlay operation
typedef enum
{
        Intersection = 0,
        Union,
        Difference,
} TOverlayType;


//Structure for storing intersections: intersections are stored according to their distance from start point
template <typename T>
struct TIntersectionsList
{
        typedef std::vector < std::map <T, Node3DCartesian <T> *> > Type;

        TIntersectionsList ( const unsigned int n ) {Type ( n ( std::map <T, Node3DCartesian <T> *> ) ); }
};


//Structure for storing segments of the same weight: inner / outer
template <typename T>
struct TSegmentsType
{
        typedef std::vector < std::vector <Node3DCartesian <T> *> > Type;
};


//List of processed points
typedef std::vector <bool> TProcessedPoints;


//Face overlay operations: union, difference, intersection, xor
class FaceOverlay
{
        public:
                template <typename T>
                static void createOverlay ( const Face <T> *f1, const Face <T> *f2, const TOverlayType overlay_type, Container <Face <T> *> &result, Container <Node3DCartesian <T> * > &intersections, Container <HalfEdge <T> *> &hl, bool & f1_simple, bool & f2_simple );

        private:
                template <typename T>
                static void computeLineSegmentsIntersections ( Container <Node3DCartesian <T> *, NonDestructable > &f1_nodes, Container <Node3DCartesian <T> *, NonDestructable > &f2_nodes, typename TIntersectionsList <T>::Type & il1,
                                typename TIntersectionsList <T>::Type & il2, Container <Node3DCartesian <T> * > &intersections );

                template <typename T>
                static void setWeights ( const Container <Node3DCartesian <T> *, NonDestructable > &f1_nodes, Container <Node3DCartesian <T> *, NonDestructable > &f2_nodes, TIndexList & weights );

                template <typename T>
                static void splitSegmentsByWeights ( const Container <Node3DCartesian <T> *, NonDestructable > &f_nodes, const TIndexList &weights, typename TSegmentsType <T>::Type & inner,  typename TSegmentsType <T>::Type & outer );

                template <typename T >
                static void createFaceOverlay ( typename TSegmentsType <T>::Type & f1,  typename TSegmentsType <T>::Type & f2, const TOverlayType & overlay_type, Container <Face <T> *> &result, Container <HalfEdge <T> *> &hl );

};

#include "FaceOverlay.hpp"

#endif
