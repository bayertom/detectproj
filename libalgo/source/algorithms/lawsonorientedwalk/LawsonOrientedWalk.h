// Description: Lawson oriented walk, find face containing point p using vector orientation test

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


#ifndef LawsonOrientedWalk_H
#define LawsonOrientedWalk_H

#include "libalgo/source/algorithms/pointlineposition/PointLinePosition.h"


//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class HalfEdge;

template <typename T>
class Face;


//Lawson oriented walk: find incident triangle / Face using walking
class LawsonOrientedWalk
{
        public:
                template <typename T>
                static HalfEdge <T> *findTriangleWalk ( const Node3DCartesian <T> *p, short * status, HalfEdge <T> *e );

                template <typename T>
                static Face <T> *findFaceWalk ( const Node3DCartesian <T> *p, short * status, Face <T> *pol );

                template <typename T>
                static HalfEdge <T> *findFaceWalk2 ( const Node3DCartesian <T> *p, short * status, HalfEdge <T> *e, const T max_steps );

};

#include "LawsonOrientedWalk.hpp"

#endif

