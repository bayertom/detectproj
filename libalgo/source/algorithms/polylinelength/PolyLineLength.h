// Description: Length of the polyline

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


#ifndef PolyLineLength_H
#define PolyLineLength_H


#include "libalgo/source/structures/list/Container.h"


//Forward declarations
template <typename T>
class Point3DCartesian;

template <typename T>
class Node3DCartesian;

class TIndexList;

//Compute length of the polyline
class PolyLineLength
{
        public:
                template <typename Point>
                static typename Point::Type getPolyLineLength ( const Container <Point> &points );

                template <typename T>
                static T getPolyLineLength ( const Container <Node3DCartesian <T> *> &points, const TIndexList & il );

};

#include "PolyLineLength.hpp"

#endif

