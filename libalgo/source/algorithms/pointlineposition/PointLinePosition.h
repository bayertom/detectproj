// Description: Position of the point and line

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


#ifndef PointLinePosition_H
#define PointLinePosition_H

#include "libalgo/source/algorithms/vectorvectororientation/VectorVectorOrientation.h"

//Forward declaration
template <typename T>
class Point3DCartesian;

//Calculatre distance point and line
class PointLinePosition
{
        public:
                template <typename T>
                static unsigned short getPointLinePosition2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 );
};

#include "PointLinePosition.hpp"

#endif
