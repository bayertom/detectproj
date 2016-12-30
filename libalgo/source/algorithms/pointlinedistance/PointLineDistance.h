// Description: Distance point-line

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


#ifndef PointLineDistance_H
#define PointLineDistance_H

//Forward declaration
template <typename T>
class Point3DCartesian;

//Point line distance
class PointLineDistance
{
        public:
                template <typename T>
                static T getPointLineDistance2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 );

                template <typename T>
                static T getPointLineDistance2D ( const T xa, const T ya, const T x1, const T y1, const T x2, const T y2 );

                template <typename T>
                static T getPointLineDistance2DSigned ( const T xa, const T ya, const T x1, const T y1, const T x2, const T y2 );

                template <typename T>
                static T getPointLineSegmentDistance2D ( const Point3DCartesian <T> * p, const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2 );

};

#include "PointLineDistance.hpp"

#endif
