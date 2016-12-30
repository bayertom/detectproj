// Description: Position of two lines

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


#ifndef LineLinePosition_H
#define LineLinePosition_H


#include "libalgo/source/algorithms/pointlineposition/PointLinePosition.h"


//Forward declaration
template <typename T>
class Point3DCartesian;


//Position/intersection of 2 lines
class LineLinePosition
{
        public:
                template <typename T>
                static unsigned short get2LineSegmentsPosition ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, const Point3DCartesian <T> * p4 , T & x_int, T & y_int );

                template <typename T>
                static unsigned short get2LinesPosition ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, const Point3DCartesian <T> * p4 , T & x_int, T & y_int );

                template <typename T>
                static unsigned short get2LineSegmentsPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int );

                template <typename T>
                static unsigned short get2ParallelsPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int );

                template <typename T>
                static bool get2LinesPosition ( const T x1, const T y1, const T x2, const T y2, const T x3, const T y3, const T x4, const T y4, T & x_int, T & y_int );
};

#include "LineLinePosition.hpp"

#endif
