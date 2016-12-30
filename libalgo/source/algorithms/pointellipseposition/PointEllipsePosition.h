// Description: Position and ellipse rotated by the angle phi (counterclockwise direction)

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


#ifndef PointEllipsePosition_H
#define PointEllispePosition_H


//Forward declaration
template <typename T>
class Point3DCartesian;

//Position of the point and ellipse
class PointEllipsePosition
{
        public:
                template <typename T>
                static unsigned short getPointEllipsePosition ( const T x, const T y, const T xc, const T yc, const T a, const T b, const T angle );
};

#include "PointEllipsePosition.hpp"

#endif
