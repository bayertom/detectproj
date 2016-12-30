// Description: Test if a point is betwwen 2 points

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


#ifndef PointBetweenPoints_H
#define PointBetweenPoints_H


class PointBetweenPoints
{
        public:

                template <typename T>
                static bool isPointBetweenPoints ( const T xp, const T yp, const T x1, const T y1, const T x2, const T y2 );

};

#include "PointBetweenPoints.hpp"

#endif
