// Description: Compute angle given by 3 points

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


#ifndef Angle3Points_H
#define Angle3Points_H

#include "libalgo/source/structures/point/Point3DCartesian.h"

//Compute angle between ordered set of three points
class Angle3Points
{

        public:
                template <typename T>
                static T getAngle3Points ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3 );
};

#include "Angle3Points.hpp"

#endif
