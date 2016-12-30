// Description: Parameters of the circle given by three points

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


#ifndef Circle3Points_H
#define Circle3Points_H

//Forward declaration
template <typename T>
class Point3DCartesian;

//Compute radius and mid point of the circle given by 3 points
class Circle3Points
{
        public:
                template <typename T>
                static void getCentreAndDiameterCircle ( const Point3DCartesian <T> * p1, const Point3DCartesian <T> * p2, const Point3DCartesian <T> * p3, T & xc, T & yc, T & r );
};

#include "Circle3Points.hpp"

#endif
