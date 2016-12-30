// Description: Compute intersection of the great circle and the plane

// Copyright (c) 2015 - 2016
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


#ifndef GreatCircleIntersection_H
#define GreatCircleIntersection_H

#include "libalgo/source/types/TTransformedLongitudeDirection.h"

#include "libalgo/source/structures/point/Point3DGeographic.h"

class GreatCircleIntersection
{
	public:
		template <typename T>
		static bool getGreatCirclePlainIntersection(const Point3DGeographic<T> &p1, const Point3DGeographic<T> &p2, const Point3DGeographic<T> &p3, const Point3DGeographic<T> &p4, 
			const Point3DGeographic<T> &p5, const Point3DGeographic<T> &p6, Point3DGeographic<T> &i1, Point3DGeographic<T> &i2, const Point3DGeographic <T> &pole, 
			const TTransformedLongitudeDirection &lon_direction);

};

#include "GreatCircleIntersection.hpp"

#endif