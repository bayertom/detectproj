// Description: Plane intersections

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


#ifndef PlaneIntersection_H
#define PlaneIntersection_H


class PlaneIntersection
{
	public:
		template <typename T>
		static bool get2PlanesIntersection(const T &x1, const T &y1, const T &z1, const T &x2, const T &y2, const T &z2, const T &x3, const T &y3, const T &z3,
			const T &x4, const T &y4, const T &z4, const T &x5, const T &y5, const T &z5, const T &x6, const T &y6, const T &z6, const T &x0, const T &y0, const T &z0,
			T &xi, T &yi, T &zi, T &ux, T &uy, T &uz);
};

#include "PlaneIntersection.hpp"

#endif