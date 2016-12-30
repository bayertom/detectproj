// Description: Sphere intersections

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


#ifndef SphereIntersection_H
#define SphereIntersection_H


class SphereIntersection
{
	public:
		template <typename T>
		static bool getSphereAndLineIntersection(const T xc, const T yc, const T zc, const T r, const T xa, const T ya, const T za,
			const T tx, const T ty, const T tz, T &xi1, T &yi1, T &zi1, T &xi2, T &yi2, T &zi2);
};

#include "SphereIntersection.hpp"

#endif