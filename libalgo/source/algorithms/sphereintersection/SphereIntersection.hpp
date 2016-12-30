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


#ifndef SphereIntersection_HPP
#define SphereIntersection_HPP


template <typename T>
bool SphereIntersection::getSphereAndLineIntersection(const T xc, const T yc, const T zc, const T r, const T xa, const T ya, const T za,
	const T ux, const T uy, const T uz, T &xi1, T &yi1, T &zi1, T &xi2, T &yi2, T &zi2)
{
	//Compute intersection of the line and sphere
	//Sphere S(C, r) given by the center and radius, line L(x, y, z) = (xa, ya, ya) + t(ux, uy, uz)
	const T A = ux * ux + uy * uy + uz * uz;
	const T B = 2.0 * (xa * ux + ya * uy + za * uz - ux * xc - uy * yc - uz * zc);
	const T C = xa * xa - 2.0 * xa * xc + xc * xc + ya * ya - 2 * ya * yc + yc * yc + za * za - 2.0 * za * zc + zc * zc - r * r;

	//Compute discriminant
	const T D = B * B - 4 * A * C;

	//Intersection of the line and sphere does not exist
	if (D < 0)
		return false;

	//Find first intersection
	double t1 = (-B - sqrt(D)) / (2.0 * A);

	//Compute first intersection point
	xi1 = xa + t1 * ux;
	yi1 = ya + t1 * uy;
	zi1 = za + t1 * uz;

	//Line is a tangent, only 1 intersecion
	if (fabs(D) < EPS)
	{
		xi2 = xi1;
		yi2 = yi1;
		zi2 = zi1;

		return true;
	}

	//Find second intersection
	double t2 = (-B + sqrt(D)) / (2.0 * A);

	//Compute second intersection point
	xi2 = xa + t2 * ux;
	yi2 = ya + t2 * uy;
	zi2 = za + t2 * uz;

	//Both intersections have a different direction according to the line: switch points
	if (fabs(t1 - 0.5) >= fabs(t2 - 0.5))
	{
		//Create temporary variables
		const T x_temp = xi1;
		const T y_temp = yi1;
		const T z_temp = zi1;

		//Switch p1 <-> p2
		xi1 = xi2;
		yi1 = yi2;
		zi1 = zi2;

		xi2 = x_temp;
		yi2 = y_temp;
		zi2 = z_temp;
	}

	return true;
}

#endif