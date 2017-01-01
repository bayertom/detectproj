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


#ifndef PlaneIntersection_HPP
#define PlaneIntersection_HPP

#include "libalgo/source/algorithms/matrixoperations2/MatrixOperations.h"

//Set namespace
using namespace MatrixOperations;

template <typename T>
bool PlaneIntersection::get2PlanesIntersection(const T &x1, const T &y1, const T &z1, const T &x2, const T &y2, const T &z2, const T &x3, const T &y3, const T &z3,
	const T &x4, const T &y4, const T &z4, const T &x5, const T &y5, const T &z5, const T &x6, const T &y6, const T &z6, const T &x0, const T &y0, const T &z0,
	T &xi, T &yi, T &zi, T &ux, T &uy, T &uz)
{
	//Find intersection of two planes Rho(P1, P2, P3), Sigma(P4, P5, P6) given by points
	//P0 =[x0, y0, z0] is the arbitrary point; start point of the intersection Pi=[xi, yi, zi] determined as closest to P0
	//Intersection: L(x, y, z) = (xi, yi, zi) + t(nx, ny, nz)
	//Solution by John Krumm, 2016

	//Compute direction vector u1 = p2 - p1  of the first plane
	const T u1x = x2 - x1;
	const T u1y = y2 - y1;
	const T u1z = z2 - z1;

	//Compute direction vector v1 = p3 - p1 of the first plane
	const T v1x = x3 - x1;
	const T v1y = y3 - y1;
	const T v1z = z3 - z1;

	//Cross product n1 = u1 x v1 of the first plane
	T n1x = u1y * v1z - v1y * u1z;
	T n1y = u1z * v1x - v1z * u1x;
	T n1z = u1x * v1y - v1x * u1y;
	const T n1_norm = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);

	//Normalize n1
	n1x /= n1_norm;
	n1y /= n1_norm;
	n1z /= n1_norm;

	//Compute direction vector u2 = p5 - p4  of the second plane
	const T u2x = x5 - x4;
	const T u2y = y5 - y4;
	const T u2z = z5 - z4;

	//Compute direction vector v2 = p6 - p4 of the second plane
	const T v2x = x6 - x4;
	const T v2y = y6 - y4;
	const T v2z = z6 - z4;

	//Cross product n1 = u1 x v1 of the second plane
	T n2x = u2y * v2z - v2y * u2z;
	T n2y = u2z * v2x - v2z * u2x;
	T n2z = u2x * v2y - v2x * u2y;
	const T n2_norm = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);

	//Normalize n2
	n2x /= n2_norm;
	n2y /= n2_norm;
	n2z /= n2_norm;

	//Create matrix M
	Matrix <T> M(5, 5);
	M(0, 0) = 2.0;
	M(0, 3) = n1x;
	M(0, 4) = n2x;
	
	M(1, 1) = 2.0;
	M(1, 3) = n1y;
	M(1, 4) = n2y;
	
	M(2, 2) = 2.0;
	M(2, 3) = n1z;
	M(2, 4) = n2z;

	M(3, 0) = n1x;
	M(3, 1) = n1y;
	M(3, 2) = n1z;

	M(4, 0) = n2x;
	M(4, 1) = n2y;
	M(4, 2) = n2z;

	//Test, if planes are colinear
	const T detM = det(M);

	//No intersection found
	if (detM < EPS)
		return false;

	//Create matrix B
	Matrix <T> b(5, 1);
	b(0, 0) = 2 * x0;
	b(1, 0) = 2 * y0;
	b(2, 0) = 2 * z0;
	b(3, 0) = x1 * n1x + y1 * n1y + z1 * n1z;
	b(4, 0) = x4 * n2x + y4 * n2y + z4 * n2z;

	//Find solution Mx = b;
	Matrix <T> x = inv(M) * b;

	//Cross product u = n1 x n2: vector of the intersection
	ux = n1y * n2z - n2y * n1z;
	uy = n1z * n2x - n2z * n1x;
	uz = n1x * n2y - n2x * n1y;
	

	//Start point of the intersection
	xi = x(0, 0);
	yi = x(1, 0);
	zi = x(2, 0);

	return true;
}


#endif