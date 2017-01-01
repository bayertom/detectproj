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


#ifndef GreatCircleIntersection_HPP
#define GreatCircleIntersection_HPP


#include "libalgo/source/algorithms/planeintersection/PlaneIntersection.h"
#include "libalgo/source/algorithms/sphereintersection/SphereIntersection.h"

template <typename T>
bool GreatCircleIntersection::getGreatCirclePlainIntersection(const Point3DGeographic<T> &p1, const Point3DGeographic<T> &p2, const Point3DGeographic<T> &p3, const Point3DGeographic<T> &p4,
	const Point3DGeographic<T> &p5, const Point3DGeographic<T> &p6, Point3DGeographic<T> &i1, Point3DGeographic<T> &i2, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection &lon_direction)
{
	//Compute intersection of the great circle given by the points p1, p2, p3 with a plane (meridian/parallel) given by points p4, p5, p6
	const T lat1 = p1.getLat();
	const T lon1 = p1.getLon();
	const T lat2 = p2.getLat();
	const T lon2 = p2.getLon();
	const T lat3 = p3.getLat();
	const T lon3 = p3.getLon();
	const T lat4 = p4.getLat();
	const T lon4 = p4.getLon();
	const T lat5 = p5.getLat();
	const T lon5 = p5.getLon();
	const T lat6 = p6.getLat();
	const T lon6 = p6.getLon();

	//P1: Convert spherical coordinates to the Cartesian
	const T x1 = cos(lat1 / RO) * cos(lon1 / RO);
	const T y1 = cos(lat1 / RO) * sin(lon1 / RO);
	const T z1 = sin(lat1 / RO);

	//P2: Convert spherical coordinates to the Cartesian
	const T x2 = cos(lat2 / RO) * cos(lon2 / RO);
	const T y2 = cos(lat2 / RO) * sin(lon2 / RO);
	const T z2 = sin(lat2 / RO);

	//P3: Convert spherical coordinates to the Cartesian
	const T x3 = cos(lat3 / RO) * cos(lon3 / RO);
	const T y3 = cos(lat3 / RO) * sin(lon3 / RO);
	const T z3 = sin(lat3 / RO);

	//P4: Convert spherical coordinates to the Cartesian
	const T x4 = cos(lat4 / RO) * cos(lon4 / RO);
	const T y4 = cos(lat4 / RO) * sin(lon4 / RO);
	const T z4 = sin(lat4 / RO);

	//P5: Convert spherical coordinates to the Cartesian
	const T x5 = cos(lat5 / RO) * cos(lon5 / RO);
	const T y5 = cos(lat5 / RO) * sin(lon5 / RO);
	const T z5 = sin(lat5 / RO);

	//P6: Convert spherical coordinates to the Cartesian
	const T x6 = cos(lat6 / RO) * cos(lon6 / RO);
	const T y6 = cos(lat6 / RO) * sin(lon6 / RO);
	const T z6 = sin(lat6 / RO);

	//Create initial point p0
	const T x0 = (x1 + x2 + x3) / 3.0;
	const T y0 = (y1 + y2 + y3) / 3.0;
	const T z0 = (z1 + z2 + z3) / 3.0;

	//Get parametric equation of the planes (p1, p2, p3) x (p4, p5, p6) intersection
	T xi, yi, zi, ux, uy, uz;
	const bool intersection1_exists = PlaneIntersection::get2PlanesIntersection(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x0, y0, z0, xi, yi, zi, ux, uy, uz);
	
	//Get coordinates of the intersection (sphere x plane intersection) in the horizontal / vertical plane
	if (intersection1_exists)
	{
		//Compute intersection of the sphere and both planes intersection
		const T xc = 0.0, yc = 0.0, zc = 0.0, r = 1.0;
		T xi1, yi1, zi1, xi2, yi2, zi2;
		
		const bool intersection2_exists = SphereIntersection::getSphereAndLineIntersection(xc, yc, zc, r, xi, yi, zi, ux, uy, uz, xi1, yi1, zi1, xi2, yi2, zi2);
		
		if (intersection2_exists)
		{
			//Convert Cartesian coordinates to the spherical
			const T ri1 = sqrt(xi1 * xi1 + yi1 * yi1 + zi1 * zi1);
			const T ri2 = sqrt(xi2 * xi2 + yi2 * yi2 + zi2 * zi2);

			//Compute latitude of intersections
			const T lati1 = asin(zi1 / ri1) * RO;
			const T lati2 = asin(zi2 / ri2) * RO;

			//Compute longitude of intersections
			const T loni1 = atan2(yi1, xi1) * RO;
			const T loni2 = atan2(yi2, xi2) * RO;

			//Set parameters to points
			i1.setLat(lati1);
			i1.setLon(loni1);
			i2.setLat(lati2);
			i2.setLon(loni2);

			return true;
		}
	}

	return false;
}


#endif