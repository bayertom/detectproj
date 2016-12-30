// Description: transformation of the geographic coordinates (lat, lon) to (lat, lon)_trans.

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

#ifndef CartTransformation_H
#define CartTransformation_H

#include <cmath>
#include <memory>

#include "libalgo/source/types/TTransformedLongitudeDirection.h"
#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/PCoordFunction.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"

#include "libalgo/source/structures/projection2/Projection.h"

//Basic cartographic transformations
class CartTransformation
{
public:

	template <typename T>
	static T latToLatTrans(const Point3DGeographic <T> &p, const Point3DGeographic <T> &pole);

	template <typename T>
	static T latTransToLat(const Point3DGeographic <T> &p, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection &lon_direction);

	template <typename T>
	static T lonToLonTrans(const Point3DGeographic <T> &p, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection &lon_direction = NormalDirection);

	template <typename T>
	static T lonTransToLon(const Point3DGeographic <T> &p, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection &lon_direction = NormalDirection2);

	template <typename T>
	static T redLon0(const T lon, const T lon0) { return (lon - lon0 < MIN_LON ? 360.0 + (lon - lon0) : (lon - lon0 > MAX_LON ? (lon - lon0) - 360 : lon - lon0)); }

	template <typename T>
	static T latToLatTrans(const T lat, const T lon, const T latp, const T lonp);

	template <typename T>
	static T latTransToLat(const T lat_trans, const T lon_trans, const T latp, const T lonp, const TTransformedLongitudeDirection &lon_direction);

	template <typename T>
	static T lonToLonTrans(const T lat, const T lon, const T latp, const T lonp, const TTransformedLongitudeDirection &lon_direction = NormalDirection2);

	template <typename T>
	static T lonTransToLon(const T lat_trans, const T lon_trans, const T latp, const T lonp, const TTransformedLongitudeDirection &lon_direction);

	template <typename T>
	static void latLontoXY(const T lat, const T lon, const std::shared_ptr<Projection<T> > &proj, const T alpha, T &X, T &Y);
	
	template <typename T>
	static void latLontoXY(const T R, const T  lat1, const T lat2, const T lat, const T lon, const T latp, const T lonp, const TTransformedLongitudeDirection lon_dir, const T lon0, const T dx, const T dy, const T c, const p_coord_function <T> getX, const p_coord_function <T> getY, const T alpha, T &X, T &Y);

	template <typename T>
	static void latsLonstoXY(const TVector <Point3DGeographic<T>> reference_points, const std::shared_ptr<Projection<T> > &proj, const T alpha, TVector <Point3DCartesian <T>> &projected_points);

	template <typename T>
	static void wgsToJTSK(const Point3DGeographic <T> &p1, Point3DCartesian <T> & p2);
};

#include "CartTransformation.hpp"

#endif
