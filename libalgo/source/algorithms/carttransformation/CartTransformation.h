// Description: transformation of the geographic coordinates (lat, lon) to (lat, lon)_trans. Conversion (lat, lon) to (x, y).

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

#ifndef CartTransformation_H
#define CartTransformation_H

#include <math.h>

#include "libalgo/source/types/TTransformedLongitudeDirection.h"

#include "libalgo/source/const/Const.h"

#include "libalgo/source/algorithms/arithmeticparser/ArithmeticParser.h"

//Forward declaration
template <typename T>
class Point3DCartesian;

template <typename T>
class Point3DGeographic;

template <typename T>
class Projection;

//Basic cartographic transformations
class CartTransformation
{
public:

	template <typename T>
	static T latToLatTrans(const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole);

	template <typename T>
	static T lonToLonTrans(const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole, const TTransformedLongitudeDirection lon_direction = NormalDirection);

	template <typename T>
	static T latLonToX(const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception = true);

	template <typename T>
	static T latLonToY(const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception = true);

	template <typename T>
	static T redLon0(const T lon, const T lon0) { return (lon - lon0 < MIN_LON ? 360.0 + (lon - lon0) : (lon - lon0 > MAX_LON ? (lon - lon0) - 360 : lon - lon0)); }

	template <typename T>
	static T latToLatTrans(const T lat, const T lon, const T latp, const T lonp);

	template <typename T>
	static T lonToLonTrans(const T lat, const T lon, const T latp, const T lonp, const TTransformedLongitudeDirection lon_direction = NormalDirection);

	template <typename T>
	static T latLonToX(TPostfixNotationDel * equation_x_postfix, TPostfixNotationDel * ftheta_equat_postfix, TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon,
		const T R, const T a, const T b, const T dx, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception = true);

	template <typename T>
	static T latLonToY(TPostfixNotationDel * equation_y_postfix, TPostfixNotationDel * ftheta_equat_postfix, TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon,
		const T R, const T a, const T b, const T dy, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception = true);

	template <typename T>
	static T latLonToCartesian(const TPostfixNotationDel * equation_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat,
		const T lon, const T R, const T a, const T b, const T shift, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception = true);

	template <typename T>
	static void wgsToJTSK(const Point3DGeographic <T> *p1, Point3DCartesian <T> * p2);
};

#include "CartTransformation.hpp"

#endif
