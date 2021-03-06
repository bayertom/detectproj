// Description: Functor, numerical differentiation of the projection equations according to R, latp, lonp, lat1, lat2, lon0, c (Method M8S, NLSP)
// Method: Stirling numerical differentiation

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


#ifndef FDiffM7_H
#define FDiffM7_H

#include "libalgo/source/types/PCoordFunction.h"

#include "libalgo/source/structures/matrix2/Matrix.h"

#include "libalgo/source/algorithms/carttransformation2/CartTransformation.h"

//Functor, numerical differentiation of the projection equations according to [R, latp, lonp, lat1, lat2, lon0, c] using the Stirling formula
template <typename T>
class FDiffM7
{
	private:
		const T lat;						//Latitude of the point
		const T lon;						//Longitude of the point
		const p_coord_function <T> equat;			//Coordinate function
		const TTransformedLongitudeDirection trans_lon_dir;	//Transformed longitude direction


	public:
		FDiffM7(const T lat_, const T lon_, const p_coord_function <T> &equat_, const TTransformedLongitudeDirection trans_lon_dir_) :
			lat(lat_), lon(lon_), equat(equat_), trans_lon_dir(trans_lon_dir_) {}

		T operator () (const Matrix <T> &XT)
		{
			//Convert ( lat, lon ) -> ( lat_trans, lon_trans)_trans
			const T lat_trans = CartTransformation::latToLatTrans(lat, lon, XT(0, 1), XT(0, 2));
			const T lon_trans = CartTransformation::lonToLonTrans(lat, lon, XT(0, 1), XT(0, 2), trans_lon_dir);

			//Reduce longitude lon0 (not lon0_trans)
			const T lon_transr = CartTransformation::redLon0(lon_trans, XT(0, 5));

			//Evaluate map projection equation; used for the partial derivative
			T res = equat(XT(0, 0), XT(0, 3), XT(0, 4), lat_trans, lon_transr, 0, 0, 0, XT(0, 6));
			
			return res;
		}
};

#endif