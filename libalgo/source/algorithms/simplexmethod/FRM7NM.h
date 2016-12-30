// Description: Functor, compute residuals between sets P,P'
// Utilized for the map projection analysis
// Method  M7 (7 determined parameters, NM)
// Calls NLS operator ()

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


#ifndef FRM7NM_H
#define FRM7NM_H

#include "libalgo/source/types/TVector.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/projection2/Projection.h"

#include "libalgo/source/structures/matrix2/Matrix.h"

template <typename T>
class FRM7NM
{
private:
	const TVector <Point3DCartesian <T> > &test_points;		//List of test points
	const TVector <Point3DGeographic <T> > &reference_points;	//List of analyzed points
	const p_coord_function <T> getX, getY;				//Pointers to the coordinate functions
	const TTransformedLongitudeDirection trans_lon_dir;		//Transformed longitude direction
	T &dx, &dy;							//Shifts between analyzed and reference maps

public:
	FRM7NM(const TVector <Point3DCartesian <T> > &test_points_, const TVector <Point3DGeographic <T> > &reference_points_, const p_coord_function <T> &pX, const p_coord_function <T> &pY,
		const TTransformedLongitudeDirection trans_lon_dir_, T &dx_, T &dy_) : test_points(test_points_), reference_points(reference_points_), getX(pX), getY(pY),
		trans_lon_dir(trans_lon_dir_), dx(dx_), dy(dy_) {}


	void operator () (const Matrix <T> &XT, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W)
	{
		FRM7 <T> function (test_points, reference_points, getX, getY, trans_lon_dir, dx, dy);

		//Transpose X
		Matrix <T> X = MatrixOperations::trans(XT);

		//Call operator ()
		function(X, Y, V, W);
	}
};

#endif