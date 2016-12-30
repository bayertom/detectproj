// Description: Functor, compute residuals between sets P,P'
// Utilized for the map projection analysis
// Method  M8S (8 determined parameters, scaled, a rotation involved, NM)
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


#ifndef FRM8NM_H
#define FRM8NM_H

#include "libalgo/source/types/TVector.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/projection2/Projection.h"

#include "libalgo/source/structures/matrix2/Matrix.h"


template <typename T>
class FRM8NM
{
private:
	const TVector <Point3DCartesian <T> > &test_points;		//List of test points
	const TVector <Point3DGeographic <T> > &reference_points;	//List of analyzed points
	const p_coord_function <T> getX, getY;				//Pointers to the coordinate functions
	const TTransformedLongitudeDirection trans_lon_dir;		//Transformed longitude direction
	T &R;								//Earth radius (will be updated)
	T &q1, &q2;							//Coefficient of  2D Helmert transformation (will be updated)
	T &dx, &dy;							//Shifts between analyzed and reference maps

public:
	FRM8NM(const TVector <Point3DCartesian <T> > &test_points_, const TVector <Point3DGeographic <T> > &reference_points_, const p_coord_function <T> &pX, const p_coord_function <T> &pY,
		const TTransformedLongitudeDirection trans_lon_dir_, T &R_, T &q1_, T &q2_, T &dx_, T &dy_) : test_points(test_points_),
		reference_points(reference_points_), getX(pX), getY(pY), trans_lon_dir(trans_lon_dir_), R(R_), q1(q1_), q2(q2_), dx(dx_), dy(dy_) {}

	void operator () (const Matrix <T> &XT, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W)
	{
		//Call NLSP functor
		FRM8 <T> function(test_points, reference_points, getX, getY, trans_lon_dir, R, q1, q2, dx, dy);
		
		//Transpose X
		Matrix <T> X = MatrixOperations::trans(XT);
		
		//Call operator ()
		function( X, Y, V, W);
		
	}


};

#endif