// Description: Functor, compute residuals between sets P,P'
// Utilized for the map projection analysis
// Method  M8S (8 determined parameters, scaled, a rotation involved, NLS)

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


#ifndef FRM8_H
#define FRM8_H

#include "libalgo/source/types/TVector.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/projection2/Projection.h"

#include "libalgo/source/structures/matrix2/Matrix.h"


template <typename T>
class FRM8
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
	FRM8(const TVector <Point3DCartesian <T> > &test_points_, const TVector <Point3DGeographic <T> > &reference_points_, const p_coord_function <T> &pX, const p_coord_function <T> &pY, 
		const TTransformedLongitudeDirection trans_lon_dir_, T &R_, T &q1_, T &q2_, T &dx_, T &dy_) : test_points(test_points_), 
		reference_points(reference_points_), getX(pX), getY(pY), trans_lon_dir(trans_lon_dir_), R(R_), q1(q1_), q2(q2_), dx(dx_), dy(dy_) {}

	void operator () (const Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W)
	{
		//Evaluate residuals between P, P', method M8, scaled
		const unsigned int m = test_points.size();

		//Apply projection proj(Q->P'): convert geographic points to the cartesian
		unsigned int index = 0;
		T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;
		TVector <Point3DCartesian <T> > reference_points_projected;

		for (const auto p : reference_points)
		{
			//(lat, lon) -> (lat_trans, lon_trans)
			Point3DGeographic <T> pole(X(0, 0), X(1, 0));
			const T lat_trans = CartTransformation::latToLatTrans(p.getLat(), p.getLon(), pole.getLat(), pole.getLon());
			const T lon_trans = CartTransformation::lonToLonTrans(p.getLat(), p.getLon(), pole.getLat(), pole.getLon(), trans_lon_dir);

			//Reduce longitude lon0 (not lon0_trans)
			const T lon_transr = CartTransformation::redLon0(lon_trans, X(4, 0));

			// (lat_trans, lon_trans) -> (X, Y)
			const T XR = getX(R, X(2, 0), X(3, 0), lat_trans, lon_transr, 0, 0, 0, X(5, 0));
			const T YR = getY(R, X(2, 0), X(3, 0), lat_trans, lon_transr, 0, 0, 0, X(5, 0));
			
			//Add point to the list
			Point3DCartesian <T> p_temp(XR, YR);
			reference_points_projected.push_back(p_temp);

			//Coordinate sums
			x_mass_test += test_points[index].getX();
			y_mass_test += test_points[index].getY();
			x_mass_reference += XR;
			y_mass_reference += YR;

			index++;
		}

		//Update weight matrix
		W = eye(2 * m, 2 * m, 1.0);

		//Compute centers of mass for both systems P, P'
		x_mass_test = x_mass_test / m;
		y_mass_test = y_mass_test / m;
		x_mass_reference = x_mass_reference / m;
		y_mass_reference = y_mass_reference / m;

		//Compute scale using the least squares adjustment: h = inv (A'WA)A'WL, 2D Helmert transformation
		T sum_xy_1 = 0, sum_xy_2 = 0, sum_xx_yy = 0;
		for (unsigned int i = 0; i < m; i++)
		{
			sum_xy_1 = sum_xy_1 + (test_points[i].getX() - x_mass_test) * (reference_points_projected[i].getX() - x_mass_reference) +
				(test_points[i].getY() - y_mass_test) * (reference_points_projected[i].getY() - y_mass_reference);
			sum_xy_2 = sum_xy_2 + (test_points[i].getY() - y_mass_test) * (reference_points_projected[i].getX() - x_mass_reference) -
				(test_points[i].getX() - x_mass_test) * (reference_points_projected[i].getY() - y_mass_reference);
			sum_xx_yy = sum_xx_yy + (reference_points_projected[i].getX() - x_mass_reference) * (reference_points_projected[i].getX() - x_mass_reference) +
				(reference_points_projected[i].getY() - y_mass_reference) * (reference_points_projected[i].getY() - y_mass_reference);
		}

		//Transformation ratios
		q1 = sum_xy_1 / sum_xx_yy;
		q2 = sum_xy_2 / sum_xx_yy;

		//Rotation
		const T alpha = atan2(q2, q1) * 180.0 / M_PI;

		//Compute coordinate differences (residuals): estimated - input
		for (unsigned int i = 0; i < m; i++)
		{
			V(i, 0) = (q1 * (reference_points_projected[i].getX() - x_mass_reference) - q2 * (reference_points_projected[i].getY() - y_mass_reference) - (test_points[i].getX() - x_mass_test));
			V(i + m, 0) = (q2 * (reference_points_projected[i].getX() - x_mass_reference) + q1 * (reference_points_projected[i].getY() - y_mass_reference) - (test_points[i].getY() - y_mass_test));
		}

		//Compute shifts dx, dy
		dx = x_mass_test - x_mass_reference * q1 + y_mass_reference * q2;
		dy = y_mass_test - x_mass_reference * q2 - y_mass_reference * q1;

		//Evaluate new radius
		R *= sqrt(q1 * q1 + q2 * q2);
	}
};

#endif