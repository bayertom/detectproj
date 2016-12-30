// Description: 2D Helmert transformation

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


#ifndef HelmertTransformation2D_HPP
#define HelmertTransformation2D_HPP

#include <algorithm>
#include <cmath>

#include "libalgo/source/algorithms/matrixoperations2/MatrixOperations.h"

#include "libalgo/source/exceptions/BadDataException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"


using namespace MatrixOperations;

template <typename T, typename Point1, typename Point2, typename Point3>
void HelmertTransformation2D::transformPoints (const TVector <Point1> &global_points, const TVector <Point2> &local_points, TVector <Point3> &transformed_points, TTransformationKeyHelmert2D <T> &key_helmert )
{
        //Compute non weighted 2D Helmert transformation
        TVector <T> weights ( global_points.size(), 1.0 );
	getTransformKey(global_points, local_points, weights, key_helmert);
	transform(global_points, local_points, transformed_points, key_helmert);
}


template <typename T, typename Point1, typename Point2>
void HelmertTransformation2D::getTransformKey (const TVector <Point1> &global_points, const TVector <Point2> &local_points, const TVector <T> weights, TTransformationKeyHelmert2D <T> &key_helmert )
{
        //Get transformation key for weighted transformation
        const unsigned int n_global = global_points.size(), n_local = local_points.size();
        T sumx_local = 0, sumy_local = 0, sumx_global = 0, sumy_global = 0;

        //Not enough points
        if ( ( n_global < 2 ) || ( n_local < 2 ) )
        {
                throw  BadDataException ( "BadDataException: not enough points. ", "Can not compute Helmert 2D transformation key. \n" );
        }

        //Less local points
        if ( n_global > n_local )
        {
                throw  BadDataException ( "BadDataException: less local points than global points. ", "Can not compute Helmert 2D transformation key. \n" );
        }

        //Compute sums of coordinates, weights
        T sum_weights = 0;

        for ( unsigned int i = 0; i < n_global; i++ )
        {
                sumx_local += weights[i] * local_points [i].getX();
                sumy_local += weights[i] * local_points [i].getY();
                sumx_global += weights[i] * global_points [i].getX();
                sumy_global += weights[i] * global_points [i].getY();

                sum_weights += weights[i];
        }

        //Compute center of mass
        key_helmert.x_mass_local = sumx_local / ( sum_weights );
        key_helmert.y_mass_local = sumy_local / ( sum_weights );
        key_helmert.x_mass_global = sumx_global / ( sum_weights );
        key_helmert.y_mass_global = sumy_global / ( sum_weights );

        //Remeber k
        key_helmert.k = sum_weights;

        //Reduction of coordinates to the center of mass
        T x_red_local, y_red_local, xred_global, yred_global, k1 = 0, k2 = 0;

        //Process all points
        key_helmert.J = 0;

        for ( unsigned int i = 0; i < n_global; i++ )
        {
                //Compute reduced coordinates
                x_red_local = local_points [i].getX() - key_helmert.x_mass_local;
                y_red_local = local_points [i].getY() - key_helmert.y_mass_local;
                xred_global = global_points [i].getX() - key_helmert.x_mass_global;
                yred_global = global_points [i].getY() - key_helmert.y_mass_global;

                //Compute coefficients of transformation
                key_helmert.J += weights[i] * ( x_red_local * x_red_local + y_red_local * y_red_local );
                k1 += weights[i] * ( xred_global * x_red_local + yred_global * y_red_local );
                k2 += weights[i] * ( yred_global * x_red_local - xred_global * y_red_local );
        }

        //Throw exception
        if ( key_helmert.J == 0 )
        {
                throw  MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute Helmert 2D transformation key, ", " divider = 0.", key_helmert.J );
        }

        //Transformation coefficients
        key_helmert.c1 = k1 / key_helmert.J;
        key_helmert.c2 = k2 / key_helmert.J;
}


template <typename T, typename Point1, typename Point2, typename Point3>
void HelmertTransformation2D::transform (const TVector <Point1> &global_points, const TVector <Point2> &local_points, TVector <Point3> &transformed_points, const TTransformationKeyHelmert2D <T> & key_helmert )

{
        //Transform all points using 2D Helmert transformation ( weighted / non-weighted )
	const unsigned int n = local_points.size();

        if ( transformed_points.size() != 0 )
        {
                throw BadDataException ( "BadDataException: list of tranformed points is not empty. ", "Can not compute Helmert 2D transformation." );
        }

        for ( unsigned int i = 0; i < n; i++ )
        {
                //Reduce coordinates
                const T x_red_local = local_points [i].getX() - key_helmert.x_mass_local,
                                            y_red_local = local_points [i].getY() - key_helmert.y_mass_local;

                //Transform point, add coordinates center of mass
                const T x_transform = key_helmert.c1 * x_red_local - key_helmert.c2 * y_red_local + key_helmert.x_mass_global,
                                            y_transform = key_helmert.c2 * x_red_local + key_helmert.c1 * y_red_local + key_helmert.y_mass_global;

                //Add point to the list
                transformed_points.push_back ( Point3 ( x_transform, y_transform ) );
        }
}

/*
template <typename T>
Matrix <T> HelmertTransformation2D::getTransformKey(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W)
{
	//Helmert transformation, matrix notation
	const unsigned int m = P.rows(), n = P.cols();

	//Compute centers of gravity

	//Set elements of A(2m, 4) matrix, applicable to each scheme
	Matrix <T> A(2 * m, 4);
	A(Q(0, m - 1, 0, 0), 0, 0);
	A(Q(0, m - 1, 1, 1), m, 0);

	A(-1.0 * Q(0, m - 1, 1, 1), 0, 1);
	A(Q(0, m - 1, 0, 0), m, 1);

	Matrix <T> I = ones(m, 1, 1.0);
	A(I, 0, 2);
	A(I, m, 3);

	//Create Y vector
	Matrix <T> Y(2 * m, 1);
	Y(P(0, m - 1, 0, 0), 0, 0);
	Y(P(0, m - 1, 1, 1), m, 0);

	//Find solution
	Matrix <T> beta = inv(trans(A) * W * A) * trans(A) * W  * Y;

	return beta;
}


template <typename T>
void HelmertTransformation2D::getTransformKey2(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W, Matrix <T> &A, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &C)
{
	//Helmert transformation, matrix notation
	const unsigned int m = P.rows(), n = P.cols();

	//Get matrices
	const Matrix <T> XG = P(0, m - 1, 0, 0);
	const Matrix <T> YG = P(0, m - 1, 1, 1);
	const Matrix <T> XL = Q(0, m - 1, 0, 0);
	const Matrix <T> YL = Q(0, m - 1, 1, 1);
	const Matrix <T> WX = W(0, m - 1, 0, m - 1);
	const Matrix <T> WY = W(m, 2 * m - 1, m, 2 * m - 1);
	
	//Sum of weights
	const T swx = sum(diag(WX));
	const T swy = sum(diag(WY));

	//Centers of gravity
	const T x_mass_glob = sum( trans( XG ) * WX ) / swx;
	const T y_mass_glob = sum( trans( YG ) * WY ) / swy;
	const T x_mass_loc = sum( trans(XL) * WX ) / swx;
	const T y_mass_loc = sum( trans(YL) * WY ) / swy;

	//Store coordinates
	C(0, 0) = x_mass_glob;
	C(0, 1) = y_mass_glob;
	C(1, 0) = x_mass_loc;
	C(1, 1) = y_mass_loc;

	//Reduce coordinates
	const Matrix <T> XLR = XL - x_mass_loc;
	const Matrix <T> YLR = YL - y_mass_loc;

	//Set elements of A(2m, 4) matrix, applicable to each scheme
	A(XLR, 0, 0);
	A(YLR, m, 0);

	A(-1.0 * YLR, 0, 1);
	A(XLR, m, 1);

	Matrix <T> I = ones(m, 1, 1.0);
	A(I, 0, 2);
	A(I, m, 3);

	//Create Y vector
	Y(XG, 0, 0);
	Y(YG, m, 0);

	//Find solution
	X = inv(trans(A) * W * A) * trans(A) * W  * Y;

}
*/

#endif
