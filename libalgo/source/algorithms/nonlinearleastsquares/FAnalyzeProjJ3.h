// Description: Functor, create Jacobi matrix J for cartometric analysis, method M7 (7 determined parameters)
// Elements are computed numerically using the Stirling method

// Copyright (c) 2010 - 2014
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


#ifndef FAnalyzeProjJ3_H
#define FAnalyzeProjJ3_H

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative6Var.h"

//Forward declarations
template <typename T>
class Projection;

template <typename T>
class Node3DCartesian;

template <typename T>
class Point3DGeographic;

template <typename T>
class FAnalyzeProjJ3
{
private:

	//Center of mass
	const T &x_mass_reference;
	const T &y_mass_reference;

	//List of points
	const Container <Node3DCartesian <T> *> &nl_test;
	const Container <Point3DGeographic <T> *> &pl_reference;
	const Container <Node3DCartesianProjected <T> *> &nl_projected;

	//Map projection, analyzed aspect, method
	const Projection <T> *proj;
	const TProjectionAspect aspect;

	const bool print_exceptions;

	unsigned int &iter;

public:

	FAnalyzeProjJ3(const Container <Node3DCartesian <T> *> &nl_test_, const Container <Point3DGeographic <T> *> &pl_reference_, const Container <Node3DCartesianProjected <T> *> &nl_projected_, Projection <T> *proj_, const TProjectionAspect aspect_, const T &x_mass_reference_, const T &y_mass_reference_, const bool print_exceptions_, unsigned int &iter_)
		: nl_test(nl_test_), pl_reference(pl_reference_), nl_projected(nl_projected_), proj(proj_), aspect(aspect_), x_mass_reference(x_mass_reference_), y_mass_reference(y_mass_reference_), print_exceptions(print_exceptions_), iter (iter_) {}


	void operator () (const Matrix <T> &X, Matrix <T> &J )
	{
		//Compute parameters of the Jacobi Matrix J
		//Jacobi matrix J = [ d_R, d_latp, d_lonp, d_lat0, d_lon0, c]
		const unsigned int m = nl_test.size();
		unsigned int correct_derivatives = m;

		//Coordinates of points
		T sumX = 0, sumY = 0;

		//Create matrix XT (1, 6) from X ( transposed )
		Matrix <T> XT = MatrixOperations::trans(X);

		//Get type of the direction
		const TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

		//Create temporary A matrix
		Matrix <T> J_T = J;

		//Get actual rotation
		const T R = X(0, 0);
		const Point3DGeographic <T> cart_pole(X(1, 0), X(2, 0));
		const T lat0 = X(3, 0);
		const T lon0 = X(4, 0);
		const T c = X(5, 0);
		const T alpha = X(6, 0);

		//Process all points: compute matrix of partial derivatives
		for (unsigned int i = 0; i < m; i++)
		{

			try
			{
				//Normal aspect: lat0, lon0
				if (aspect == NormalAspect)
				{
					//Upper part of the matrix: R, latp=90, lonp=0, lat0, lon0, alpha
					J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);


					//Lower part of the matrix: lat0, lon0: R, latp=90, lonp=0, lat0, lon0, alpha
					J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
				}

				//Transverse aspect: lonp, lat0
				else  if (aspect == TransverseAspect)
				{
					//Upper part of the matrix: R, latp=0, lonp, lat0, lon0=0, alpha
					J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);

					//Lower part of the matrix: R, latp=0, lonp, lat0, lon0=0, alpha
					J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
				}

				//Oblique aspect: latp, lonp, lat0
				else
				{
					//Upper part of the matrix: R, latp, lonp, lat0, lon0=0, alpha
					J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
					//double test = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>( pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), FirstDerivative, XT, VariableX4, NUM_DERIV_STEP, print_exceptions);
					J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					//J_T ( i, 4 ) = NumDerivative::getDerivative ( FProjEquationDerivative6Var <T> (  pl_reference [i]->getLat(), pl_reference [i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir ), XT, VariableX5, NUM_DERIV_STEP, print_exceptions );
					J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);

					//Lower part of the matrix: R, latp, lonp, lat0, lon0=0, alpha
					J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
					J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
					//J_T ( i + m, 4 ) = NumDerivative::getDerivative ( FProjEquationDerivative6Var <T> ( , pl_reference [i]->getLat(), pl_reference [i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir ), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions );
					J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
				}
			}

			//Decrease amount of corrected points
			catch (Exception & error)
			{
				correct_derivatives--;
			}
		}

		//J_T.print();

		//Compute column sums
		T sum0X = 0.0, sum1X = 0.0, sum2X = 0.0, sum3X = 0.0, sum4X = 0.0, sum5X = 0.0,
			sum0Y = 0.0, sum1Y = 0.0, sum2Y = 0.0, sum3Y = 0.0, sum4Y = 0.0, sum5Y = 0.0;

		for (unsigned int i = 0; i < m; i++)
		{
			//Sums of X derivatives
			sum0X += J_T(i, 0);
			sum1X += J_T(i, 1);
			sum2X += J_T(i, 2);
			sum3X += J_T(i, 3);
			sum4X += J_T(i, 4);
			sum5X += J_T(i, 5);

			//Sums of Y derivatives
			sum0Y += J_T(i + m, 0);
			sum1Y += J_T(i + m, 1);
			sum2Y += J_T(i + m, 2);
			sum3Y += J_T(i + m, 3);
			sum4Y += J_T(i + m, 4);
			sum5Y += J_T(i + m, 5);
		}

		//Compute Jacobi matrix
		for (unsigned int i = 0; i < m; i++)
		{

			//X derivatives
			J(i, 0) = (J_T(i, 0) - sum0X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 0) - sum0Y / m) * sin(alpha * M_PI / 180);
			J(i, 1) = (J_T(i, 1) - sum1X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 1) - sum1Y / m) * sin(alpha * M_PI / 180);
			J(i, 2) = (J_T(i, 2) - sum2X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 2) - sum2Y / m) * sin(alpha * M_PI / 180);
			J(i, 3) = (J_T(i, 3) - sum3X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 3) - sum3Y / m) * sin(alpha * M_PI / 180);
			J(i, 4) = (J_T(i, 4) - sum4X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 4) - sum4Y / m) * sin(alpha * M_PI / 180);
			J(i, 5) = (J_T(i, 5) - sum5X / m) * cos(alpha * M_PI / 180) - (J_T(i + m, 5) - sum5Y / m) * sin(alpha * M_PI / 180);
			J(i, 6) = (-(nl_projected[i]->getX() - x_mass_reference) * sin(alpha * M_PI / 180) - (nl_projected[i]->getY() - y_mass_reference) * cos(alpha * M_PI / 180)) * M_PI / 180; //Directly computed

			//Y derivatives
			J(i + m, 0) = (J_T(i, 0) - sum0X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 0) - sum0Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 1) = (J_T(i, 1) - sum1X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 1) - sum1Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 2) = (J_T(i, 2) - sum2X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 2) - sum2Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 3) = (J_T(i, 3) - sum3X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 3) - sum3Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 4) = (J_T(i, 4) - sum4X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 4) - sum4Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 5) = (J_T(i, 5) - sum5X / m) * sin(alpha * M_PI / 180) + (J_T(i + m, 5) - sum5Y / m) * cos(alpha * M_PI / 180);
			J(i + m, 6) = ((nl_projected[i]->getX() - x_mass_reference) * cos(alpha * M_PI / 180) - (nl_projected[i]->getY() - y_mass_reference) * sin(alpha * M_PI / 180)) *M_PI / 180; //Directly computed
		}

		//J.print();
		//Not enough points
		if (correct_derivatives < 3)
		{
			throw BadDataException("BadDataException: not enough correct partial derivatives, maybe error in equation. ", "Can not compute Jacobi matrix.");
		}

		//J.print();

		iter++;
	}
};

#endif
