// Description: Non-Linear Least Squares algorithms

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


#ifndef NonLinearLeastSquares_HPP
#define NonLinearLeastSquares_HPP

#include <cmath>

#include "libalgo/source/const2/Const.h"


//Set namespace
using namespace MatrixOperations;


template <typename T, typename FunctionJ, typename FunctionV>
T NonLinearLeastSquares::BFGSH(FunctionJ function_j, FunctionV function_v, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned int &iterations,
	const T alpha, const T nu, const T max_error, const unsigned int max_iterations, const T max_diff, std::ostream * output)
{
	//Solving Non-linear Least Squares using the hybrid BFGS algorithm
	//Combination of the Gauss-Newton and BFGS method, algorithm by L Luksan
	//Default values: nu = 0.0001, alpha = 0.0001;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), V2(m, 1), dX(n, 1), E(n, n, 0.0, 1.0);

	//Assign matrix
	Matrix <T> Y2 = Y;

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);
	//X.print();

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> H_new = H;
	Matrix <T> G_new = G;
	//H.print();
	//J.print();

	//Compute objective function
	T fx = norm(trans(V) * W * V), fx_new = fx;

	//Set iterations to 0
	iterations = 0;

	//Perform iterations
	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute new dX
		dX = pinv1(H) * G * (-1.0);

		//Too long step, reduction
		const T ndX = norm(dX);
		if (ndX > MAX_NLS_STEP_LENGTH) 
			dX = dX * 100 / ndX;

		//Compute new trial X
		Matrix <T> X2 = X + dX;

		//Reflection into the search space
		reflection(X2, A, B);

		//Compute new trial V matrix (residuals)
		function_v(X2, Y2, V2, W);

		//Apply back-step using a bisection
		const T t_min = 1.0e-10;
		T t = 1.0;
		
		while ((sum2(V2) > sum2(V) + sum(trans(V) * J *  dX) * t * alpha * 2.0) && (t > t_min))
		{
			//Step t bisection
			t /= 2;

			//Compute new X2
			X2 = X + dX * t;

			//Reflection into the search space
			reflection(X2, A, B);

			//Compute new V matrix: do not change parameters in one iteration step
			function_v(X2, Y2, V2, W);
		}
		
		//Compute new X using back-step method
		X = X + dX * t;

		//Reflection into the search space
		reflection(X, A, B);

		//Compute new V matrix and residual matrix V
		function_v(X, Y, V, W);

		//Compute new J matrix
		function_j(X, J);

		//std::cout << iterations;
		//X.print();

		//Compute new residuals and gradient
		G_new = trans(J) * W * V;
		fx_new = norm(trans(V) * W * V);
		
		//Terminal condition
		if ((norm(G) < max_error) /*|| (fabs(fx_new - fx) < 1.0 * max_diff * std::min(1.0, fx))*/ || (fx < max_error) || ( norm(dX) < 1.0e-10 ) )
		{
			break;
		}

		//Compute selection criterium
		const T df = (fx - fx_new) / fx;

		//Compute Hessian matrix as H=J*W*J (Gauss-Newton)
		if (df > nu)
		{
			H_new = trans(J) * W * J;
		}

		//Compute Hessian matrix from BFGS
		else
		{
			//Compute solution difference
			Matrix <T> s = dX * t;

			//Compute gradient difference
			Matrix <T> y = G_new - G;

			//Compute denominators
			const T ys = sum(trans(y) * s);
			const T shs = sum(trans(s) * H * s);

			//Compute update, if y * s > 0 (symmetric positive definite update)
			if (ys > 0)
			{
				Matrix <T> dH(n, n);

				if (fabs(ys) > MIN_FLOAT && fabs(shs) > MIN_FLOAT)
				{
					dH += y * trans(y) / ys - H * s * trans(H * s) / shs;
				}

				//Compute Hessian matrix using quasi-Newton update: H = J*W*J + dH
				H_new = trans(J) * W * J + dH;
			}

			//Do not update, if y * s < = 0
			else
				H_new = trans(J) * W * J;
		}

		//Assign values
		H = H_new;
		G = G_new;
		fx = fx_new;

		//X.print(); 
		//dX.print(output);
		//dX.print();
		//RES.print(output);*/
		//*output << dx << '\n';
		//*output << RES(0, 0) << '\n';
	}

	//Compute final values in V

	function_v(X, Y, V, W);

	//Evaluate minimum
	const T fxmin = norm(trans(V) * W * V);

	std::cout << " [" << iterations << " it., fmin = " << fxmin << "]" << '\n';

	return fxmin;
}



template <typename T>
void NonLinearLeastSquares::reflection(Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX)
{
	//Reflect elements of vectors into the search space represented by the n-dimensional cuboid
	const unsigned int n = X.rows();
	const T max_multiplier = 100;

	for (unsigned int i = 0; i < n; i++)
	{
		//int iter = 0;

		//Process each element of the vector
		while ((X(i, 0) < XMIN(i, 0)) || (X(i, 0) > XMAX(i, 0)))
		{
			//std::cout << X(i, 0) << " ";

			//XMIN == XMAX
			if (XMAX(i, 0) - XMIN(i, 0) < EPS)
			{
				X(i, 0) = XMIN(i, 0);
				break;
			}

			//Left form the lower bound
			else if (X(i, 0) > XMAX(i, 0))
			{

				X(i, 0) = 2 * XMAX(i, 0) - X(i, 0);
			}

			//Right to the upper bound
			else if (X(i, 0) < XMIN(i, 0))
			{

				 X(i, 0) = 2 * XMIN(i, 0) - X(i, 0);
			}
		}
	}
}


#endif
