// Description: Non Linear Least Squares algorithms

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


#ifndef NonLinearLeastSquares_HPP
#define NonLinearLeastSquares_HPP

#include <cmath>

#include "libalgo/source/exceptions/BadDataException.h"


//Set namespace
using namespace MatrixOperations;


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::GN ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations )
{
        //Solving non Linear Least Squares using the Gaussian-Newton method
        //Faster than backtracking method
        unsigned short iterations = 0;
        T cost_old = MAX_FLOAT;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V ( m, 1 ), dX ( n, 1 ), dX2 ( n, 1 );

        //Assign matrix
        Matrix <T> X_Old = X;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Perform iterations
        while ( iterations < max_iterations )
        {
		//Increment iterations
		iterations++;

                //Compute new J matrix: linearize solution using the Taylor approximation
                function_j ( X, J );

                //Stop computation
                if ( sum2 ( trans ( J ) * W * V ) < max_error )
                        break;

                //Compute Minimum Weighteed Least Squares using qr decomposition
                //Jacobian J = [ d_R, d_latp, d_lonp, d_lat0, d_lon0, d_dx, d_dy]
                //dX = mlsqr ( J, W, V ) * ( -1.0 );
                dX = pinv1 ( trans ( J ) * W * J ) * trans ( J ) * W * V * ( -1.0 );

                //Compute new X
                X = X + dX;

                //Compute new V matrix
                function_v ( X, Y, V, W );

                //X.print();
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::GND(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, const T alpha, const T max_error, const unsigned short max_iterations, const T max_diff, std::ostream * output)
{
        //Solving Non linear Least Squares using the damped Gaussian-Newton method
        T cost_old = MAX_FLOAT;

        //Set iterations to 0
        iterations = 0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V2 ( m, 1 ), dX ( n, 1 ), J_new ( m, n );

        //Assign matrix
        Matrix <T> Y2 = Y;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute matrices
        function_j ( X, J );
        Matrix <T> G = trans ( J ) * W * V;
        Matrix <T> F = trans ( V ) * W * V;
	//V.print();

        //Perform iterations
        while ( iterations < max_iterations )
        {
		//Increment iterations
		iterations++;

                //Compute new J matrix: linearize solution using the Taylor approximation
                //function_j ( X, J );

                //Compute the direction 
		dX = pinv1(trans(J) * W * J) * G * (-1.0);

		//dX.print();

                //Compute new trial X
                Matrix <T> X2 = X + dX;
		//Matrix <T> X2 = addStep(X, dX, A, B);

		//Reflection into the search space
		reflection(X2, A, B);

                //Compute new trial V matrix (residuals)
                function_v ( X2, Y2, V2, W );

		//V2.print();

		//Apply back-step using a bisection
		const T t_min = 1.0e-10;
		T t = 1.0;

		while ( (sum2(V2) > sum2(V) + sum(trans(V) * J *  dX) * t * alpha * 2.0) && (t > t_min) )
                {
                        //Step t bisection
                        t /= 2;

                        //Compute new X2
                        X2 = X + dX * t;
			//X2 = addStep(X, dX * t, A, B);

			//Reflection into the search space
			reflection(X2, A, B);

                        //Compute new V matrix: do not change parameters in one iteration step
                        function_v ( X2, Y2, V2, W, false );
                }

                //Compute new X using damping factor t
                X = X + dX * t;
		//X = addStep(X, dX * t, A, B);

		//Reflection into the search space
		reflection(X, A, B);

                //Compute V matrix
                function_v ( X, Y, V, W );

                //Compute new J matrix
                function_j ( X, J );

                //Compute new residuals and gradient
                Matrix <T> F_new = trans ( V ) * W * V;
                Matrix <T> G_new = trans ( J ) * W * V;
		
                //Terminal condition
                if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < max_diff * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                ( F ( 0, 0 ) <  max_error ) )
			break;

                //Assign old values
                G = G_new;
                F = F_new;

                //Increment iterations
                iterations++;

		//X.print();
		//J.print();
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //X.print();
        //X.print ( output );
        //T rr = norm ( trans (V) * V );
        //std::cout << "res =" << rr << '\n';

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}



template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::BFGS(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
	const T alpha, const T max_error, const unsigned short max_iterations, const T max_diff, std::ostream * output)
{
	//Solving non Linear Least Squares with using the BFGS algorithm
	//Default value: alpha = 0.0001;
	T cost_old = MAX_FLOAT;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), V2(m, 1), dX(n, 1), E(n, n, 0.0, 1.0);

	//Assign matrix
	Matrix <T> Y2 = Y;

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> H_new = H;
	Matrix <T> G_new = G;
	
	//Compute objective function
	T fx = norm(trans(V) * W * V), fx_new = fx, fx_dk = fx;

	//Set iterations to 0
	iterations = 0;

	//Perform iterations
	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute new dX
		dX = pinv1 ( H ) * G * ( -1.0 );
		//dX = mlsqr(J, W, V) * (-1.0);
		
		//dX.print();

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
			//X2 = addStep(X, dX * t, A, B);

			//Reflection into the search space
			reflection(X2, A, B);

			//Compute new V matrix: do not change parameters in one iteration step
			function_v(X2, Y2, V2, W, false);
		}


		//Compute new X using back-step method
		X = X + dX * t;

		//Reflection into the search space
		reflection(X, A, B);

		//Compute new V matrix and residual matrix V
		function_v(X, Y, V, W);

		//Compute new J matrix
		function_j(X, J);

		//Compute new residuals and gradient
		G_new = trans(J) * W * V;
		fx_new = norm(trans(V) * W * V);
		
		//Terminal condition
		if ((norm(G) < max_error) || (fabs(fx_new - fx) < 1.0 * max_diff * std::min(1.0, fx)) || (fx < max_error))
		{
			break;
		}

		//Compute solution difference
		Matrix <T> s = dX * t;

		//Compute gradient difference
		Matrix <T> y = G_new - G;

		//Compute denominators
		const T ys = norm(trans(y) * s);
		const T shs = norm(trans(s) * H * s);

		//Compute update, if y * s > 0 (symmetric positive definite update)
		if (ys > 0)
		{
			Matrix <T> dH(n, n);

			if ((ys > MIN_FLOAT) && (shs > MIN_FLOAT))
			{
				dH = y * trans(y) / ys - H * s * trans(H * s) / shs;
			}

			//Compute Hessian matrix using quasi-Newton update: H = J*W*J + dH
			H_new = trans(J) * W * J + dH;
		}

		//Do not update, if y * s < = 0
		else H_new = H;

		//Assign values
		H = H_new;
		G = G_new;
		//F = F_new;
	}

	//Compute final values in V
	function_v(X, Y, V, W);

	return norm(trans(V) * W * V);
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::BFGSH(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
	const T alpha, const T nu, const T max_error, const unsigned short max_iterations, const T max_diff, std::ostream * output)
{
        //Solving Non-linear Least Squares using the hybrid BFGS algorithm
        //Combination of the Gauss-Newton and BFGS method, algorithm by L Luksan
	//Default values: nu = 0.0001, alpha = 0.0001;

        //Create matrices
        const unsigned short k = 5, m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V2 ( m, 1 ),  dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Assign matrix
        Matrix <T> Y2 = Y;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );
	
        //Compute matrices
        Matrix <T> H = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;
	Matrix <T> H_new = H;
	Matrix <T> G_new = G;
	//H.print();
	//J.print();

	//Compute objective function
	T fx = norm(trans(V) * W * V), fx_new = fx;

        //Set iterations to 0
        iterations = 0;

        //Perform iterations
        while ( iterations < max_iterations )
        {
		//Increment iterations
		iterations++;

                //Compute new dX
		dX = pinv1(H) * G * (-1.0);

		//reflection2(X, dX, A, B);
		//H.print();
		//dX.print();

                //Compute new trial X
                Matrix <T> X2 = X + dX;
		
		//Reflection into the search space
		reflection(X2, A, B);

                //Compute new trial V matrix (residuals)
                function_v ( X2, Y2, V2, W, false );

                //Apply back-step using a bisection
                const T t_min = 1.0e-10; 
		T t = 1.0;

                while ( ( sum2 ( V2 ) > sum2 ( V ) + sum ( trans ( V ) * J *  dX ) * t * alpha * 2.0 ) && ( t > t_min ) )
                {
                        //Step t bisection
                        t /= 2;

			//Compute new X2
			X2 = X + dX * t;

			//Reflection into the search space
			reflection(X2, A, B);
			//reflection2(X, dX*t, A, B);

                        //Compute new V matrix: do not change parameters in one iteration step
                        function_v ( X2, Y2, V2, W, false );
                }

		//Compute new X using back-step method
                X = X + dX * t;

		//Reflection into the search space
		reflection(X, A, B);
		//reflection2(X, dX, A, B);

                //Compute new V matrix and residual matrix V
                function_v ( X, Y, V, W, true );

                //Compute new J matrix
                function_j ( X, J );

		//X.print();

                //Compute new residuals and gradient
		G_new = trans(J) * W * V;
		fx_new = norm( trans(V) * W * V );
		
		//if (iterations >= 100) break;

		//Terminal condition
		if ((norm(G) < max_error) || (fabs(fx_new - fx) < 1.0 * max_diff * std::min(1.0, fx)) || (fx < max_error))
		{
			break;
		}
		
		
                //Compute selection criterium
                const T df = ( fx - fx_new ) / fx;

                //Compute Hessian matrix as H=J*W*J (Gauss-Newton)
                if ( df > nu )
                {
                        H_new = trans ( J ) * W * J;
                }
		
                //Compute Hessian matrix from BFGS
                else
                {
			//Compute solution difference
			Matrix <T> s = dX * t;

			//Compute gradient difference
			Matrix <T> y = G_new - G;

			//Compute denominators
			const T ys = sum ( trans( y ) * s );
			const T shs = sum ( trans( s ) * H * s);

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


		//Test
		Matrix <T> X0(7, 1);

		//ECK5
		/*
		X0(0, 0) = 6378.0;
		X0(1, 0) = 90.0000000;
		X0(2, 0) = 0.0000000;
		X0(3, 0) = 50.0000000;
		X0(4, 0) = 0.0000085;
		X0(5, 0) = 1.0000000;
		X0(6, 0) = 90.0;
		*/

		//EQDC
		/*
		X0(0, 0) = 6378.0;
		X0(1, 0) = 90.0000000;
		X0(2, 0) = 0.0000000;
		X0(3, 0) = 50.0000000;
		X0(4, 0) = 24.263880643;
		X0(6, 0) = 108.587213657;
		
		//LAEA
		/*
		X0(0, 0) = 6378.0;
		X0(1, 0) = 52.0000000;
		X0(2, 0) = 10.0000000;
		X0(3, 0) = 0.0000000;
		X0(4, 0) = 0.0000085;
		X0(5, 0) = 1.0000000;
		X0(6, 0) = 90.0;
		*/
		/*
		//Merc
		X0(0, 0) = 6378.0;
		X0(1, 0) = 0.0000000;
		X0(2, 0) = -90.9120171;
		X0(3, 0) = 47.0000000;
		X0(4, 0) = 0;
		X0(5, 0) = 1.0000000;
		X0(6, 0) = -90.0;
		*/
		//Sinusoidal
		/*
		X0(0, 0) = 6378.0;
		X0(1, 0) = 90.0000000;
		X0(2, 0) = 0.0000000;
		X0(3, 0) = 50.0000000;
		X0(4, 0) = 0.0000085;
		X0(5, 0) = 1.0000000;
		X0(6, 0) = 90.0;
		*/
		//Werner-Staab
		/*
		X0(0, 0) = 6378.0;
		X0(1, 0) = 90.0000000;
		X0(2, 0) = 0.0000000;
		X0(3, 0) = 0.0000000;
		X0(4, 0) = 0.0000085;
		X0(5, 0) = 1.0000000;
		X0(6, 0) = 90.0;
		*/
		Matrix <T> RES = trans(V)*W*V;
		T dx = 0;
		if (X.rows() == 7)
		{
			dx = norm(X - X0) / norm(X0)*100.0;
		}
		else
		{
			//Matrix <T> XX = 
		}

		//X.print(output); 
		//dX.print(output);
		//dX.print();
		//RES.print(output);*/
		//*output << dx << '\n';
		//*output << RES(0, 0) << '\n';
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        std::cout << "iter:" << iterations << '\n';

        return norm ( trans ( V ) * W * V ) ;
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations,
	const T nu1, const T nu2, const T nu3, const T gamma1, const T gamma2, const T lambda_min, const T lambda_max, const T max_error, const unsigned short max_iterations, const T max_diff, T lambda, std::ostream * output)
{
	//Solving non Linear Least Squares with using the current Levenberq-Marquard algorithm
	//The original version with the Cauchy point initialization
	//Algorithm by Jorge Nocedal
	T cost_old = MAX_FLOAT;

	//Set iterations to 0
	iterations = 0;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), dX(n, 1), I(n, n, 0.0, 1.0);

	//Initialize matrices
	Matrix <T> Y2 = Y, V2 = V, W2 = W, J2 = J;

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> F = trans(V) * W * V * 0.5;
	Matrix <T> G2 = trans(G) * H * G;

	//Compute �nitial parameters
	T ng = norm(G);
	//lambda = 0.001 * ng;
	
	//lambda = std::max(norm(trans(G) * G) / norm(trans(G) * H  * G), 0.001);
	//T lambda =  0.01 *max ( diag ( H ) );
	Matrix <T> GB = trans(G) * H * G;
	//T lambda = std::max(1.0, ng * ng / fabs(GB(0, 0)));
	//lambda = 30 * ng;
	//lambda = 3000000;

	//Perform iterations
	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute Minimum Weighteed Least Squares using qr decomposition
		dX = pinv1(H + I * lambda) * G * (-1.0);

		//dX.print();

		//Compute new trial X
		//Matrix <T> X2 = X + dX;
		Matrix <T> X2 = addStep(X, dX, A, B);

		//Reflection into the search space
		reflection(X2, A, B);

		//Compute new V matrix and residual matrix V
		function_v(X2, Y2, V2, W2);

		//Compute new J matrix
		function_j(X2, J2);

		//Compute dQ: model function
		Matrix <T> dQM = trans(dX) * H * dX * 0.5 + trans(G) * dX;
		const T dQ = -dQM(0, 0);

		//Compute dF
		const Matrix <T> F2 = trans(V2) * W2 * V2 * 0.5;
		const T dF = F(0, 0) - F2(0, 0);
		const T rho = dF / dQ;

		//std::cout << "nu1= " << nu1 << "   nu2= " << nu2 << "   gam1 = " << gamma1 << "   gam2=" << gamma2 << "   lmin=" << lambda_min << "   l_max=" << lambda_max << '\n';
		//dX.print();
		//std::cout << "rho=" << rho;
		//std::cout << lambda << "  ";

		//Terminal condition
		if ((norm(G) < max_error) || (fabs(F2(0, 0) - F(0, 0)) < 1.0 * max_diff * std::max(1.0, F(0, 0))) ||
			(F(0, 0) <  max_error))
			break;

		//Worse agrrement between the model and functi, decrease lambda
		if (rho < nu1)
		{
			lambda = std::min(lambda * gamma1, lambda_max);
		}

		//Good agreement between the model and function, increase lambda
		else if (rho > nu2)
		{
			lambda = std::max(lambda * gamma2, lambda_min);
		}

		//Accept the step
		if (rho > nu3)
		{
			//Accept the trial solution
			X = X2;

			//Update solution
			F = F2;
			J = J2;
			V = V2;
			W = W2;

			//Actualize matrices
			H = trans(J) * W * J;
			G = trans(J) * W * V;
			ng = norm(G);
		}

		//X2.print();
		//X2.print(output);
	}

	//Compute final values in V
	function_v(X, Y, V, W);

	//X.print();
	//std::cout << "residuals";
	//V.print();
	//X.print(output);

	//Return squares of residuals
	return norm(trans(V) * W * V);
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::TR(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations,
	const T nu1, const T nu2, const T nu3, const T gamma1, const T gamma2, const T delta_min, const T delta_max, const T max_error, const unsigned short max_iterations, const T max_diff, T delta, std::ostream * output)
{
	//Solving non Linear Least Squares with using the trust region algorithm
	//Lambda solved using the More-Sorensen method
	//Algorithm by Conn
	T cost_old = MAX_FLOAT;

	//Set iterations to 0
	iterations = 0;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), dX(n, 1), I(n, n, 0.0, 1.0);

	//Initialize matrices
	Matrix <T> Y2 = Y, V2 = V, W2 = W, J2 = J;
	
	//Reflection into the search space
	reflection(X, A, B);

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> F = trans(V) * W * V * 0.5;
	Matrix <T> G2 = trans(G) * H * G;

	//Compute �nitial parameters
	T ng = norm(G);

	//Initialize trust region size
	//T delta = 0.00001 * ng;
	//T delta = 30 * ng;
	Matrix <T> GB = trans(G) * H * G;
	//T delta = std::max(15000.0, ng * ng / fabs(GB(0, 0)));
	//T delta = std::max(10.0, ng * ng / fabs(GB(0, 0)));
	//
	//delta = 3000000;
	//delta = 0.001 * ng;
	std::cout << ng << " ";
	//delta = 10000;
	//Perform iterations
	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute optimal step inside the trust region
		dX = optTRStep2(J, V, W, delta, 0.2 * max_iterations);

		//dX.print();

		//Compute new trial X
		Matrix <T> X2 = X + dX;
		//Matrix <T> X2 = addStep(X, dX, A, B);

		//Reflection into the search space
		reflection(X2, A, B);
		//Matrix <T> dX2 = X2 - X;
		//bool out = false;
		//if (norm(dX - dX2) > 0.1)
		//	out = true;

		//Compute new V matrix and residual matrix V
		function_v(X2, Y2, V2, W2);

		//Compute new J matrix
		function_j(X2, J2);

		//X2.print();

		//Compute new residuals
		const Matrix <T> F2 = trans(V2) * W2 * V2 * 0.5;

		//Terminal condition
		if ((norm(G) < max_error) || (fabs(F2(0, 0) - F(0, 0)) < 1.0 * max_diff * std::max(1.0, F(0, 0))) ||
			(F(0, 0) < max_error))
			break;

		//Compute dQ: model function
		Matrix <T> dQM = trans(dX) * H * dX * 0.5 + trans(G) * dX;
		const T dQ = -dQM(0, 0);

		//Compute dF
		const T dF = F(0, 0) - F2(0, 0);
		const T rho = dF / dQ;

		//std::cout << "F = " << F(0, 0) << "   FN = " << F2(0, 0) << "   dQ = " << dQ << '\n';
		//std::cout << "rho = " << rho << '\n';

		//Norm of the gradient to asses the new solution
		const T ndx = norm(dX);

		//Bad estimation of nu
		const T gdx = norm(trans(G) * dX);
		const T numer = (1 - 0.9) * gdx;
		const T denom = (1 - 0.9) * (F(0, 0) + gdx + 0.9 * (F(0, 0) + norm(dQM)) - F2(0, 0));
		const T nu_bad = numer / denom;

		//Accept the step
		if (rho >= 0.01 )
		{
			//Accept the trial solution
			X = X2;

			//Update solution
			F = F2;
			J = J2;
			V = V2;
			W = W2;

			//Actualize matrices
			H = trans(J) * W * J;
			G = trans(J) * W * V;
			ng = norm(G);
		}

		//std::cout << "F=" << F(0, 0) << "  F2=" << F2(0, 0) << "  Rho = " << rho << "  Delta = " << delta << '\n';
		
		/*
		//Terrible agreemenat of the model and function, rapidly shrink the trust region
		if (rho < 0 || out)
		{
			//std::cout << "minus   ";
			//std::cout << delta;
			delta = std::min(0.25 * ndx, std::max(0.0625, nu_bad)* delta);
			//std::cout << "   " << delta << '\n';
		}
		else
		*/
		//Worse agreement of the model and function shrink the trust region
		if (rho < 0.01 )
		{
			delta = 0.25 * ndx;
		}

		//Good agreement of the model and function, expand the trust region
		else if (rho > 0.9 )
		{
			delta = std::max(gamma2 * ndx, delta);
		}

		//X.print();
		//X2.print(output);
	}

	//Compute final values in V
	function_v(X, Y, V, W);

	//X.print();
	//std::cout << "cost = " << F(0, 0);
	//V.print();
	//X.print(output);

	//Return squares of residuals
	return norm(trans(V) * W * V);
}



template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::TR2(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
	const T nu1, const T nu2, const T nu3, const T gamma1, const T gamma2, const T delta_min, const T delta_max, const T max_error, const unsigned short max_iterations, const T max_diff, T delta, std::ostream * output)
{
	//Solving non Linear Least Squares with using the trust region algorithm
	//Lambda solved using the More-Sorensen method
	//Algorithm by Conn
	T cost_old = MAX_FLOAT;

	//Set iterations to 0
	iterations = 0;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), dX(n, 1), I(n, n, 0.0, 1.0);

	//Initialize matrices
	Matrix <T> Y2 = Y, V2 = V, W2 = W, J2 = J;

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> F = trans(V) * W * V * 0.5;
	Matrix <T> G2 = trans(G) * H * G;

	//Compute �nitial parameters
	T ng = norm(G);

	//Initialize trust region size
	//T delta = 30 * ng;
	Matrix <T> GB = trans(G) * H * G;
	//T delta = std::max(100.0, ng * ng / fabs(GB(0, 0)));
	delta = 10000;

	//Perform iterations
	unsigned short l = 0;
	T fmin = F(0, 0), fr = F(0, 0), fc = F(0, 0), sigma_r = 0, sigma_c = 0;

	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute optimal step inside the trust region
		dX = optTRStep2(J, V, W, delta, 0.2 * max_iterations);

		//		dX.print();

		//Compute new trial X
		//Matrix <T> X2 = X + dX;
		Matrix <T> X2 = addStep(X, dX, A, B);

		//Reflection into the search space
		//reflection(X2, A, B);

		//Compute new V matrix and residual matrix V
		function_v(X2, Y2, V2, W2);

		//Compute new J matrix
		function_j(X2, J2);

		//X2.print();

		//Compute new residuals
		const Matrix <T> F2 = trans(V2) * W2 * V2 * 0.5;

		//Terminal condition
		if ((norm(G) < max_error) || (fabs(F2(0, 0) - F(0, 0)) < 1.0 * max_diff * std::max(1.0, F(0, 0))) ||
			(F(0, 0) < max_error))
			break;


		//Compute dF
		const T dFh = fr - F2(0, 0);
		const T dFc = F(0,0) - F2(0, 0);

		//Compute dQ: model function
		Matrix <T> dQM = trans(dX) * H * dX * 0.5 + trans(G) * dX;
		const T dQ = -dQM(0, 0);

		//Compute Rho
		const T rho_h = dFh / ( dQ + sigma_r);
		const T rho_c = dFc / dQ;
		const T rho = std::max(rho_h, rho_c);

		//std::cout << "F = " << F(0, 0) << "   FN = " << F2(0, 0) << "   dQ = " << dQ << '\n';
		//std::cout << "rho = " << rho << '\n';

		//Accept the step
		if (rho >= 0.01)
		{
			//Accept the trial solution
			X = X2;

			//Update sigma
			sigma_c += dQ;
			sigma_r += dQ;

			//Update solution
			F = F2;
			J = J2;
			V = V2;
			W = W2;

			//Actualize matrices
			H = trans(J) * W * J;
			G = trans(J) * W * V;
			ng = norm(G);

			//Update the best value
			if (F2(0, 0) < fmin)
			{
				fc = F2(0, 0);
				fmin = F2(0, 0);
				sigma_c = 0;
				l = 0;
			}

			else
			{
				l++;

				//Update the candidate
				if (F2(0, 0) > fc)
				{
					fc = F2(0, 0);
					sigma_c = 0;
				}

				//Reset the reference value
				if (l == 45)
				{
					fr = fc;
					sigma_r = sigma_c;
				}
			}
		}

		//Norm of the gradient to asses the new solution
		const T ndx = norm(dX);

		//Bad estimation of nu
		const T gdx = norm(trans(G) * dX);
		const T numer = (1 - 0.9) * gdx;
		const T denom = (1 - 0.9) * (F(0, 0) + gdx + 0.9 * (F(0, 0) + norm(dQM)) - F2(0, 0));
		const T nu_bad = numer / denom;

		//std::cout << "F=" << F(0, 0) << "  F2=" << F2(0, 0) << "  Rho = " << rho << "  Delta = " << delta << '\n';
		/*
		//Terrible agreemenat of the model and function, rapidly shrink the trust region
		if (rho < 0)
		{
		//std::cout << "minus   ";
		//std::cout << delta;
		delta = std::min(0.25 * ndx, std::max(0.0625, nu_bad)* delta);
		//std::cout << "   " << delta << '\n';
		}

		//Worse agreement of the model and function shrink the trust region
		else */
		if (rho < 0.01)
		{
			delta = 0.25 * ndx;
		}

		//Good agreement of the model and function, expand the trust region
		else if (rho > 0.9)
		{
			delta = std::max(gamma2 * ndx, delta);
		}

		//X.print();
		//X2.print(output);

	}

	//Compute final values in V
	function_v(X, Y, V, W);

	//X.print();
	//std::cout << "cost = " << F(0, 0);
	//V.print();
	//X.print(output);

	//Return squares of residuals
	return norm(trans(V) * W * V);
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::TR3(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations,
	const T nu1, const T nu2, const T nu3, const T gamma1, const T gamma2, const T delta_min, const T delta_max, const T max_error, const unsigned short max_iterations, const T max_diff, T delta, std::ostream * output)
{
	//Solving non Linear Least Squares with using the trust region algorithm
	//Lambda solved using the More-Sorensen method
	//Algorithm by Conn
	T cost_old = MAX_FLOAT;

	//Set iterations to 0
	iterations = 0;

	//Create matrices
	const unsigned short m = W.rows(), n = X.rows();
	Matrix <T> J(m, n), dX(n, 1), I(n, n, 0.0, 1.0);

	//Initialize matrices
	Matrix <T> Y2 = Y, V2 = V, W2 = W, J2 = J;

	//Compute initial V matrix (residuals)
	function_v(X, Y, V, W);

	//Compute initial J matrix
	function_j(X, J);

	//Compute matrices
	Matrix <T> H = trans(J) * W * J;
	Matrix <T> G = trans(J) * W * V;
	Matrix <T> F = trans(V) * W * V * 0.5;
	Matrix <T> G2 = trans(G) * H * G;

	//Compute �nitial parameters
	T ng = norm(G);

	//Initialize trust region size
	delta = 3000 * ng;
	Matrix <T> GB = trans(G) * H * G;
	//delta = std::max(10.0, ng * ng / fabs(GB(0, 0)));
	//std::cout << "ng = " << ng;
	//delta = 300000000;

	//Perform iterations
	T rho_old = 0;
	unsigned short j = 0;
	while (iterations < max_iterations)
	{
		//Increment iterations
		iterations++;

		//Compute optimal step inside the trust region
		dX = optTRStep2(J, V, W, delta, 0.2 * max_iterations);

				//dX.print()

		//Compute new trial X
		//Matrix <T> X2 = X + dX;
		Matrix <T> X2 = addStep(X, dX, A, B);

		//Reflection into the search space
		//reflection(X2, A, B);

		//Compute new V matrix and residual matrix V
		function_v(X2, Y2, V2, W2);

		//Compute new J matrix
		function_j(X2, J2);

		//X2.print();

		//Compute new residuals
		Matrix <T> F2 = trans(V2) * W2 * V2 * 0.5;

		//Terminal condition
		if ((norm(G) < max_error) || (fabs(F2(0, 0) - F(0, 0)) < 1.0 * max_diff * std::max(1.0, F(0, 0))) ||
			(F(0, 0) < max_error))
			break;

		//Compute dQ: model function
		Matrix <T> dQM = trans(dX) * H * dX * 0.5 + trans(G) * dX;
		const T dQ = -dQM(0, 0);

		//Compute dF
		const T dF = F(0, 0) - F2(0, 0);
		T rho = dF / dQ;

		//std::cout << "F = " << F(0, 0) << "   FN = " << F2(0, 0) << "   dQ = " << dQ << '\n';
		//std::cout << "rho = " << rho << '\n';

		if ((rho >= 0.9) && (rho_old >= 0.9))
		{
			bool cont_iter = true;

			//Assign step and trust region size
			Matrix <T> dXj1 = dX;
			T deltaj1 = delta;
			T rho_max = 0;

			while (cont_iter)
			{
				Matrix <T> Y2j = Y, V2j = V, W2j = W, J2j = J, F2j = F2;

				//Actualize size of the trust region
				const T deltaj = 2.0 * deltaj1;

				//Compute new trust region step
				Matrix <T> dXj = optTRStep2(J, V, W, deltaj, 0.2 * max_iterations);

				//Create new trial point
				//Matrix <T> X2j = X + dXj;
				Matrix <T> X2j = addStep(X, dXj, A, B);

				//Reflection into the search space
				//reflection(X2j, A, B);

				//X2j.print();

				//Compute new V matrix and residual matrix V
				function_v(X2j, Y2j, V2j, W2j);

				//Compute new J matrix
				function_j(X2j, J2j);

				//X2.print();

				//Compute new residuals
				F2j = trans(V2j) * W2j * V2j * 0.5;

				//Compute new rho
				Matrix <T> dQMj = trans(dXj) * H * dXj * 0.5 + trans(G) * dXj;
				const T dQj = -dQMj(0, 0);

				//Compute dF
				const T dFj = F(0, 0) - F2j(0, 0);
				const T rhoj = dFj / dQj;

				//Accept the step, update values
				//if (rhoj >= rho_max)
				if (rhoj >= 0.01)
				{
					rho_max = rhoj;

					//Assign determined parameters
					dX = dXj;
					rho = rhoj;
					delta = deltaj;
					X2 = X2j;
					F2 = F2j;
					J2 = J2j;
					V2 = V2j;
					W2 = W2j;
				}

				//Continue with iterations
				cont_iter = (rhoj >= 0.9) && (norm(dXj) > norm(dXj1));

				//Assign old values
				dXj1 = dXj;
				deltaj1 = deltaj;
			} 
		}


		//Norm of the gradient to asses the new solution
		const T ndx = norm(dX);

		//Bad estimation of nu
		const T gdx = norm(trans(G) * dX);
		const T numer = (1 - 0.9) * gdx;
		const T denom = (1 - 0.9) * (F(0, 0) + gdx + 0.9 * (F(0, 0) + norm(dQM)) - F2(0, 0));
		const T nu_bad = numer / denom;

		//Accept the step
		if (rho >= 0.01)
		{
			//Accept the trial solution
			X = X2;

			//Update solution
			F = F2;
			J = J2;
			V = V2;
			W = W2;

			//Actualize matrices
			H = trans(J) * W * J;
			G = trans(J) * W * V;
			ng = norm(G);
		}

		//std::cout << "F=" << F(0, 0) << "  F2=" << F2(0, 0) << "  Rho = " << rho << "  Delta = " << delta << '\n';
		/*
		//Terrible agreemenat of the model and function, rapidly shrink the trust region
		if (rho < 0)
		{
			//std::cout << "minus   ";
			//std::cout << delta;
			delta = std::min(0.25 * ndx, std::max(0.0625, nu_bad)* delta);
			//std::cout << "   " << delta << '\n';
		}

		//Worse agreement of the model and function shrink the trust region
		else */
			if (rho < 0.01)
		{
			delta = 0.25 * ndx;
		}

		//Good agreement of the model and function, expand the trust region
		else if (rho > 0.9)
		{
			delta = std::max(gamma2 * ndx, delta);
		}

		//Assign rho
		rho_old = rho;
		//X.print();
		//X2.print(output);
	}

	//Compute final values in V
	function_v(X, Y, V, W);

	//X.print();
	//std::cout << "cost = " << F(0, 0);
	//V.print();
	//X.print(output);

	//Return squares of residuals
	return norm(trans(V) * W * V);
}



template <typename T>
Matrix <T> NonLinearLeastSquares::optTRStep(const Matrix <T> &J, const Matrix <T> &V, const Matrix <T> &W, const T delta, const T max_iter, const T max_error)
{
	//More-Sorensen method: determine optimal trust region step
	//Algorithm by L. Luksan
	unsigned short iterations = 0;
	const unsigned short n = J.cols();
	bool stop = false;

	//Determined step inside the trust region
	Matrix <T> s(n, 1), I(n, n, 0.0, 1.0), L(n, n), d(n, 1), e(n, 1), v(n, 1);

	//Parameters of the algorithm
	const T beta_min = 0.001, delta_min = 0.9, delta_max = 1.1;

	//Compute initial matrices
	const Matrix <T> B = trans(J) * W * J;
	const Matrix <T> G = trans(J) * W * V;

	//Eigenvalues decomposition
	Matrix <T> EV(n, n, 0, 1), EL(n, 1);
	eig(B, EV, EL);

	//Get the smalles and largest eigenvalues
	const T l_min = EL(0, 0);
	const T l_max = EL(n - 1, 0);
	//std::cout << "lmax" << l_max;

	//Initialize bounds for lambda
	T mju_min = max(diag(B * (-1.0)));
	T lambda_min = norm(G) / delta - l_max;
	T lambda_max = norm(G) / delta - l_min;

	//Initialize lambda
	T lambda = std::max(0.0, std::max(mju_min, lambda_min));

	//Perform iteration
	while (iterations < max_iter)
	{
		//***Trust region, Step 2***

		//Set up the threshold
		lambda_min = std::max(0.0, std::max(mju_min, lambda_min));
		if ((lambda <= mju_min) && (iterations > 0))
			lambda = std::max(sqrt(lambda_min * lambda_max), lambda_min + beta_min * (lambda_max - lambda_min));

		//std::cout << "iter = " << iterations << "  lam_min = " << lambda_min << "  lam_max = " << lambda_max << "   lam = " << lambda << '\n';
		
		//Increment iterations
		iterations++;
		//***Trust region, Step 3***

		//Gill-Murray decomposition
		Matrix <double> R(n, n), E(n, n);
		bool indefinite = false;
		gill(B + I * lambda, R, E, indefinite);
		gill(B + I * lambda, L, d, e, v);

		//Find maximum diagonal element of E
		const T max_diag_E = max(abs(diag(e)));

		//If (B + lambda * I) is positive definite diag|E| = 0
		if (max_diag_E < 1.0e-3)
		{
			trStep4(B, R, s, G, v, lambda_min, lambda_max, lambda, delta_min, delta_max, delta, mju_min, stop);
		}

		// (B + lambda * I) is not positive definite: diag|E| != 0
		else
		{
			//Normalize direction
			v = v * (1.0 / norm(v));

			//Actualize mju_min
			const Matrix <T> vb = trans(v) * (B + I * lambda) * v;

			//mju_min = lambda - vb(0, 0);
			mju_min = std::max(mju_min, lambda - vb(0, 0));
		}

		//std::cout << "iter = " << iterations << "   lam = " << lambda <<  "  delta = " << delta << '\n';

		//Return results
		if ( stop )
			return s;
	}

	return s;
}



template <typename T>
Matrix <T> NonLinearLeastSquares::optTRStep2(const Matrix <T> &J, const Matrix <T> &V, const Matrix <T> &W, const T delta, const T max_iter, const T max_error)
{
	//More-Sorensen method: determine optimal trust region step
	//Algorithm by L. Luksan
	unsigned short iterations = 0;
	const unsigned short n = J.cols();
	bool stop = false;

	//Determined step inside the trust region
	Matrix <T> s(n, 1), I(n, n, 0.0, 1.0), L(n, n), d(n, 1), e(n, 1), v(n, 1);

	//Parameters of the algorithm
	const T beta_min = 0.001, beta_max = 0.999, delta_min = 0.9, delta_max = 1.1;

	//Compute initial matrices
	const Matrix <T> B = trans(J) * W * J;
	const Matrix <T> G = trans(J) * W * V;

	//Eigenvalues decomposition
	Matrix < T> EV(n, n, 0, 1), EL(n, 1);
	eig(B, EV, EL);

	//Get the smalles and largest eigenvalues
	const T l_min = EL(0, 0);
	const T l_max = EL(n - 1, 0);
	//std::cout << "lmax" << l_max;

	Matrix <T> SB = sumCols(B);
	const T max_el = max(SB);

	//Initialize bounds for lambda
	T mju_min = max(diag(B * (-1.0)));
	T lambda_min = norm(G) / delta - max_el;
	T lambda_max = norm(G) / delta + max_el;

	//Initialize lambda
	T lambda = std::max(0.0, std::max(mju_min, lambda_min));

	//Perform iteration
	while (iterations < max_iter)
	{
		//***Trust region, Step 2***
		Matrix <T> R(n, n);

		//Set up the threshold
		lambda_min = std::max(0.0, std::max(mju_min, lambda_min));

		//Go to step 3
		if ( lambda > mju_min /*|| iterations == 0*/)
		{
			trStep32(B, R, s, G, v, lambda_min, lambda_max, lambda, delta_min, delta_max, delta, mju_min, stop);
		}

		else if (lambda <= mju_min )
		{
			const T lambda_plus = sqrt(lambda_min * lambda_max);

			//Upper and lower bound
			const T lambda_low = lambda_min + beta_min * (lambda_max - lambda_min);
			const T lambda_high = lambda_min + beta_max * (lambda_max - lambda_min);
			
			//Set new value of lambda
			if (lambda_plus < lambda_low)
			{
				lambda = lambda_low;
			}

			else if (lambda_plus > lambda_high)
			{
				lambda = lambda_high;
			}

			else
			{
				lambda = lambda_plus;
			}
		}

		//Go to step 3
		//trStep32(B, R, s, G, v, lambda_min, lambda_max, lambda, delta_min, delta_max, delta, mju_min, stop);
			
		//std::cout << "iter = " << iterations << "  lam_min = " << lambda_min << "  lam_max = " << lambda_max << "   lam = " << lambda << '\n';

		//Increment iterations
		iterations++;
	
		//Return results
		if (stop)
			return s;
	}

	return s;
}



template <typename T>
void NonLinearLeastSquares::trStep32(const Matrix <T> &B, Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min, bool &stop )		    
{
	//Gill-Murray decomposition
	const unsigned int n = B.rows();
	Matrix <T> E(n, n), I(n, n, 0.0, 1.0), L(n, n), d(n, 1), e(n, 1);

	bool indefinite = false;
	gill(B + I * lambda, R, E, indefinite);
	gill(B + I * lambda, L, d, e, v);

	//Find maximum diagonal element of E
	const T max_diag_E = max(abs(diag(e)));

	//If (B + lambda * I) is positive definite diag|E| = 0
      	if (max_diag_E < 1.0e-3)
	{
		trStep42(B, R, s, G, v, lambda_min, lambda_max, lambda, delta_min, delta_max, delta, mju_min, stop);
	}

	// (B + lambda * I) is not positive definite: diag|E| != 0
	else
	{
		//Normalize direction
		v = v * (1.0 / norm(v));

		//Actualize mju_min
		const Matrix <T> vb = trans(v) * (B + I * lambda) * v;
		
		mju_min = lambda - vb(0, 0);
		//mju_min = std::max(mju_min, lambda - vb(0, 0));
	}
}



template <typename T>
void NonLinearLeastSquares::trStep4(const Matrix <T> &B, const Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min,  bool &stop)
{
	//Step 4 of the trust region subproblem

	//Compute S, so as R'RS + g = 0
	//Matrix<T> kkk =B + I *lambda;
	//kkk.print();
	//B.print();
	s = inv(trans(R) * R) * G * (-1.0);
	//R.print();
	//S.print();
	//S = pinv1(B) * G * (-1.0);
	//dX.print();

	//Actualize lambda, go to Step 6
	const T ns = norm(s);
	if (ns > delta_max * delta)
	{
		//Actualize lambda
		lambda_min = lambda;

		//Go to step 6
		lambda += trStep6(R, s, v, delta);
	}

	//Stop computation
	else if ((delta_min * delta <= ns) && (ns <= delta_max * delta) || (ns < delta_min * delta) && (fabs(lambda) < 0.0001))
	{
		stop = true;
	}

	//Eigenvalues  decomposition, go to Step 5
	else if ((ns < delta_min * delta) && (lambda > 0))
	{
		// Actualize threshold
		lambda_max = lambda;

		//Goto step 5
		trStep5(B, R, s, v, lambda, delta_min, delta, mju_min, stop);
	}
}


template <typename T>
void NonLinearLeastSquares::trStep42(const Matrix <T> &B, const Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min, bool &stop)
{
	//Step 4 of the trust region subproblem

	//Compute S, so as R'RS + g = 0
	//Matrix<T> kkk =B + I *lambda;
	//kkk.print();
	//B.print();
	s = pinv1(trans(R) * R) * G * (-1.0);
	//R.print();
	//S.print();
	//S = pinv1(B) * G * (-1.0);
	//dX.print();

	//Actualize lambda, go to Step 6
	const T ns = norm(s);
	if (ns > delta_max * delta)
	{
		//Actualize lambda
		lambda_min = lambda;

		//Go to step 6
		trStep62(R, s, G, v, lambda_min, lambda_max, lambda, delta, mju_min);
	}

	//Stop computation
	else if ((delta_min * delta <= ns) && (ns <= delta_max * delta) || (ns < delta_min * delta) && (fabs(lambda) < MIN_FLOAT))
	{
		stop = true;
	}

	//Eigenvalues  decomposition, go to Step 5
	else if ((ns < delta_min * delta) && (lambda > 0))
	{
		// Actualize threshold
		lambda_max = lambda;

		//Goto step 5
		trStep52(B, R, s, G, v, lambda_min, lambda_max, lambda, delta_min, delta, mju_min, stop);
	}
}




template <typename T>
void NonLinearLeastSquares::trStep5(const Matrix <T> &B, const Matrix <T> &R, Matrix<T> &s, Matrix<T> &v, T &lambda, const T delta_min, const T delta, T &mju_min, bool &stop)
{
	// Step 5 of the trust region subproblem
	const unsigned short n = B.cols();

	Matrix <T> EV(n, n, 0, 1), EL(n, 1), I(n, n, 0.0, 1.0);

	//Eigenvalue decomposition
	//eig(B + I * lambda, EV, EL);
	//eig(B, EV, EL);
	eig(B, EV, EL);

	//Get the smallest eigenvector and norm
	v = trans(EV(0, 0, 0, n - 1));
	v = v * (1.0 / norm(v));

	//Compute alpha for vn'*S > 0
	const T s2 = norm(s)  * norm(s);
	const T vts = norm(trans(v) * s);
	std::cout << "vts = " << vts;

	const T vtsd_sqr = sqrt(vts * vts + delta* delta - s2);
	const T alpha = (vts <= 0 ? vtsd_sqr - vts : (delta * delta - s2) / (vtsd_sqr + vts));

	//Test condition
	const T rvn2 = norm(R * v) * norm(R * v);
	const T rs2 = norm(R * s) * norm(R * s);
	const T l = alpha * alpha * rvn2;
	const T r = (1 - delta_min * delta_min) * (rs2 + lambda * delta * delta);

	//Terminate step
	if (l <= r)
	{
		//Actualize S
		s = s + v * alpha;

		//Terminate step
		stop = true;
	}

	//Update mju_min and lambda
	else

	{
		//mju_min = lambda - rvn2;
		mju_min = std::max(mju_min, lambda - rvn2);

		//Goto Step 6, compute new lambda
		lambda += trStep6(R, s, v, delta);
	}
}


template <typename T>
void NonLinearLeastSquares::trStep52(const Matrix <T> &B, const Matrix <T> &R, Matrix<T> &s, const Matrix <T> &G, Matrix<T> &v, const T lambda_min, const T lambda_max, T &lambda, const T delta_min, const T delta, T &mju_min, bool &stop)
{
	// Step 5 of the trust region subproblem
	const unsigned short n = B.cols();

	Matrix <T> EV(n, n, 0, 1), EL(n, 1), I(n, n, 0.0, 1.0);

	//Eigenvalue decomposition
	//eig(B + I * lambda, EV, EL);
	//eig(B, EV, EL);
	eig(B, EV, EL);

	//Get the smallest eigenvector and norm
	v = trans(EV(0, 0, 0, n - 1));
	v = v * (1.0 / norm(v));

	//Compute alpha for vn'*S > 0
	const T s2 = norm(s)  * norm(s);
	const T vts = norm(trans(v) * s);
	std::cout << "vts = " << vts;

	const T vtsd_sqr = sqrt(vts * vts + delta* delta - s2);
	const T alpha = (delta * delta - s2) / (vtsd_sqr + vts);

	//Test condition
	const T rvn2 = norm(R * v) * norm(R * v);
	const T rs2 = norm(R * s) * norm(R * s);
	const T l = alpha * alpha * rvn2;
	const T r = (1 - delta_min * delta_min) * (rs2 + lambda * delta * delta);

	//Terminate step
	if (l <= r)
	{
		//Actualize S
		s = s + v * alpha;

		//Terminate step
		stop = true;
	}

	//Update mju_min and lambda
	else
	{
		mju_min = lambda - rvn2;
		//mju_min = std::max(mju_min, lambda - rvn2);

		//Goto Step 6, compute new lambda
		trStep62(R, s, G, v, lambda_min, lambda_max, lambda, delta, mju_min);
	}
}


template <typename T>
T NonLinearLeastSquares::trStep6(const Matrix <T> &R, const Matrix <T> &s, Matrix<T> &v, const T delta)
{

	// Step 6 of the trust region subproblem
	v = inv(trans(R)) * s;

	//Actualize lambda solving the Newton iteration step
	const T ns = norm(s);
	const T nv = norm(v);
	const T s2 = ns * ns;
	const T v2 = nv * nv;
	
	const T d_lambda = s2 / v2 * (norm(s) - delta) / delta;

	return d_lambda;
}


template <typename T>
void NonLinearLeastSquares::trStep62(const Matrix <T> &R, const Matrix <T> &s, const Matrix <T> &G, Matrix<T> &v, const T lambda_min, const T lambda_max, T &lambda, const T delta, const T mju_min )
{

	// Step 6 of the trust region subproblem
	if (norm(G) < 0.000001)
	{
		lambda = mju_min;
		return;
	}
	
	//Compute step
	v = pinv1(trans(R)) * s;

	//Actualize lambda solving the Newton iteration step
	const T ns = norm(s);
	const T nv = norm(v);
	const T s2 = ns * ns;
	const T v2 = nv * nv;

	const T lambda_plus = lambda + s2 / v2 * (ns - delta) / delta;

	//Actualize lambda
	if (lambda_plus < lambda_min)
	{
		lambda = lambda_min;
	}

	else if (lambda_plus > lambda_max)
	{
		lambda = lambda_max;
	}

	else
	{
		lambda = lambda_plus;
	}
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
			if (XMAX(i, 0) - XMIN(i, 0) < MAX_FLOAT_OPER_ERROR)
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


template <typename T>
void NonLinearLeastSquares::reflection2(Matrix <T> &X, Matrix <T> &dX, const Matrix <T> &XMIN, const Matrix <T> &XMAX)
{
	//Reflect elements of vectors into the search space represented by the n-dimensional cuboid
	const unsigned int n = X.rows();
	T beta = 1.0;
	
	for (int i = 0; i < n; i++)
	{
		//Denominator
		const T den = std::max(fabs(dX(i, 0)), 0.0000001);

		if (X(i, 0) + dX(i, 0) > XMAX(i, 0))
		{
			const T beta_max = fabs((XMAX(i, 0) - X(i, 0)) / den);
			beta = std::min(beta, min(beta_max, 1.0));
			std::cout << "beta_max = " << beta_max  << '\n';
		}

		else if (X(i, 0) + dX(i, 0) < XMIN(i, 0))
		{
			const T beta_min = fabs((XMIN(i, 0) - X(i, 0)) / den);
			beta = std::min(beta, min(beta_min, 1.0));
			std::cout << "beta_min = " << beta_min << '\n';
		}

		//std::cout << "beta_min = " << beta_min << "   "  << "beta_max" <<beta_max << '\n';

		//beta = std::min(beta, min(min(beta_min, beta_max), 1.0));
	}

	dX = dX * beta;
	std::cout << "beta = " << beta << '\n';
}



template <typename T>
void NonLinearLeastSquares::reflection3(Matrix <T> &X, Matrix <T> &dX, const Matrix <T> &XMIN, const Matrix <T> &XMAX)
{
	//Reflect elements of vectors into the search space represented by the n-dimensional cuboid
	const unsigned int n = X.rows();
	Matrix <T> beta(n, 1, 1);

	for (int i = 0; i < n; i++)
	{
		//Denominator
		const T den = std::max(fabs(dX(i, 0)), 0.0000001);

		if (X(i, 0) + dX(i, 0) > XMAX(i, 0))
		{
			const T beta_max = fabs((XMAX(i, 0) - X(i, 0)) / den);
			beta (i, 0) =   0.25 * min(beta_max, 1.0);
			std::cout << "beta_max = " << beta_max << '\n';
		}

		else if (X(i, 0) + dX(i, 0) < XMIN(i, 0))
		{
			const T beta_min = fabs((XMIN(i, 0) - X(i, 0)) / den);
			beta(i, 0) = 0.25 * min(beta_min, 1.0);
			std::cout << "beta_min = " << beta_min << '\n';
		}

		//std::cout << "beta_min = " << beta_min << "   "  << "beta_max" <<beta_max << '\n';

		//beta = std::min(beta, min(min(beta_min, beta_max), 1.0));
	}

	dX = dX % beta;
	//std::cout << "beta = " << beta << '\n';
}

/*
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
			if (XMAX(i, 0) - XMIN(i, 0) < MAX_FLOAT_OPER_ERROR)
			{
				X(i, 0) = XMIN(i, 0);
				break;
			}

			//Left form the lower bound
			else if (X(i, 0) > XMAX(i, 0))
			{
				
				//X(i) too far from the upper bound (due to the pseudoinverse)
				if (X(i, 0) - XMAX(i, 0) > max_multiplier * (XMAX(i, 0) - XMIN(i, 0)))
				{
					X(i, 0) = 2 * XMAX(i, 0) - fmod(X(i, 0), XMAX(i, 0) - XMIN(i, 0));
				}
				
				//X(i) close to the upper bound
				else
					
				{
					X(i, 0) = 2 * XMAX(i, 0) - X(i, 0);
				}
			}

			//Right to the upper bound
			else if (X(i, 0) < XMIN(i, 0))
			{
				
				//X(i) too far from the lower bound (due to the pseudoinverse)
				if (XMIN(i, 0) - X(i, 0) > max_multiplier * (XMAX(i, 0) - XMIN(i, 0)))
				{
					X(i, 0) = 2 * XMIN(i, 0) - fmod(X(i, 0), XMAX(i, 0) - XMIN(i, 0));
				}
				
				//X(i) close to the lower bound
				else
				
				{
					X(i, 0) = 2 * XMIN(i, 0) - X(i, 0);
				}
			}
		}
	}
}

*/

/*
template <typename T>
void NonLinearLeastSquares::reflection(Matrix <T> &X, const Matrix <T> &A, const Matrix <T> &B)
{
	//Reflection into the search space
	const unsigned int n = X.rows();
	const T eps = 1.0e-5;

	//Reflect elements of vectors into the search space represented by the n-dimensional cuboid
	for (unsigned int i = 0; i < n; i++)
	{
		if ((A(i, 0) < B(i, 0)) && ( X(i, 0) < A(i, 0)))
		{
			double x = X(i, 0);
			double rem = fmod(x - A(i, 0), (B(i, 0) - A(i, 0))) - A(i, 0);
			if (i > 0 && rem < -180)
				std::cout << "error" << x << "  " << rem << "  " << A(i, 0) << "  " << B(i, 0) << '\n';

			X(i, 0) = fmod(X(i, 0) - A(i, 0), (B(i, 0) - A(i, 0))) - A(i, 0);

		}

		else if ( (A(i, 0) < B(i, 0)) && (X(i, 0) > B(i, 0)))
		{
			double x = X(i, 0);
			double rem = fmod(x - B(i, 0), (B(i, 0) - A(i, 0))) - B(i, 0);
			if (i > 0 && rem > 180)
				std::cout << "error" << x << "  " << rem << "  " << A(i, 0) << "  " << B(i, 0) << '\n';

			X(i, 0) = fmod(X(i, 0) - B(i, 0), (B(i, 0) - A(i, 0))) - B(i, 0);

		}
	}
}
*/

template <typename T>
bool NonLinearLeastSquares::checkIntervals(Matrix <T> &X, const Matrix <T> &A, const Matrix <T> &B)
{
	//Reflection into the search space
	const unsigned int n = X.rows();

	//Reflect elements of vectors into the search space represented by the n-dimensional cuboid
	for (unsigned int i = 0; i < n; i++)
	{
		//Check, if item is inside the search space
		if ((A(i, 0) < B(i, 0)) && (X(i, 0) < A(i, 0) || X(i, 0) > B(i, 0)))
		{
			return false;
		}
	}

	return true;
}


template <typename T>
Matrix <T> NonLinearLeastSquares::addStep(const Matrix <T> &X, const Matrix <T> &dX, const Matrix <T> &A, const Matrix <T> &B)
{
	//Add steps, elements ouside the interval are replaced with the maximum value
	unsigned int n = X.rows();
	
	Matrix <T> X2 = X;

	for ( unsigned int i = 0; i < n; i++)
	{
		//Inside the interval
		if ((X(i, 0) + dX(i, 0) >= A(i, 0)) && (X(i, 0) + dX(i, 0) <= B(i, 0)))
		{
			X2(i, 0) = X(i, 0) + dX(i, 0);
		}
		
		//Bellow the lower bound
		else if (X(i, 0) + dX(i, 0) < A(i, 0))
		{
			X2(i, 0) = A(i, 0);
		}

		//Over the upper bound
		else if ( X(i, 0) + dX(i, 0) > B(i, 0))
		{
			X2(i, 0) = B(i, 0);
		}
		
	}

	return X2;
}

#endif
