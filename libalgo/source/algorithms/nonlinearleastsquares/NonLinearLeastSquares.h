// Description: Non Linear Least Squares algorithms

// Copyright (c) 2010 - 2016
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


#ifndef NonLinearLeastSquares_H
#define NonLinearLeastSquares_H

#include <set>
#include <ostream>
#include <iostream>

#include "libalgo/source/algorithms/matrixoperations/MatrixOperations.h"



// Description: Non LinearLeast Squares algorithms
class NonLinearLeastSquares
{
        public:
                template <typename T, typename FunctionA, typename FunctionV, typename FunctionC>
                static T GN ( FunctionA function_a, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations );

                template <typename T, typename FunctionA, typename FunctionV, typename FunctionC>
		static T GND(FunctionA function_a, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T alpha = 0.0001, const T max_error = 1.0e-10, const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, std::ostream * output = &std::cout);

		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T BFGS(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T alpha = 0.0001, const T max_error = 1.0e-10, const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, std::ostream * output = &std::cout);

		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T BFGSH(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T alpha = 0.0001, const T nu = 0.0001, const T max_error = 1.0e-10, const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, std::ostream * output = &std::cout);

		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T LM(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T nu1 = 0.25, const T nu2 = 0.75, const T nu3 = 0.0001, const T gamma1 = 0.5, const T gamma2 = 2.0, const T lambda_min = 1.0e-6, const  T lambda_max = 1.0e6, const T max_error = 1.0e-10, 
			const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, T lambda = 1, std::ostream * output = &std::cout);
		
		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T TR(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T nu1 = 0.25, const T nu2 = 0.75, const T nu3 = 0.0001, const T gamma1 = 0.25, const T gamma2 = 2.5, const T delta_min = 1.0e-6, const  T delta_max = 1.0e6, const T max_error = 1.0e-10, 
			const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, T delta = 1.0, std::ostream * output = &std::cout);

		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T TR2(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T nu1 = 0.25, const T nu2 = 0.75, const T nu3 = 0.0001, const T gamma1 = 0.25, const T gamma2 = 2.5, const T delta_min = 1.0e-6, const  T delta_max = 1.0e6, const T max_error = 1.0e-10, 
			const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, T delta = 1.0, std::ostream * output = &std::cout);

		template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
		static T TR3(FunctionJ function_j, FunctionV function_v, FunctionC function_c, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned short &iterations, 
			const T nu1 = 0.25, const T nu2 = 0.75, const T nu3 = 0.0001, const T gamma1 = 0.25, const T gamma2 = 2.5, const T delta_min = 1.0e-6, const  T delta_max = 1.0e6, const T max_error = 1.0e-10, 
			const unsigned short max_iterations = 100, const T max_diff = 1.0e-12, T delta = 1.0, std::ostream * output = &std::cout);

		template <typename T>
		static void reflection(Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

	private:
		template <typename T>
		static Matrix <T> optTRStep(const Matrix <T> &J, const Matrix <T> &V, const Matrix <T> &W, const T delta, const T max_ter = 50, const T max_error = ARGUMENT_ROUND_ERROR );

		template <typename T>
		static Matrix <T> optTRStep2(const Matrix <T> &J, const Matrix <T> &V, const Matrix <T> &W, const T delta, const T max_iter, const T max_error = ARGUMENT_ROUND_ERROR);

		template <typename T>
		static void trStep32(const Matrix <T> &B, Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min, bool &stop);

		template <typename T>
		static void trStep4(const Matrix <T> &B, const Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min, bool &stop);

		template <typename T>
		static void trStep42(const Matrix <T> &B, const Matrix <T> &R, Matrix <T> &s, const Matrix <T> &G, Matrix <T> &v, T & lambda_min, T &lambda_max, T &lambda, const T delta_min, const T delta_max, const T delta, T &mju_min, bool &stop);

		template <typename T>
		static void trStep5 (const Matrix <T> &B, const Matrix <T> &R, Matrix<T> &s, Matrix<T> &v, T &lambda, const T delta_min, const T delta, T &mju_min, bool &stop);
		
		template <typename T>
		static void trStep52(const Matrix <T> &B, const Matrix <T> &R, Matrix<T> &s, const Matrix <T> &G, Matrix<T> &v, const T lambda_min, const T lambda_max, T &lambda, const T delta_min, const T delta, T &mju_min, bool &stop);
		
		template <typename T>
		static T trStep6(const Matrix <T> &R, const Matrix <T> &s, Matrix<T> &v, const T delta);
		
		template <typename T>
		static void trStep62(const Matrix <T> &R, const Matrix <T> &s, const Matrix <T> &G, Matrix<T> &v, const T lambda_min, const T lambda_max, T &lambda, const T delta, const T mju_min);

		template <typename T>
		static void reflection2(Matrix <T> &X, Matrix <T> &dX, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

		template <typename T>
		static void reflection3(Matrix <T> &X, Matrix <T> &dX, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

		template <typename T>
		static bool checkIntervals(Matrix <T> &X, const Matrix <T> &A, const Matrix <T> &B);

		template <typename T>
		static Matrix <T> addStep(const Matrix <T> &X, const Matrix <T> &dX, const Matrix <T> &A, const Matrix <T> &B);


};

#include "NonLinearLeastSquares.hpp"

#endif
