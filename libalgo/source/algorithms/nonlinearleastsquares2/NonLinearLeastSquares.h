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


#ifndef NonLinearLeastSquares_H
#define NonLinearLeastSquares_H

#include <set>
#include <ostream>
#include <iostream>

#include "libalgo/source/algorithms/matrixoperations2/MatrixOperations.h"


// Description: Non LinearLeast Squares algorithms
class NonLinearLeastSquares
{
        public:
               	template <typename T, typename FunctionJ, typename FunctionV>
		static T BFGSH(FunctionJ function_j, FunctionV function_v, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &A, const Matrix <T> &B, unsigned int &iterations, 
			const T alpha = 0.0001, const T nu = 0.0001, const T max_error = 1.0e-10, const unsigned int max_iterations = 50, const T max_diff = 1.0e-12, std::ostream * output = &std::cout);
		
		template <typename T>
		static void reflection(Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

};

#include "NonLinearLeastSquares.hpp"

#endif
