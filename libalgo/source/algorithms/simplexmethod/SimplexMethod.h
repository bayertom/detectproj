// Description: Downhill simplex optimization method

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


#ifndef SimplexMethod_H
#define SimplexMethod_H

//Include  C++98/C++11 version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/algorithms/matrixoperations/MatrixOperations.h"
#else
	#include "libalgo/source/algorithms/matrixoperations2/MatrixOperations.h"
#endif


//Downhill simplex optimization method
class SimplexMethod
{
        public:

                template <typename T, typename Function>
                static T NelderMead ( Function function, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &XMIN, const Matrix <T> &XMAX, unsigned int &iterations, const T max_error, const unsigned int max_iterations = 600, const bool add_x0 = false, std::ostream * output = &std::cout );

        private:

                template <typename T>
                static void createRandSimplex ( const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> & XX, const bool add_x0 = false);

                template <typename T, typename Function>
                static void shrink ( Function function, Matrix <T> &W, Matrix <T> &X, const Matrix<T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &Y,  Matrix <T> &V, const T Sigma );

                template <typename T>
                static void reflection (Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

};

#include "SimplexMethod.hpp"

#endif
