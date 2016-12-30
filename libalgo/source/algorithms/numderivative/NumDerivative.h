// Description: Compute numeric derivative using Stirling method, the function is defined using the functor
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


#ifndef NumDerivative_H
#define NumDerivative_H

#include "libalgo/source/types/TDerivativeType.h"
#include "libalgo/source/types/TDerivativeVariable.h"

//Include  C++98/C++11 version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/structures/matrix/Matrix.h"
#else
	#include "libalgo/source/structures/matrix2/Matrix.h"
#endif


//Compute numerical derivative using Stirling method
class NumDerivative
{
	
        public:

                template <typename T, typename Function>
		static T getDerivative(Function function, const Matrix <T> &args, const TDerivativeType deriv_type, const TDerivativeVariable deriv_var, const T deriv_step, const bool print_exceptions = false);


        private:

                template <typename T, typename Function>
		static void computeFunctionValues(Function function, const Matrix <T> &args, const TDerivativeVariable deriv_var, const T deriv_step, T * values, T * fvalues, const bool print_exceptions = false);


                template <typename T>
		static T computeStirlingFormula(T * values, T * fvalues, const TDerivativeType deriv_type, const T deriv_step);

};

#include "NumDerivative.hpp"

#endif
