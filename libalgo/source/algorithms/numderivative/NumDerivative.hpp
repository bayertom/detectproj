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


#ifndef NumDerivative_HPP
#define NumDerivative_HPP

#include "libalgo/source/exceptions/BadDataException.h"


template <typename T, typename Function>
T NumDerivative::getDerivative(Function function, const Matrix <T> &args, const TDerivativeType deriv_type, const TDerivativeVariable deriv_var, const T deriv_step, const bool print_exceptions)
{
        //Values and function values
        T values [7], fvalues[7];

        //Compute function values
        computeFunctionValues ( function, args, deriv_var, deriv_step, values, fvalues, print_exceptions );

        //Compute differences and Stirling formula
        return computeStirlingFormula ( values, fvalues, deriv_type, deriv_step );
}


template <typename T, typename Function>
void NumDerivative::computeFunctionValues(Function function, const Matrix <T> &args, const TDerivativeVariable deriv_var, const T deriv_step, T * values, T * fvalues, const bool print_exceptions)
{
        //Compute values and function values of the function having one variable z = f(lat, lon)
        if ( ( unsigned int ) deriv_var > args.cols() -  1 )
                throw BadDataException ( "BadDataException: not enough arguments (deriv_index > arg.size() - 1).", "Can not compute the partial derivative..." );

        for ( unsigned short i = 0; i < 7; i++ )
        {
                //Compute value from argument
                values[i] = T ( args ( 0, deriv_var ) + ( i - 2 ) * deriv_step );

                //Create temporary argument, replace arg[deriv_index] with the value
                Matrix <T> args_temp = args;
                args_temp ( 0, deriv_var ) = values[i];

                //Compute the function value
                fvalues[i] =  function ( args_temp );

		//std::cout << values[i] << " " << fvalues[i] << '\n';
		
        }
}


template <typename T>
T NumDerivative::computeStirlingFormula(T * values, T * fvalues, const TDerivativeType deriv_type, const T deriv_step)
{
        // First differences
        T dif1[6];

        for ( unsigned short i = 0; i < 6; i++ )
        {
                dif1[i] = fvalues[i + 1] - fvalues[i];
        }

        // Second differences
        T dif2[5];

        for ( unsigned short i = 0; i < 5; i++ )
        {
                dif2[i] = dif1[i + 1] - dif1[i];
        }

        // Third differences
        T dif3[4];

        for ( unsigned short i = 0; i < 4; i++ )
        {
                dif3[i] = dif2[i + 1] - dif2[i];
        }

	// Fourth differences
	T dif4[3];

	for (unsigned short i = 0; i < 3; i++)
	{
		dif4[i] = dif3[i + 1] - dif3[i];
	}

	// Fifth differences
	T dif5[2];

	for (unsigned short i = 0; i < 2; i++)
	{
		dif5[i] = dif4[i + 1] - dif4[i];
	}

        //Stirling formula for numeric derivative
	if (deriv_type == FirstDerivative)

		//First derivative
		return ( (dif1[2] + dif1[3]) / 2 - (dif3[1] + dif3[2]) / 12 +  (dif5[0] + dif5[1]) / 60 ) / deriv_step;
	else

		//Second derivative
		return (dif2[2]  - dif4[1] / 12 + (dif5[1] - dif5[0]) / 90) / ( deriv_step * deriv_step );
}


#endif
