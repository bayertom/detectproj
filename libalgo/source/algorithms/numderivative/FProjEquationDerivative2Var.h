// Description: Functor, compute numeric derivative of projection equation using the Stirling formula
// Derivatives are: lat, lon

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


#ifndef FProjEquationDerivative2Var_H
#define FProjEquationDerivative2Var_H


#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"



//Functor, performs partial derivative of the projection equation using the Stirling formula
// Derivatives are: lat lonp
template <typename T>
class FProjEquationDerivative2Var
{

        private:
                //Map projection parameters
		const TPostfixNotationDel * equation_postfix;
		const TPostfixNotationDel * ftheta_equat_postfix;
		const TPostfixNotationDel  *theta0_equat_postfix;
                const T R;
                const T a;
                const T b;
                const T dx;
                const T dy;
                const T c;
                const T lat0;
                const T lat1;
                const T lat2;
                const T lon0;

        public:

		FProjEquationDerivative2Var(const TPostfixNotationDel * equation_postfix_, const TPostfixNotationDel * ftheta_equat_postfix_, const TPostfixNotationDel * theta0_equat_postfix_, const T R_, const T a_, const T b_, const T dx_, const T dy_, const T c_, const T lat0_, const T lat1_, const T lat2_, const T lon0_) :
			equation_postfix(equation_postfix_), ftheta_equat_postfix(ftheta_equat_postfix_), theta0_equat_postfix(theta0_equat_postfix_), R(R_), a(a_), b(b_), dx(dx_), dy(dy_), c(c_), lat0(lat0_), lat1(lat1_), lat2(lat2_), lon0(lon0_) {}

                T operator () ( const Matrix <T> &arg )
                {
                        //Compute partial derivative of the map projection equation
			return CartTransformation::latLonToCartesian(equation_postfix, ftheta_equat_postfix, theta0_equat_postfix, arg(0, 0), arg(0, 1), R, a, b, 0.0, c, lat0, lat1, lat2, false);
                        //return ArithmeticParser::parseEquation ( equation, arg ( 0, 0 ), arg ( 0, 1 ), R, a, b, c, lat0, lat1, lat2, false );
                }

};

#endif
