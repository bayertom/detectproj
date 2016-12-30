// Description: Functor, compute numeric derivative of additional parameter p of the projection equation using the Stirling formula
// Derivatives is: p

// Copyright (c) 2010 - 2015
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


#ifndef FProjEquationDerivative1Var_H
#define FProjEquationDerivative1Var_H


#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/arithmeticparser/ArithmeticParser.h"


//Functor, performs partial derivative of the projection equation using the Stirling formula
// Derivatives are: lat lonp
template <typename T>
class FProjEquationDerivative1Var
{

private:
	//Map projection parameters
	const TPostfixNotationDel *ftheta_equat_postfix;
	const T lat;
	const T lon;
	const T R;
	const T a;
	const T b;
	const T dx;
	const T dy;
	const T c;
	const T lat0;
	const T lat1;
	const T lat2;

public:

	FProjEquationDerivative1Var(const TPostfixNotationDel * ftheta_equat_postfix_, const T lat_, const T lon_, const T R_, const T a_, const T b_, const T dx_, const T dy_, const T c_, const T lat0_, const T lat1_, const T lat2_) :
		ftheta_equat_postfix(ftheta_equat_postfix_), lat(lat_), lon(lon_), R(R_), a(a_), b(b_), dx(dx_), dy(dy_), c(c_), lat0(lat0_), lat1(lat1_), lat2(lat2_) {}

                T operator () ( const Matrix <T> &arg )
                {
			return ArithmeticParser::parseEquation(ftheta_equat_postfix, lat, lon, R, a, b, c, lat0, lat1, lat2, arg(0, 0), false);
                }

};

#endif
