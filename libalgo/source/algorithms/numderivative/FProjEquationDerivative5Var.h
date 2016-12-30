// Description: Functor, compute numeric derivative of the projection equation using the Stirling formula
// Derivatives are: latp, lonp, lat0, lon0 (used for NLPS solution)

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


#ifndef FProjEquationDerivative5Var_H
#define FProjEquationDerivative5Var_H


#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"


//Functor, performs partial derivative of the projection equation using the Stirling formula
// Derivatives are: R, latp, lonp, lat0, lon0, c
template <typename T>
class FProjEquationDerivative5Var
{

        private:
                //Map projection parameters
		const TPostfixNotationDel * equation_postfix;
		const TPostfixNotationDel * ftheta_equat_postfix;
		const TPostfixNotationDel * theta0_equat_postfix;
		const T R;
                const T lat;
                const T lon;
                const T a;
                const T b;
                const T lat1;
                const T lat2;
                const TTransformedLongtitudeDirection trans_lon_dir;


        public:

		FProjEquationDerivative5Var(const TPostfixNotationDel * equation_postfix_, const TPostfixNotationDel * ftheta_equat_postfix_, const TPostfixNotationDel * theta0_equat_postfix_,  const T R_, const T lat_, const T lon_, const T a_, const T b_, const T lat1_, const T lat2_, const TTransformedLongtitudeDirection trans_lon_dir_) :
			equation_postfix(equation_postfix_), ftheta_equat_postfix(ftheta_equat_postfix_), theta0_equat_postfix(theta0_equat_postfix_),  R(R_), lat(lat_), lon(lon_), a(a_), b(b_), lat1(lat1_), lat2(lat2_), trans_lon_dir(trans_lon_dir_) {}

                T operator () ( const Matrix <T> &arg )
                {
                        //Reduce lon
                        const T lon_red = CartTransformation::redLon0 ( lon, arg ( 0, 3 ) );

                        //Convert ( lat, lon ) -> ( lat, lon)_trans
			const T lat_trans = CartTransformation::latToLatTrans(lat, lon_red, arg(0, 0), arg(0, 1));
                        const T lon_trans = CartTransformation::lonToLonTrans (lat, lon_red, arg ( 0, 0 ), arg ( 0, 1 ), trans_lon_dir );

                        //Compute partial derivative of the map projection equation
			T res = CartTransformation::latLonToCartesian(equation_postfix, ftheta_equat_postfix, theta0_equat_postfix, lat_trans, lon_trans, R, a, b, 0.0, arg(0, 4), arg(0, 2), arg(0, 2), arg(0, 4), false);
			//T res = CartTransformation::latLonToCartesian(equation, ftheta_equat, theta0_equat, lat_trans, lon_trans, R, a, b, 0.0, arg(0, 4), arg(0, 2), lat1, lat2, false);
			
			//T res =  ArithmeticParser::parseEquation ( equation, lat_trans, lon_trans, R, a, b, arg ( 0, 4 ), arg ( 0, 2 ), lat1, lat2, false );
                        return res;
                }

};

#endif
