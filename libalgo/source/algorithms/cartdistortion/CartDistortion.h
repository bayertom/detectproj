// Description: Compute cartographic distortions, length, area and angle distortions, Tissot�s indicatrix

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


#ifndef CartDistortion_H
#define CartDistortion_H

#include "libalgo/source/algorithms/arithmeticparser/ArithmeticParser.h"

//Forward declaration
template <typename T>
class Point3DGeographic;

template <typename T>
class Projection;

template <typename T>
class FProjEquationDerivative6Var;

template <typename T>
class FProjEquationDerivative2Var;

//Structure storing parameters of the Tissot Indicatrix
template <typename T>
struct TTissotIndicatrix
{
        T a_tiss, b_tiss, Ae, Ae_proj, b_mer;

        TTissotIndicatrix () : a_tiss ( 1.0 ), b_tiss ( 1.0 ), Ae ( 0.0 ), Ae_proj ( 0.0 ), b_mer ( 0.0 ) {}
        TTissotIndicatrix ( const T a_tiss_, const T b_tiss_, const T Ae_, const T Ae_proj_, const T b_mer_ ) :
                a_tiss ( a_tiss_ ), b_tiss ( b_tiss_ ), Ae ( Ae_ ), Ae_proj ( Ae_proj_ ), b_mer ( b_mer_ ) {}
};


//Projectedgraphic distortion
class CartDistortion
{
        public:
                template <typename T>
                static T H ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T K ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T Theta ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T S ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T P ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T BM ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static T BP ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

                template <typename T>
                static TTissotIndicatrix <T> Tiss ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions = false );

        public:
                template <typename T>
		static T H(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T K(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T Theta(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T S(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T P(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T BM(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T BP(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static TTissotIndicatrix <T> Tiss(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T W(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a_, const T b_, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

                template <typename T>
		static T Airy(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions = false);

        private:

                template <typename T>
                static T H ( const T R, const T dx_dlat, const T dy_dlat );

                template <typename T>
                static T K ( const T lat, const T R, const T dx_dlon, const T dy_dlon );

                template <typename T>
                static T S ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon );

                template <typename T>
                static T P ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon );

                template <typename T>
                static T BM ( const T dx_dlat, const T dy_dlat );

                template <typename T>
                static T BP ( const T dx_dlon, const T dy_dlon );

                template <typename T>
                static T Theta ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon );

                template <typename T>
                static TTissotIndicatrix <T> Tiss ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon );
};

#include "CartDistortion.hpp"

#endif
