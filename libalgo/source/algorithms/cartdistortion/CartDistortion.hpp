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


#ifndef CartDistortion_HPP
#define CartDistortion_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

//Include  correct version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/structures/projection/Projection.h"
#else
	#include "libalgo/source/structures/projection2/Projection.h"
#endif

#include "libalgo/source/algorithms/numderivative/NumDerivative.h"
#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative2Var.h"
#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative6Var.h"

#include "libalgo/source/exceptions/MathInvalidArgumentException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
T CartDistortion::H ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute distortion h (length distortion in meridian) for point p = [lat, lon] in cartographic projection
	return H(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                   proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::K ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute distortion k (length distortion in parallel) for point p = [lat, lon] in cartographic projection
	return K(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),   p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                   proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::Theta ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute angle between meridian and parallel in point p = [lat, lon]
	return Theta(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                       proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::S ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute distortion S (aerial distortion) for point p = [lat, lon] in cartographic projection
	return S(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                   proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::P ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute P for point p = [lat, lon] in cartographic projection
	return P(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                   proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::BM ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute bearing of the meridian given by the point p = [lat, lon] in cartographic projection
	return BM(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                    proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::BP ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute bearing of the parallel given by the point p = [lat, lon] in cartographic projection
	return BP(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),   p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                    proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
TTissotIndicatrix <T> CartDistortion::Tiss ( const T step, const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exceptions )
{
        //Compute major, semi-major axis, rotation angle of Tissot indicatrix and bearing of the meridian in point p = [lat, lon]
	return Tiss(step, proj->getXEquatPostfix(), proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
                      proj->getDy(), proj->getC(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), print_exceptions );
}


template <typename T>
T CartDistortion::H(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute distortion h (length distortion in meridian) for point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);

        //Distortion
        return H ( R, dx_dlat, dy_dlat );
}


template <typename T>
T CartDistortion:: H ( const T R, const T dx_dlat, const T dy_dlat )
{
        //Compute distortion h (length distortion in meridian) for point p = [lat, lon] in cartographic projection
        if ( R == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion h,", "R = ", R );
        }

        //Distortion
        return sqrt ( ( dx_dlat * dx_dlat + dy_dlat * dy_dlat ) ) / R;
}


template <typename T>
T CartDistortion::K(const T step, const TPostfixNotationDel * equation_x_postfix, const TPostfixNotationDel * equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute distortion k (length distortion in parallel) for point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinates functions
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        //Compute distortion
        return K ( lat, R, dx_dlon, dy_dlon );
}


template <typename T>
T CartDistortion:: K ( const T lat, const T R, const T dx_dlon, const T dy_dlon )
{
        //Compute distortion k (length distortion in parallel) for point p = [lat, lon] in cartographic projection

        //Do not divide by 0
        if ( R == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion k, ", "R = ", R );
        }

        //Throw exception
        if ( ( lat > MAX_LAT - ANGLE_ROUND_ERROR ) || ( lat < MIN_LAT + ANGLE_ROUND_ERROR ) )
        {
                return MAX_FLOAT;
        }

        //Compute distortion
        return sqrt ( ( dx_dlon * dx_dlon  + dy_dlon * dy_dlon ) ) / ( R * cos ( lat * M_PI / 180 ) );
}


template <typename T>
T CartDistortion::Theta(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute angle between meridian and parallel in point p = [lat, lon]
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        return Theta ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );
}


template <typename T>
T CartDistortion:: Theta ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon )
{
        //Compute angle between meridian and parallel in point p = [lat, lon]

        //Do not divide by 0
        if ( R == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion theta, ", "R = ", R );
        }

        //Throw exception
        if ( ( lat > MAX_LAT - ANGLE_ROUND_ERROR ) || ( lat < MIN_LAT + ANGLE_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute cartographic distortion theta, ", "fabs(lat) > 90, lat=", lat );
        }

        //Compute denominator
        const T denom = ( dx_dlat * dx_dlat + dy_dlat * dy_dlat ) * ( dx_dlon * dx_dlon + dy_dlon * dy_dlon );

        //Throw exception
        if ( denom == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion theta, ", "fabs (denom) = 0", denom );
        }

        //Compute asin
        T theta_asin = ( dx_dlat * dy_dlon - dx_dlon * dy_dlat ) / sqrt ( denom );

        //Throw exception
        if ( ( theta_asin > 1.0 + ARGUMENT_ROUND_ERROR ) || ( theta_asin < -1.0 - ARGUMENT_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute cartographic distortion theta, ", "asin(arg), arg = ", theta_asin );
        }

        //Correct longitude
        else if ( theta_asin > 1.0 )
        {
                theta_asin = 1.0;
        }

        //Correct longitude
        else if ( theta_asin < -1.0 )
        {
                theta_asin = -1.0;
        }

        //Compute theta
        return  fabs ( asin ( theta_asin ) * 180.0 / M_PI );
}


template <typename T>
T CartDistortion::S(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute distortion S (aerial distortion) for point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        //Compute distortion S
        return S ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );
}



template <typename T>
T CartDistortion::S ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon )
{
        //Compute distortion S (aerial distortion) for point p = [lat, lon] in cartographic projection

        //Throw exception
        if ( R == 0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion theta, ", "R = 0" );
        }

        //Throw exception
        if ( ( lat > MAX_LAT - ANGLE_ROUND_ERROR ) || ( lat < MIN_LAT + ANGLE_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute cartographic distortion s, ", "fabs(lat) >= 90, lat=", lat );
        }

        //Compute aerial distortion s
        return fabs ( ( dx_dlat * dy_dlon - dx_dlon * dy_dlat ) / ( R * R * cos ( lat * M_PI / 180 ) ) );
}


template <typename T>
T CartDistortion::P(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute P for point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        //Comnpute P
        return P ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );
}


template <typename T>
T CartDistortion::P ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon )
{
        //Compute P for point p = [lat, lon] in cartographic projection

        //Throw exception
        if ( R == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortion theta, ", "R = ", R );
        }

        //Throw exception
        if ( ( lat > MAX_LAT - ANGLE_ROUND_ERROR ) || ( lat < MIN_LAT + ANGLE_ROUND_ERROR ) )
        {
                return MAX_FLOAT;
        }

        //Compute aerial distortion s
        return  2 * ( dx_dlat * dx_dlon + dy_dlat * dy_dlon ) / ( R * R * cos ( lat * M_PI / 180.0 ) );
}


template <typename T>
T CartDistortion::BM(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute bearing of the meridian given by the point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);

        return BM ( dx_dlat, dy_dlat );
}


template <typename T>
T CartDistortion:: BM ( const T dx_dlat, const T dy_dlat )
{
        //Compute bearing of the meridian given by the point p = [lat, lon] in cartographic projection

        //Throw exception
        if ( dx_dlat == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic bearing of the meridian, ", "dx_dlat = ", dx_dlat );
        }

        if ( dy_dlat == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic bearing of the meridian, ", "dx_dlat = ", dy_dlat );
        }

        T bm = atan2 ( dy_dlat, dx_dlat ) * 180.0 / M_PI;

        //Convert interval (-Pi, Pi) to (0, 2 * PI)
        if ( bm < 0.0 )
        {
                return bm + 360.0;
        }

        return bm;
}


template <typename T>
T CartDistortion::BP(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute bearing of the parallel given by the point p = [lat, lon] in cartographic projection
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        //Throw exception
        if ( ( dx_dlon == 0.0 ) && ( dy_dlon == 0.0 ) )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic bearing of the meridian, ", "dx_dlat = 0, dy_dlat = 0." );
        }

        return BP ( dx_dlon, dy_dlon );
}


template <typename T>
T CartDistortion:: BP ( const T dx_dlon, const T dy_dlon )
{
        //Compute bearing of the parallel given by the point p = [lat, lon] in cartographic projection

        //Throw exception
        if ( ( dx_dlon == 0.0 ) && ( dy_dlon == 0.0 ) )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic bearing of the meridian, ", "dx_dlat = 0, dy_dlat = 0." );
        }

        T bp = atan2 ( dy_dlon, dx_dlon ) * 180 / M_PI;

        //Convert interval (-Pi, Pi) to (0, 2 * PI)
        if ( bp < 0 )
        {
                return bp + 360.0;
        }

        return bp;
}


template <typename T>
TTissotIndicatrix <T> CartDistortion::Tiss(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel * equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute major, semi-major axis, rotation angle of Tissot indicatrix and bearing of the meridian in point p = [lat, lon]
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        return Tiss ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );
}


template <typename T>
TTissotIndicatrix <T> CartDistortion:: Tiss ( const T lat, const T R, const T dx_dlat, const T dx_dlon, const T dy_dlat, const T dy_dlon )
{
        //Compute major, semi-major axis, rotation angle of Tissot indicatrix and bearing of the meridian in point p = [lat, lon]

        //Get distortions in analyzed point
        const T h = H ( R, dx_dlat, dy_dlat );
        const T k = K ( lat, R, dx_dlon, dy_dlon );
        const T p = P ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );

        //Get angle meridian-parallel in analyzed point
        const T theta = Theta ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );

        //Get bearing of the meridian in analyzed point
        const T b_mer = BM ( dx_dlat, dy_dlat );

        //Compute c^2, d^2
        T c_sqr = h * h + k * k + 2 * h * k * sin ( theta * M_PI / 180 );
        T d_sqr = h * h + k * k - 2 * h * k * sin ( theta * M_PI / 180 );

        //Throw exception
        if ( c_sqr < - ARGUMENT_ROUND_ERROR )
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute Tissot Indicatrix parameters, ", "sqrt(c) < 0, c = ", c_sqr );

        if ( d_sqr < - ARGUMENT_ROUND_ERROR )
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute Tissot Indicatrix parameters, ", "sqrt(d) < 0, d = ", d_sqr );

        //Correct values: round errors
        if ( c_sqr < 0.0 )
                c_sqr = 0.0;

        if ( d_sqr < 0.0 )
                d_sqr = 0.0;

        //Compute extremal azimuth
        const T denom = h * h - k * k;

        T Ae = 0.0;

        if ( denom == 0.0 )
        {
                if ( p > 0.0 ) Ae = 90.0;
                else if ( p < 0.0 ) Ae = -90.0;
                else MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute Tissot Indicatrix parameters, ", "atan(0/0) < 0, d = ", p );
        }

        else
        {
                Ae = 0.5 * ( atan ( p / denom ) * 180 / M_PI );
        }

        //Test if found extremal azimuth belongs to max distortion or min distortion
        const T m1 = h * h * cos ( Ae * M_PI / 180 ) * cos ( Ae * M_PI / 180 ) + k * k * sin ( Ae * M_PI / 180 ) * sin ( Ae * M_PI / 180 ) +
                     p * sin ( Ae * M_PI / 180 ) * cos ( Ae * M_PI / 180 );
        const T m2 = h * h * cos ( ( Ae + 90 ) * M_PI / 180 ) * cos ( ( Ae + 90 ) * M_PI / 180 ) + k * k * sin ( ( Ae + 90 ) * M_PI / 180 ) * sin ( ( Ae + 90 ) * M_PI / 180 ) +
                     p * sin ( ( Ae + 90 ) * M_PI / 180 ) * cos ( ( Ae + 90 ) * M_PI / 180 );

        //If belongs to min distortion, correct azimuth: add Pi/2
        if ( m2 > m1 )
                Ae += 90.0;

        //Compute major and minor axis ot the Tissot Indicatrix
        const T a_tiss = ( sqrt ( c_sqr ) + sqrt ( d_sqr ) ) / 2.0;
        const T b_tiss = ( sqrt ( c_sqr ) - sqrt ( d_sqr ) ) / 2.0;

        //Compute extremal azimuth projected
        if ( a_tiss == 0.0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute azimuth of Tissot Indicatrix, ", "fabs (denom) = 0", a_tiss );
        }

        const T Ae_proj = atan ( b_tiss / a_tiss * tan ( Ae * M_PI / 180.0 ) ) * 180.0 / M_PI;

        //Create struct [a, b, Ae, Ae_proj, bm]
        TTissotIndicatrix <T> tiss ( a_tiss, b_tiss, Ae, Ae_proj, b_mer );

        return tiss;
}


template <typename T>
T CartDistortion::W(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a_, const T b_, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute max angular distortion
	const TTissotIndicatrix <T> tiss = Tiss(step, equation_x_postfix, equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix, lat, lon, R, a_, b_, dx, dy, lat0, lat1, lat2, lon0);

        //Throw exception
        if ( ( tiss.a_tiss == 0.0 ) && ( tiss.b_tiss == 0.0 ) )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute cartographic distortions a, b,  ", " a = 0, b = 0" );
        }

        //Max angle distortion
        T omega_sin = fabs ( tiss.b_tiss - tiss.a_tiss ) / ( tiss.a_tiss + tiss.b_tiss );

        //Throw exception
        if ( ( omega_sin > 1 + ARGUMENT_ROUND_ERROR ) || ( omega_sin < -1 - ARGUMENT_ROUND_ERROR ) )
        {
                throw MathInvalidArgumentException <T> ( "MathInvalidArgumentException: can not compute cartographic distortion omega, ", "asin(arg)", omega_sin );
        }

        //Correct value
        if ( omega_sin > 1.0 )
        {
                omega_sin = 1.0;
        }

        //Correct value
        else if ( omega_sin < -1.0 )
        {
                omega_sin = -1.0;
        }

        //Compute max angle distortion
        return  asin ( omega_sin ) * 2 * M_PI / 180;
}


template <typename T>
T CartDistortion::Airy(const T step, const TPostfixNotationDel *equation_x_postfix, const TPostfixNotationDel *equation_y_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T dy, const T c, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exceptions)
{
        //Compute local Airy criterion
        Matrix <double> args ( 1, 2 );
        args ( 0, 0 ) = lat; args ( 0, 1 ) = lon;

        //Partial derivative of coordinate functions
	const T dx_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  lat, lon, R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dy_dlat = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  lat, lon, R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX1, step, print_exceptions);
	const T dx_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix,  lat, lon, R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);
	const T dy_dlon = 180.0 / M_PI * NumDerivative::getDerivative(FProjEquationDerivative2Var <T>(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix,  lat, lon, R, a, b, dx, dy, c, lat0, lat1, lat2, lon0), args, FirstDerivative, VariableX2, step, print_exceptions);

        //Compute parameters of a Tissot indicatrix
        TTissotIndicatrix <T> tissot;
        tissot =  Tiss ( lat, R, dx_dlat, dx_dlon, dy_dlat, dy_dlon );

        return sqrt ( 0.5 * ( ( tissot.a_tiss - 1.0 ) * ( tissot.a_tiss - 1.0 ) + ( tissot.b_tiss - 1.0 ) * ( tissot.b_tiss - 1.0 ) ) );
}



#endif
