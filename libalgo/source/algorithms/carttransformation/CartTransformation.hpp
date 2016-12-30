// Description: transformation of the geographic coordinates (lat, lon) to (lat, lon)_trans. Conversion (lat, lon) to (x, y).

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

#ifndef CartTransformation_HPP
#define CartTransformation_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/point/Point3DCartesianProjected.h"
#include "libalgo/source/structures/projection/Projection.h"

#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative1Var.h"
#include "libalgo/source/algorithms/numderivative/NumDerivative.h"

#include "libalgo/source/exceptions/MathInvalidArgumentException.h"
#include "libalgo/source/exceptions/MathOverflowException.h"
#include "libalgo/source/exceptions/ParseException.h"
#include "libalgo/source/exceptions/BadDataException.h"

template <typename T>
T CartTransformation::latToLatTrans(const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole)
{
	//Transform latitude  ( lat, lon ) -> ( lat_tans ) using a cartographic pole (latp, lonp)
	return latToLatTrans(p->getLat(), p->getLon(), pole->getLat(), pole->getLon());
}

template <typename T>
T CartTransformation::lonToLonTrans(const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole, const TTransformedLongitudeDirection lon_direction)
{
	//Transform longitude  ( lat, lon, lat_trans ) -> ( lon_tans ) using a cartographic pole (latp, lonp)
	return lonToLonTrans(p->getLat(), p->getLon(), pole->getLat(), pole->getLon(), lon_direction);
}

template <typename T>
T CartTransformation::latLonToX(const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception)
{
	//Compute x coordinate of the point P[lat, lon] in specific projection
	//Projection equations in the postfix notation
	return latLonToCartesian(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getC(),
		proj->getLat0(), proj->getLat1(), proj->getLat2(), print_exception);
}

template <typename T>
T CartTransformation::latLonToY(const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception)
{
	//Compute y coordinate of the point P[lat, lon] in specific projection
	//Projection equations in the postfix notation
	return latLonToCartesian(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDy(), proj->getC(),
		proj->getLat0(), proj->getLat1(), proj->getLat2(), print_exception);
}

template <typename T>
T CartTransformation::latLonToX(TPostfixNotationDel * equation_x_postfix, TPostfixNotationDel * ftheta_equat_postfix, TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dx, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception)
{
	//Compute x coordinate of the point P[lat, lon] in specific projection
	//Projection equations in the postfix notation
	return latLonToCartesian(equation_x_postfix, ftheta_equat_postfix, equation_theta0_postfix, lat, lon, R, a, b, dx, c, lat0, lat1, lat2, print_exception);
}

template <typename T>
T CartTransformation::latLonToY(TPostfixNotationDel * equation_y_postfix, TPostfixNotationDel * ftheta_equat_postfix, TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T dy, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception)
{
	//Compute y coordinate of the point P[lat, lon] in specific projection
	//Projection equations in the postfix notation
	return latLonToCartesian(equation_y_postfix, ftheta_equat_postfix, equation_theta0_postfix, lat, lon, R, a, b, dy, c, lat0, lat1, lat2, print_exception);
}

template <typename T>
T CartTransformation::latToLatTrans(const T lat, const T lon, const T latp, const T lonp)
{
	//Transform latitude  ( lat, lon ) -> ( lat_tans ) using a cartographic pole (latp, lonp)

	//Throw exception: bad lat
	if (fabs(lat) > MAX_LAT)
	{
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: ", "can not convert lat to lat_trans, lat > +- Pi/2", lat);
	}

	//Throw exception: bad lon
	if (fabs(lon) > MAX_LON)
	{
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: ", "can not convert lon to lon_trans, lon > +- Pi", lon);
	}

	//Projection in normal aspect
	if (fabs(MAX_LAT - latp) < ANGLE_ROUND_ERROR)
	{
		return lat;
	}

	//Same coordinates as the cartographic pole, singular point
	if ((fabs(lon - lonp) < ANGLE_ROUND_ERROR) && (fabs(lat - latp) < ANGLE_ROUND_ERROR))
	{
		return MAX_LAT;
	}

	//Compute latitude
	T lat_trans_asin = sin(lat * M_PI / 180.0) * sin(latp * M_PI / 180.0) + cos(lat * M_PI / 180.0) * cos(latp * M_PI / 180.0) * cos((lonp - lon) * M_PI / 180.0);

	//Throw exception
	if ((lat_trans_asin > 1.0 + ARGUMENT_ROUND_ERROR) || (lat_trans_asin < -1.0 - ARGUMENT_ROUND_ERROR))
	{
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: ", "can not convert lat to lat_trans, asin(arg), arg = ", lat_trans_asin);
	}

	//Correct latitude
	if (lat_trans_asin > 1.0)
	{
		return MAX_LAT;
	}

	//Correct latitude
	if (lat_trans_asin < -1.0)
	{
		return MIN_LAT;
	}

	//Compute transformed latitude
	return  asin(lat_trans_asin) * 180.0 / M_PI;
}

template <typename T>
T CartTransformation::lonToLonTrans(const T lat, const T lon, const T latp, const T lonp, const TTransformedLongitudeDirection lon_direction)
{
	//Transform longitude  ( lat, lon, lat_trans ) -> ( lon_tans ) using a cartographic pole (latp, lonp)

	//Throw exception: bad lat
	if (fabs(lat) > MAX_LAT)
	{
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: ", " can not convert lon to lon_trans, lat > +- Pi/2", lat);
	}

	//Throw exception: bad lon
	if (fabs(lon) > MAX_LON)
	{
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: ", " can not convert lon to lon_trans, lon > +- Pi", lon);
	}

	//Projection in normal position
	if (fabs(MAX_LAT - latp) < ANGLE_ROUND_ERROR)
	{
		return lon;
	}

	//Same coordinates as the cartographic pole: singular point
	if ((fabs(lon - lonp) < ANGLE_ROUND_ERROR) && (fabs(lat - latp) < ANGLE_ROUND_ERROR))
	{
		return lon;
	}

	//Compute lon_trans: Normal direction
	T lon_trans = atan2(cos(lat * M_PI / 180) * sin((lon - lonp) * M_PI / 180), cos((lon - lonp) * M_PI / 180) * sin(latp * M_PI / 180) * cos(lat * M_PI / 180) - sin(lat * M_PI / 180) * cos(latp * M_PI / 180)) * 180 / M_PI;

	//Reversed direction (JTSK)
	if (lon_direction == ReversedDirection)
		lon_trans = -lon_trans;

	//Normal direction 2
	else if (lon_direction == NormalDirection2)
	{
		if (lon_trans < 0)
			lon_trans = lon_trans + 180;
		else
			lon_trans = lon_trans - 180;
	}

	//Reversed direction 2
	else if (lon_direction == ReversedDirection2)
	{
		if (lon_trans < 0)
			lon_trans = -180 - lon_trans;
		else
			lon_trans = 180 - lon_trans;
	}

	return lon_trans;
}

template <typename T>
T CartTransformation::latLonToCartesian(const TPostfixNotationDel * equation_postfix, const TPostfixNotationDel * ftheta_equat_postfix, const TPostfixNotationDel * equation_theta0_postfix, const T lat, const T lon, const T R, const T a, const T b, const T shift, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception)
{
	//Compute cartesian coordinate x / y of the point P[lat, lon] in specific projection,
	//Projection equations in the postfix notation
	//Additional equations solved by Newton-Raphson method are included
	T theta = lat;

	//Apply Newton-Raphson method to evaluate p
	if ((ftheta_equat_postfix != NULL) && (equation_theta0_postfix != NULL))
	{
		bool convergence = false;
		unsigned int iterations = 0;
		const unsigned int max_iterations = 10;
		const T max_diff = 1.0e-5;

		//Determine p[0]
		theta = ArithmeticParser::parseEquation(equation_theta0_postfix, lat, lon, R, a, b, c, lat0, lat1, lat2, theta, print_exception);

		//Apply Newton-Raphson method
		while ((!convergence) && (iterations <= max_iterations))
		{
			Matrix <T> X(1, 1);
			X(0, 0) = theta;

			//Compute F(theta) and F'(theta)
			const T ftheta = ArithmeticParser::parseEquation(ftheta_equat_postfix, lat, lon, R, a, b, c, lat0, lat1, lat2, theta, print_exception);
			const T ftheta_der = NumDerivative::getDerivative(FProjEquationDerivative1Var <T>(ftheta_equat_postfix, lat, lon, R, a, b, 0, 0, c, lat0, lat1, lat2), X, FirstDerivative, VariableX1, 0.1 * NUM_DERIV_STEP, print_exception);

			//New value of theta
			T theta_n = 0;
			if (fabs(ftheta_der) < MIN_FLOAT)
				theta_n = theta;
			else
				theta_n = theta - ftheta / ftheta_der;

			//Test terminal condition
			if (fabs(theta_n - theta) < max_diff)
			{
				convergence = true;
			}

			//Assign new p
			theta = theta_n;

			iterations++;
		}

		//Convert to deg
		theta = theta * 180 / M_PI;
	}

	//Compute equations including determined theta, if necessary
	const T coord = ArithmeticParser::parseEquation(equation_postfix, lat, lon, R, a, b, c, lat0, lat1, lat2, theta, print_exception) + shift;

	return coord;
}

template <typename T>
void CartTransformation::wgsToJTSK(const Point3DGeographic <T> *p1, Point3DCartesian <T> * p2)
{
	//Deg => rad
	const T PI = 4.0 * atan(1.0), Ro = 57.29577951;

	//WGS-84
	const T A_WGS = 6378137.0000, B_WGS = 6356752.3142;
	const T E2_WGS = (A_WGS * A_WGS - B_WGS * B_WGS) / (A_WGS * A_WGS);

	//Bessel
	const T A_Bes = 6377397.1550, B_Bes = 6356078.9633;
	const T E2_Bes = (A_Bes * A_Bes - B_Bes * B_Bes) / (A_Bes * A_Bes), E_Bes = sqrt(E2_Bes);

	//Scale, Translation, Rotation
	const T m = -3.5623e-6;
	const T X0 = -570.8285, Y0 = -85.6769, Z0 = -462.8420;
	const T OMX = 4.9984 / 3600 / Ro, OMY = 1.5867 / 3600 / Ro, OMZ = 5.2611 / 3600 / Ro;

	//JTSK
	const T  FI0 = 49.5 / Ro;
	const T U0 = (49.0 + 27.0 / 60 + 35.84625 / 3600) / Ro;
	const T ALFA = 1.000597498372;
	const T LA_FERRO = (17.0 + 40.0 / 60) / Ro;
	const T K = 0.9965924869;
	const T R = 6380703.6105;
	const T UK = (59.0 + 42.0 / 60 + 42.6969 / 3600) / Ro;
	const T VK = (42.0 + 31.0 / 60 + 31.41725 / 3600) / Ro;
	const T Ro0 = 1298039.0046;
	const T S0 = 78.5 / Ro;

	//Point parameters
	T W = sqrt(1 - E2_WGS * sin(p1->getLat() / Ro) * sin(p1->getLat() / Ro));
	T M = A_WGS * (1 - E2_WGS) / (W * W * W);
	T N = A_WGS / W;

	//Transformation (B,L,H)WGS => (X,Y,Z)WGS
	T X_WGS = (N + p1->getH()) * cos(p1->getLat() / Ro) * cos(p1->getLon() / Ro);
	T Y_WGS = (N + p1->getH()) * cos(p1->getLat() / Ro) * sin(p1->getLon() / Ro);
	T Z_WGS = (N * (1 - E2_WGS) + p1->getH()) * sin(p1->getLat() / Ro);

	//Transformation (X,Y,Z)WGS =>(X,Y,Z)Bes
	T X_Bes = X0 + (m + 1) * (X_WGS + Y_WGS * OMZ - Z_WGS * OMY);
	T Y_Bes = Y0 + (m + 1) * (-X_WGS * OMZ + Y_WGS + Z_WGS * OMX);
	T Z_Bes = Z0 + (m + 1) * (X_WGS * OMY - Y_WGS * OMX + Z_WGS);

	//Transformation (X,Y,Z)Bes => (BLH)Bess
	T rad = sqrt(X_Bes * X_Bes + Y_Bes * Y_Bes);
	T la = 2 * atan(Y_Bes / (rad + X_Bes));
	T p = atan((A_Bes * Z_Bes) / (B_Bes * rad));
	T t = (Z_Bes + E2_Bes * A_Bes * A_Bes / B_Bes * pow(sin(p), 3)) / (rad - E2_Bes * A_Bes * pow(cos(p), 3));
	T fi = atan(t);
	T H = sqrt(1 + t * t) * (rad - A_Bes / sqrt(1 + (1 - E2_Bes) * t * t));

	//Transformation (fi,la)Bes => (u,v)sphere  (Gauss conformal projection)
	la = la + LA_FERRO;
	T ro = 1 / K * pow(tan(fi / 2 + PI / 4) * pow((1 - E_Bes * sin(fi)) / (1 + E_Bes * sin(fi)), E_Bes / 2), ALFA);
	T u = 2 * atan(ro) - PI / 2;
	T v = ALFA * la;

	//Transformation (u,v)sphere => (s,d)sphere
	T s = asin(sin(UK) * sin(u) + cos(UK) * cos(u) * cos(VK - v));
	T d = asin(sin(VK - v) * cos(u) / cos(s));

	//Transformation (s,d)sphere => (Ro,Eps)plane (Lambert conformal projection)
	T n = sin(S0);
	T Ro_JTSK = Ro0 * pow((tan((S0 / 2 + PI / 4)) / (tan(s / 2 + PI / 4))), n);
	T eps_JTSK = n * d;

	//(Ro, eps) => (x,y)
	T X_JTSK = Ro_JTSK * cos(eps_JTSK);
	T Y_JTSK = Ro_JTSK * sin(eps_JTSK);

	//Set computed parameters to the projected point
	if (p1->getPointLabel() != NULL)
	{
		p2->setPointLabel(p1->getPointLabel());
	}

	p2->setX(X_JTSK);
	p2->setY(Y_JTSK);
	p2->setZ(H);
}

#endif