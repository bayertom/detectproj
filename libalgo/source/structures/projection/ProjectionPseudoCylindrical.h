// Description: Pseudocylindrical projection, derived from Projection
// Supported equations in the non-closed form
// For evaluation, stored in the postfix notation


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

#ifndef ProjectionPseudoCylindrical_H
#define ProjectionPseudoCylindrical_H


#include "Projection.h"


//Pseudocylindrical projection
template <typename T>
class ProjectionPseudoCylindrical : virtual public Projection <T>
{
protected:

	T lat0;

	Point3DGeographic <T> cart_pole;
	TTransformedLongtitudeDirection lon_dir;

	//If projection equations in the non-closed form contain additional parameters, they are solved by Newton-Raphson method.
	//Functions f(theta), and theta0 are stored
	//New solution theta[i+1] = theta[i] - f(theta[i])/f'(theta[i]) 
	char * ftheta_equat;				//Parameter function of the variable theta
	char * theta0_equat;				//Initial value of parameter theta, theta = theta0
	TPostfixNotationDel ftheta_equat_postfix;	//Parameter function of the variable theta in the postfix notation
	TPostfixNotationDel theta0_equat_postfix;	//Initial value of parameter theta, theta = theta0 in the postfix notation


public:
	ProjectionPseudoCylindrical() : Projection<T>(), lat0(45), cart_pole(MAX_LAT, 0.0), lon_dir(NormalDirection2), ftheta_equat(NULL), theta0_equat(NULL), ftheta_equat_postfix(0), theta0_equat_postfix(0) {}
	ProjectionPseudoCylindrical(const T R_, const T lat0_, const T latp_, const T lonp_, const TTransformedLongtitudeDirection lon_dir_, const char * ftheta_equat_, const char * theta0_equat_, const TPostfixNotationDel * ftheta_equat_postfix_, const TPostfixNotationDel * theta0_equat_postfix_,
		const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_,
		const char * projection_family_, const char * projection_name_);
	ProjectionPseudoCylindrical(const ProjectionPseudoCylindrical <T> &proj);
	virtual ~ProjectionPseudoCylindrical();

public:

	virtual Point3DGeographic <T> getCartPole() const { return cart_pole; }
	virtual T getLat0() const { return lat0; }
	virtual T getLat1() const { return 0.0; }
	virtual T getLat2() const { return 0.0; }
	virtual T getA() const { return this->R; }
	virtual T getB() const { return this->R; }

	//virtual TMinMax <T> getLatPInterval () const {return TMinMax <T> ( MAX_LAT, MAX_LAT );}
	//virtual TMinMax <T> getLonPInterval () const {return TMinMax <T> ( 0.0, 0.0 );}
	virtual TMinMax <T> getLatPInterval() const { return TMinMax <T>(MIN_LAT, MAX_LAT); }
	virtual TMinMax <T> getLonPInterval() const { return TMinMax <T>(MIN_LON, MAX_LON); }
	virtual TMinMax <T> getLat0Interval() const { return TMinMax <T>(0.0, MAX_LAT0); }
	virtual TMinMax <T> getLatPIntervalH(const TMinMax <T> &lat) const { return getLatPInterval(); }
	virtual TMinMax <T> getLonPIntervalH(const TMinMax <T> &lon) const { return getLonPInterval(); }
	virtual TTransformedLongtitudeDirection getLonDir() const { return lon_dir; }
	virtual const char * getFThetaEquat() const { return ftheta_equat; }
	virtual const char * getTheta0Equat() const { return theta0_equat; }
	virtual const TPostfixNotationDel * getFThetaEquatPostfix() const { return ftheta_equat_postfix.size() > 0 ? & ftheta_equat_postfix : NULL; }
	virtual const TPostfixNotationDel * getTheta0EquatPostfix() const { return theta0_equat_postfix.size() > 0 ? &theta0_equat_postfix : NULL; }
	virtual TPostfixNotationDel * getFThetaEquatPostfix() { return ftheta_equat_postfix.size() > 0 ? &ftheta_equat_postfix : NULL; }
	virtual TPostfixNotationDel * getTheta0EquatPostfix() { return theta0_equat_postfix.size() > 0 ? &theta0_equat_postfix : NULL; }

	virtual void setCartPole(const Point3DGeographic <T> & cart_pole_) { cart_pole = cart_pole_; }
	virtual void setLat0(const T lat0_) { lat0 = lat0_; }
	virtual void setLat1(const T lat1) {}
	virtual void setLat2(const T lat2) {}
	virtual void setA(const T a) {}
	virtual void setB(const T b) {}
	virtual void setLonDir(const TTransformedLongtitudeDirection lon_dir_) { lon_dir = lon_dir_; }
	virtual void setFThetaEquat(const char * ftheta_equat_);
	virtual void setTheta0Equat(const char * theta0_equat_);
	virtual void setFThetaEquatPostfix (const TPostfixNotationDel & ftheta_equat_postfix_) { ftheta_equat_postfix = ftheta_equat_postfix_; }
	virtual void setTheta0EquatPostfix (const TPostfixNotationDel & theta0_equat_postfix_) { theta0_equat_postfix = theta0_equat_postfix_; }
	virtual void FThetaEquatToPostfix();
	virtual void Theta0EquatToPostfix();
	virtual void getShortCut(char * shortcut) const { strcpy(shortcut, "PsCyli"); }
	virtual ProjectionPseudoCylindrical <T> *clone() const { return new ProjectionPseudoCylindrical <T>(*this); }
	virtual void print(std::ostream * file = &std::cout) const {}

};

#include "ProjectionPseudoCylindrical.hpp"

#endif
