// Description: General cartographic projection with pure wirtual methods
// Supported equations in the non-closed form,

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

#ifndef Projection_H
#define Projection_H

#include <cstring>

#include "libalgo/source/types/TInterval.h"
#include "libalgo/source/types/TTransformedLongitudeDirection.h"
#include "libalgo/source/types/PCoordFunction.h"

#include "libalgo/source/structures/point/Point3DGeographic.h"

//Common cartographic projection: abstract class
template <typename T>
class Projection 
{
	protected:

		T R;						//Sphere radius
		T lon0;						//Projection central parallel
		T dx;						//Aditive constant dx
		T dy;						//Aditive constant dy
		T c;						//Additional constant of the projection
		std::string name;				//Projection name (In accordance with Proj.4)

		p_coord_function <T> X;				//Pointer to the coordinate function X = F(lat, lon)
		p_coord_function <T> Y;				//Pointer to the coordinate function Y = Y(lat, lon)

	public:
		Projection();					//The remaining part of the constructor is in Projections.h

		Projection(const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const p_coord_function <T> &pX, const p_coord_function <T> &pY, const std::string name_) :
			R(R_), lon0(lon0_), dx(dx_), dy(dy_), c(c_), name(name_), X(pX), Y(pY) 	{};

		virtual ~Projection() = 0;

		//Methods common for all derived classes (methods are not virtual)
		T getR() const	{return R;}
		T getLon0() const {return lon0;}
		T getDx() const	{return dx;}
		T getDy() const {return dy;}
		T getC() const	{return c;}
		T getX(const T lat, const T lon) const { return X(this->getR(), this->getLat1(), this->getLat2(), lat, lon, this->getLon0(), this->getDx(), this->getDy(), this->getC()); }
		T getX(const T R_, const T lat1_, const T lat2_, const T lat_, const T lon_, const T lon0_, const T dx_, const T dy_, const T c_) {return X(R_, lat1_, lat2_, lat_, lon_, lon0_, dx_, dy_, c_);}
		const p_coord_function <T> getX() const { return X; }
		T getY(const T lat, const T lon) const { return Y(this->getR(), this->getLat1(), this->getLat2(), lat, lon, this->getLon0(), this->getDx(), this->getDy(), this->getC()); }
		T getY(const T R_, const T lat1_, const T lat2_, const T lat_, const T lon_, const T lon0_, const T dx_, const T dy_, const T c_) { return Y(R_, lat1_, lat2_, lat_, lon_, lon0_, dx_, dy_, c_); }
		const p_coord_function <T> getY() const { return Y; }
		std::string getName() const {return name;}

		void setR(const T R_)	{R = R_;}
		void setLon0(const T lon0_)	{lon0 = lon0_;}
		void setDx(const T dx_)	{dx = dx_;}
		void setDy(const T dy_)	{dy = dy_;}
		void setC(const T c_) {c = c_;}
		void setName(const std::string &name_){name = name_;}

	public:
		//Methods different for all derived classes ( methods are virtual )
		virtual Point3DGeographic <T> getCartPole() const = 0;
		virtual T getLat1() const = 0;
		virtual T getLat2() const = 0;
		virtual T getA() const = 0;
		virtual T getB() const = 0;
		virtual TTransformedLongitudeDirection getLonDir() const = 0;
		virtual std::string getFamily() const = 0;
		virtual TInterval <T> getLatPInterval() const = 0;
		virtual TInterval <T> getLonPInterval() const = 0;
		virtual TInterval <T> getLat1Interval() const = 0;

		virtual void setCartPole(const Point3DGeographic <T> & pole) = 0;
		virtual void setLat1(const T lat1) = 0;
		virtual void setLat2(const T lat2) = 0;
		virtual void setA(const T a) = 0;
		virtual void setB(const T b) = 0;
		virtual void setLonDir(const TTransformedLongitudeDirection &lon_dir_) = 0;
};

//Implementation of the pure virtual destructor
template <typename T>
Projection <T>::~Projection()
{
};


#endif
