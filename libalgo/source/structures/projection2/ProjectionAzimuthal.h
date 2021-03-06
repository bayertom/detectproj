// Description: Azimuthal projection, derived from Projection

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


#ifndef ProjectionAzimuthal_H
#define ProjectionAzimuthal_H

#include "Projection.h"

#include "libalgo/source/const2/Const.h"


//Azimuthal projection
template <typename T>
class ProjectionAzimuthal : public Projection <T>
{
protected:
	Point3DGeographic <T> cart_pole;						//Meta-pole position
	TTransformedLongitudeDirection lon_dir;					//Transformed londitude mode

public:
	ProjectionAzimuthal() : Projection <T>(), cart_pole(MAX_LAT, 0.0), lon_dir(NormalDirection2) {}
	ProjectionAzimuthal(const T R_, const T latp_, const T lonp_, const TTransformedLongitudeDirection lon_dir_, const T lon0_, const T dx_, const T dy_, const T c_, const p_coord_function <T> &pX, const p_coord_function <T> &pY,
		const std::string name_) : Projection <T>(R_, lon0_, dx_, dy_, c_, pX, pY, name_), cart_pole(latp_, lonp_), lon_dir(lon_dir_) {}
	virtual ~ProjectionAzimuthal() {}

public:
	virtual Point3DGeographic <T> getCartPole() const { return cart_pole; }
	virtual T getLat1() const { return 0.0; }
	virtual T getLat2() const { return 0.0; }
	virtual T getA() const { return this->R; }
	virtual T getB() const { return this->R; }
	virtual TTransformedLongitudeDirection getLonDir() const { return lon_dir; }
	virtual std::string getFamily() const { return "Azim"; }
	virtual TInterval <T> getLatPInterval() const { return TInterval <T> {MIN_LAT, MAX_LAT}; }
	virtual TInterval <T> getLonPInterval() const { return TInterval <T> {MIN_LON, MAX_LON}; }
	virtual TInterval <T> getLat1Interval() const { return TInterval <T> {0, 0}; }
	
	virtual void setCartPole(const Point3DGeographic <T> & cart_pole_) { cart_pole = cart_pole_; }
	virtual void setLat1(const T lat1) {};
	virtual void setLat2(const T lat2) {};
	virtual void setA(const T a) {};
	virtual void setB(const T b) {};
	virtual void setLonDir(const TTransformedLongitudeDirection &lon_dir_) { lon_dir = lon_dir_; }
};

#endif
