// Description: Ellipsoidal projection, derived from Projection

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


#ifndef ProjectionEllipsoidal_H
#define ProjectionEllipsoidal_H

#include "Projection.h"

#include "libalgo/source/const2/Const.h"


//Ellipsoidal projection
template <typename T>
class ProjectionEllipsoidal : public Projection <T>
{
protected:

	T a, b;				//Ellpsoid semi-axes
	T lat1;				//True parallel

public:
	ProjectionEllipsoidal() : Projection <T>(), a(1), b(1), lat1(0) {}
	ProjectionEllipsoidal(const T R_, const T a_, const T b_, const T lat1_, const T lon0_, const T dx_, const T dy_, const T c_, const p_coord_function <T> &pX, const p_coord_function <T> &pY,
		const std::string name_) : Projection <T>(R_, lon0_, dx_, dy_, c_, pX, pY, name_), a(a_), b(b_), lat1(lat1_) {}
	virtual ~ProjectionEllipsoidal() {}

public:

	virtual Point3DGeographic <T> getCartPole() const { return Point3DGeographic <T>(MAX_LAT, 0.0); }
	virtual T getLat1() const { return lat1; }
	virtual T getLat2() const { return lat1; }
	virtual T getA() const { return a; }
	virtual T getB() const { return b; }
	virtual TTransformedLongitudeDirection getLonDir() const { return NoDirection; }
	virtual std::string getFamily() const { return "Elips"; }
	virtual TInterval <T> getLatPInterval() const { return TInterval <T> {MAX_LAT, MAX_LAT}; }
	virtual TInterval <T> getLonPInterval() const { return TInterval <T> {0.0, 0.0}; }
	virtual TInterval <T> getLat1Interval() const { return TInterval <T> {0, MAX_LAT1}; }
	
	virtual void setCartPole(const Point3DGeographic <T> & pole) {}
	virtual void setLat1(const T lat1_) { lat1 = lat1_; };
	virtual void setLat2(const T lat2) {}
	virtual void setA(const T a_) { a = a_; }
	virtual void setB(const T b_) { b = b_; }
	virtual void setLonDir(const TTransformedLongitudeDirection &lon_dir_) {}
};

#endif
