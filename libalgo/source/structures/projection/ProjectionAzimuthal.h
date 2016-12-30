// Description: Azimuthal projection, derived from Projection

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


#ifndef ProjectionAzimuthal_H
#define ProjectionAzimuthal_H


#include "Projection.h"

#include "libalgo/source/structures/point/Point3DGeographic.h"


//Azimuthal projection
template <typename T>
class ProjectionAzimuthal : virtual public Projection <T>
{
        protected:
                Point3DGeographic <T> cart_pole;
                TTransformedLongtitudeDirection lon_dir;

        public:
                ProjectionAzimuthal() : Projection <T> (), cart_pole ( MAX_LAT, 0.0 ), lon_dir ( NormalDirection ) {}
		ProjectionAzimuthal(const T R_, const T latp_, const T lonp_, const TTransformedLongtitudeDirection lon_dir_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_, 
			const char * projection_family_, const char * projection_name_) : Projection <T>(R_, lon0_, dx_, dy_, c_, x_equat_, y_equat_, x_equat_postfix_, y_equat_postfix_, projection_family_, 
			projection_name_), cart_pole(latp_, lonp_), lon_dir(lon_dir_) {}
                virtual ~ProjectionAzimuthal() {}

        public:
                virtual Point3DGeographic <T> getCartPole() const {return cart_pole; }
                virtual T getLat0() const {return 0.0;}
                virtual T getLat1() const {return 0.0;}
                virtual T getLat2() const {return 0.0;}
                virtual T getA() const {return this->R;}
                virtual T getB() const {return this->R;}

                virtual TMinMax <T> getLatPInterval () const {return TMinMax <T> ( MIN_LAT, MAX_LAT );}
                virtual TMinMax <T> getLonPInterval () const {return TMinMax <T> ( MIN_LON, MAX_LON );}
                virtual TMinMax <T> getLat0Interval () const {return TMinMax <T> ( 0, 0 );}
                virtual TMinMax <T> getLatPIntervalH ( const TMinMax <T> &lat ) const { return ( lat.max_val <= 0 ? TMinMax <T> ( MIN_LAT, 0.0 ) : ( lat.max_val >= 0 ? TMinMax <T> ( 0.0, MAX_LAT ) : getLatPInterval() ) );}
                virtual TMinMax <T> getLonPIntervalH ( const TMinMax <T> &lon ) const { return lon; }
                virtual TTransformedLongtitudeDirection getLonDir () const { return lon_dir; }
		virtual const char * getFThetaEquat() const { return NULL; }
		virtual const char * getTheta0Equat() const { return NULL; }
		virtual const TPostfixNotationDel * getFThetaEquatPostfix() const { return NULL; }
		virtual const TPostfixNotationDel * getTheta0EquatPostfix() const { return NULL; }
		virtual TPostfixNotationDel * getFThetaEquatPostfix() { return NULL;}
		virtual TPostfixNotationDel * getTheta0EquatPostfix() { return NULL;}

                virtual void setCartPole ( const Point3DGeographic <T> & cart_pole_ ) { cart_pole = cart_pole_; }
                virtual void setLat0 ( const T lat0 ) {};
                virtual void setLat1 ( const T lat1 ) {};
                virtual void setLat2 ( const T lat2 ) {};
                virtual void setA ( const T a ) {};
                virtual void setB ( const T b ) {};
                virtual void setLonDir ( const TTransformedLongtitudeDirection lon_dir_ ) { lon_dir = lon_dir_; }
		virtual void setFThetaEquat(const char * ftheta_equat_) {};
		virtual void setTheta0Equat(const char * theta0_equat_) {};
		virtual void setFThetaEquatPostfix (const TPostfixNotationDel & ftheta_equat_postfix_) {};
		virtual void setTheta0EquatPostfix (const TPostfixNotationDel & theta0_equat_postfix_) {};
		virtual void FThetaEquatToPostfix() {};
		virtual void Theta0EquatToPostfix() {};

                virtual void getShortCut ( char * shortcut ) const { strcpy ( shortcut, "Azim" ); }
                virtual ProjectionAzimuthal <T> *clone() const {return new ProjectionAzimuthal ( *this );}
                virtual void print ( std::ostream * file = &std::cout ) const {}
};

#endif

