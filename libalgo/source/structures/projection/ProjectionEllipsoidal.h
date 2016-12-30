// Description: Ellipsoidal projection, derived from Projection

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

#ifndef ProjectionEllipsoidal_H
#define ProjectionEllipsoidal_H


#include "Projection.h"


//Ellipsoidal projection
template <typename T>
class ProjectionEllipsoidal : virtual public Projection <T>
{
        protected:

                T a;
                T b;
                T lat0;

        public:
                ProjectionEllipsoidal() : Projection <T> (), a ( 1 ), b ( 1 ), lat0 ( 0 ) {}
		ProjectionEllipsoidal(const T R_, const T a_, const T b_, const T lat0_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_, 
			const char * projection_family_, const char * projection_name_) : Projection <T>(R_, lon0_, dx_, dy_, c(c_), x_equat_, y_equat_, x_equat_postfix_, 
			y_equat_postfix_, projection_family_, projection_name_), a(a_), b(b_), lat0(lat0_) {}
                virtual ~ProjectionEllipsoidal() {}

        public:

                virtual Point3DGeographic <T> getCartPole() const {return Point3DGeographic <T> ( MAX_LAT, 0.0 );}
                virtual T getLat0() const {return lat0;}
                virtual T getLat1() const {return lat0;}
                virtual T getLat2() const {return lat0;}
                virtual T getA() const {return a;}
                virtual T getB() const {return b;}

                virtual TMinMax <T> getLatPInterval () const {return TMinMax <T> ( MAX_LAT, MAX_LAT );}
                virtual TMinMax <T> getLonPInterval () const {return TMinMax <T> ( 0.0, 0.0 );}
                virtual TMinMax <T> getLat0Interval () const {return TMinMax <T> ( 0.0, MAX_LAT0 );}
                virtual TMinMax <T> getLatPIntervalH ( const TMinMax <T> &lat ) const {return getLatPInterval();}
                virtual TMinMax <T> getLonPIntervalH ( const TMinMax <T> &lon ) const {return getLonPInterval();}
                virtual TTransformedLongtitudeDirection getLonDir () const { return NoDirection;}
		virtual const char * getFThetaEquat() const { return NULL; }
		virtual const char * getTheta0Equat() const { return NULL; }
		virtual const TPostfixNotationDel * getFThetaEquatPostfix() const { return NULL; } 
		virtual const TPostfixNotationDel * getTheta0EquatPostfix() const { return NULL; } 
		virtual TPostfixNotationDel * getFThetaEquatPostfix() { return NULL; }
		virtual TPostfixNotationDel * getTheta0EquatPostfix() { return NULL; }

                virtual void setCartPole ( const Point3DGeographic <T> & pole )  {}
                virtual void setLat0 ( const T lat0_ ) {lat0 = lat0_;}
                virtual void setLat1 ( const T lat1 ) {};
                virtual void setLat2 ( const T lat2 ) {}
                virtual void setA ( const T a_ ) {a = a_;}
                virtual void setB ( const T b_ ) {b = b_;}
                virtual void setLonDir ( const TTransformedLongtitudeDirection lon_dir_ ) {}
		virtual void setFThetaEquat(const char * ftheta_equat_) {};
		virtual void setTheta0Equat(const char * theta0_equat_) {};
		virtual void setFThetaEquatPostfix (const TPostfixNotationDel & ftheta_equat_postfix_) {};
		virtual void setTheta0EquatPostfix (const TPostfixNotationDel & theta0_equat_postfix_) {};
		virtual void FThetaEquatToPostfix() {};
		virtual void Theta0EquatToPostfix() {};

                virtual void getShortCut ( char * shortcut ) const { strcpy ( shortcut, "Elli" ); }
                virtual ProjectionEllipsoidal <T> *clone() const {return new ProjectionEllipsoidal <T> ( *this );}
                virtual void print ( std::ostream * file = &std::cout ) const {}
};

#endif
