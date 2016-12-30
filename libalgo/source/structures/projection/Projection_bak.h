// Description: Cartographic projection

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


#ifndef Projection_H
#define Projection_H

#include <algorithm>
#include <iostream>
#include <ostream>
#include <string.h>


#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"

#include "libalgo/source/io/File.h"


//Structure stores min-max values
template <typename T>
struct TMinMax
{
        T min_val, max_val;
        typedef T value_type;

        TMinMax <T> ( ) : min_val ( 0.0 ), max_val ( 0.0 ) {}
        TMinMax <T> ( const T min_val_, const T max_val_ ) : min_val ( min_val_ ), max_val ( max_val_ ) {}
};


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
		T alpha;					//Rotatiton of the projection
                char * x_equat;					//Equation x
                char * y_equat;					//Equation y
                char * projection_name;				//Projection name

		////
		T x_mass_test;
		T y_mass_test;
		T x_mass_reference;
		T y_mass_reference;

        public:
                Projection() : R ( 1.0 ), lon0 ( 0.0 ), dx ( 0.0 ),  dy ( 0.0 ), c ( 1.0 ), alpha ( 0.0 ), x_equat ( NULL ), y_equat ( NULL ), projection_name ( NULL ) {}
                Projection ( const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_ );
                Projection ( const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const char * projection_name_ );
                Projection ( const Projection <T> *proj );
                Projection ( const Projection <T> &proj );
                virtual ~Projection() = 0;

        public:
                //Methods common for all derived classes (methods are not virtual)
                T getR() const  {return R;}
                T getLon0() const  {return lon0;}
                T getDx() const  {return dx;}
                T getDy() const  {return dy;}
                T getC() const  {return c;}
		T getAlpha() const {return alpha;}
		
		///////////
		T getXMassTest() const {return x_mass_test;}
		T getYMassTest() const {return y_mass_test;}
		T getXMassReference() const {return x_mass_reference;}
		T getYMassReference() const {return y_mass_reference;}
		//////////

                const char * getXEquat () const {return x_equat;}
                const char * getYEquat () const {return y_equat;}
                const char * getProjectionName () const {return projection_name;}

                void setR ( const T R_ ) {R = R_;}
                void setLon0 ( const T lon0_ ) {lon0 = lon0_;}
                void setDx ( const T dx_ ) {dx = dx_;}
                void setDy ( const T dy_ ) {dy = dy_;}
                void setC ( const T c_ ) {c = c_;}
		void setAlpha( const T alpha_) { alpha = alpha_;}
                void setXEquat ( const char * x_equat_ );
                void setYEquat ( const char * y_equat_ );
                void setProjectionName ( const char * projection_name_ );

		////////////////
		void setXMassTest ( const T x_mass_test_ ) {x_mass_test = x_mass_test_;}
		void setYMassTest ( const T y_mass_test_ ) {y_mass_test = y_mass_test_;}
		void setXMassReference ( const T x_mass_reference_ ) {x_mass_reference = x_mass_reference_;}
		void setYMassReference ( const T y_mass_reference_ ) {y_mass_reference = y_mass_reference_;}
		///////////////


        public:
                //Methods different for all derived classes ( methods are virtual )
                virtual Point3DGeographic <T> getCartPole() const = 0;
                virtual T getLat0() const = 0;
                virtual T getLat1() const = 0;
                virtual T getLat2() const = 0;
                virtual T getA() const = 0;
                virtual T getB() const = 0;
                virtual TMinMax <T> getLatPInterval () const = 0;
                virtual TMinMax <T> getLonPInterval () const = 0;
                virtual TMinMax <T> getLat0Interval () const = 0;
                virtual TMinMax <T> getLatPIntervalH ( const TMinMax <T> &lat ) const = 0;
                virtual TMinMax <T> getLonPIntervalH ( const TMinMax <T> &lon ) const = 0;
                virtual TTransformedLongtitudeDirection getLonDir () const = 0;

                virtual void setCartPole ( const Point3DGeographic <T> & pole )  = 0;
                virtual void setLat0 ( const T lat0 ) = 0;
                virtual void setLat1 ( const T lat1 ) = 0;
                virtual void setLat2 ( const T lat2 ) = 0;
                virtual void setA ( const T a ) = 0;
                virtual void setB ( const T b ) = 0;
                virtual void setLonDir ( const TTransformedLongtitudeDirection lon_dir_ ) = 0;

                virtual void getShortCut ( char * shortcut ) const = 0;
                virtual Projection <T> *clone() const = 0;
                virtual void print ( std::ostream * file = &std::cout ) const = 0;
};

#include "Projection.hpp"

#endif
