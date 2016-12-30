// Description: General cartographic projection with pure wirtual methods
// For evaluations, its equations in the postfix notation are stored
// Supported equations in the non-closed form,

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


#ifndef Projection_H
#define Projection_H

#include <algorithm>
#include <iostream>
#include <ostream>
#include <string.h>


#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"
#include "libalgo/source/algorithms/arithmeticparser/ArithmeticParser.h"

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

                char * x_equat;					//Equation X
				char * y_equat;					//Equation Y
				TPostfixNotationDel x_equat_postfix;		//Equation X converted to the postfix notation
				TPostfixNotationDel y_equat_postfix;		//Equation Y converted to the postfix notation
				char * projection_family;			//Projection family
				char * projection_name;				//Projection name

		public:
			Projection() : R(1.0), lon0(0.0), dx(0.0), dy(0.0), c(0.5), x_equat(NULL), y_equat(NULL), x_equat_postfix(0), y_equat_postfix(0), projection_family(NULL), projection_name(NULL) {}
			Projection(const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_);
			Projection(const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_, const char * projection_family_, const char * projection_name_);
			Projection(const Projection <T> *proj);
			Projection(const Projection <T> &proj);
			virtual ~Projection() = 0;

        public:
                //Methods common for all derived classes (methods are not virtual)
                T getR() const  {return R;}
                T getLon0() const  {return lon0;}
                T getDx() const  {return dx;}
                T getDy() const  {return dy;}
                T getC() const  {return c;}


				const char * getXEquat() const { return x_equat; }
				const char * getYEquat() const { return y_equat; }
				const TPostfixNotationDel * getXEquatPostfix() const { return &x_equat_postfix; }
				const TPostfixNotationDel * getYEquatPostfix() const { return &y_equat_postfix; }
				TPostfixNotationDel * getXEquatPostfix() { return &x_equat_postfix; }
				TPostfixNotationDel * getYEquatPostfix() { return &y_equat_postfix; }
				const char * getFamily() const { return projection_family; }
				const char * getName() const { return projection_name; }

				void setR(const T R_) { R = R_; }
				void setLon0(const T lon0_) { lon0 = lon0_; }
				void setDx(const T dx_) { dx = dx_; }
				void setDy(const T dy_) { dy = dy_; }
				void setC(const T c_) { c = c_; }
				void XEquatToPostfix();
				void YEquatToPostfix();

				void setXEquat(const char * x_equat_);
				void setYEquat(const char * y_equat_);
				void setXEquatPostfix(const TPostfixNotationDel &x_equat_postfix_) { x_equat_postfix = x_equat_postfix_; }
				void setYEquatPostfix(const TPostfixNotationDel &y_equat_postfix_) { y_equat_postfix = y_equat_postfix_; }
				void setFamily(const char * projection_family_);
				void setName(const char * projection_name_);


			public:
				//Methods different for all derived classes ( methods are virtual )
				virtual Point3DGeographic <T> getCartPole() const = 0;
				virtual T getLat0() const = 0;
				virtual T getLat1() const = 0;
				virtual T getLat2() const = 0;
				virtual T getA() const = 0;
				virtual T getB() const = 0;
				virtual TMinMax <T> getLatPInterval() const = 0;
				virtual TMinMax <T> getLonPInterval() const = 0;
				virtual TMinMax <T> getLat0Interval() const = 0;
				virtual TMinMax <T> getLatPIntervalH(const TMinMax <T> &lat) const = 0;
				virtual TMinMax <T> getLonPIntervalH(const TMinMax <T> &lon) const = 0;
				virtual TTransformedLongitudeDirection getLonDir() const = 0;
				virtual const char * getFThetaEquat() const = 0;
				virtual const char * getTheta0Equat() const = 0;
				virtual const TPostfixNotationDel * getFThetaEquatPostfix() const = 0;
				virtual const TPostfixNotationDel * getTheta0EquatPostfix() const = 0;
				virtual TPostfixNotationDel * getFThetaEquatPostfix() = 0;
				virtual TPostfixNotationDel * getTheta0EquatPostfix() = 0;

				virtual void setCartPole(const Point3DGeographic <T> & pole) = 0;
				virtual void setLat0(const T lat0) = 0;
				virtual void setLat1(const T lat1) = 0;
				virtual void setLat2(const T lat2) = 0;
				virtual void setA(const T a) = 0;
				virtual void setB(const T b) = 0;
				virtual void setLonDir(const TTransformedLongitudeDirection lon_dir_) = 0;
				virtual void setFThetaEquat(const char * ftheta_equat_) = 0;
				virtual void setTheta0Equat(const char * theta0_equat_) = 0;
				virtual void setFThetaEquatPostfix(const TPostfixNotationDel & ftheta_equat_postfix_) = 0;
				virtual void setTheta0EquatPostfix(const TPostfixNotationDel & theta0_equat_postfix_) = 0;
				virtual void FThetaEquatToPostfix() = 0;
				virtual void Theta0EquatToPostfix() = 0;

				virtual void getShortCut(char * shortcut) const = 0;
				virtual Projection <T> *clone() const = 0;
				virtual void print(std::ostream * file = &std::cout) const = 0;
			};

#include "Projection.hpp"

#endif
