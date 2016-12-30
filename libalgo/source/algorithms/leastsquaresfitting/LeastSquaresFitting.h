// Description: Least squares fitting, 2D line

// Copyright (c) 2010 - 2011
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


#ifndef LeastSquaresFitting_H
#define LeastSquaresFitting_H

#include <set>
#include <ostream>
#include <iostream>

#include "libalgo/source/structures/list/Container.h"

//Forward declaration
template <typename T>
struct TIndexSet;


//New user type definition: parameters of the regression line
template <typename Point>
struct TFittingLine
{
        unsigned int hash_val;                                                  //Hash of b and alpha

        //Fit line equation: c * y = tan( alpha ) * x + b
        typename Point::Type alpha,						//Angle of the regression line with the horizontal axis x
                 b,								//Shift of the regression line
                 c,								//Coefficient, c = 0 for vertical line
                 xt, yt,							//Centroids of the set
                 error;								//Fit error

        typename TIndexSet <Point>::Type points_indices;			//Set of indices of the regression line points sorted by y-coordinate and subsequently by x-coordinate

        TFittingLine ( const Container <Point> &nl ) : hash_val ( 0.0 ),
                alpha ( 0.0 ), b ( 0.0 ), c ( 0.0 ), xt ( 0.0 ),
                yt ( 0.0 ), error ( 0.0 ), points_indices ( nl ) {}

        void printIndices()
        {
                for ( typename TIndexSet <Point>::Type::iterator i_set = points_indices.begin(); i_set != points_indices.end(); i_set++ )
                        std::cout << *i_set << " "; std::cout << '\n';
        }
};

//Find best fitting curve using Least squares
class LeastSquaresFitting
{
        public:
                template <typename Point>
                static void fitLine ( const Container <Point> &nl, TFittingLine <Point> &regression_line, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );
};


#include "LeastSquaresFitting.hpp"

#endif
