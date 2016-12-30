// Description: RANSAC algorithm implementation

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

#ifndef Ransac_H
#define Ransac_H

#include <set>
#include <stdio.h>


#include "libalgo/source/structures/list/Container.h"
#include "libalgo/source/algorithms/leastsquaresfitting/LeastSquaresFitting.h"


//Forward declarations
template <typename T>
class Node3DCartesianProjected;

template <typename T>
class Point3DGeographic;

template <typename T>
struct TFittingLine;

template <typename T>
class sortFittingLinesByHash;

template <typename T>
struct TMeridiansList;

template <typename T>
struct TParallelsList;


//New user type
template <typename Point>
struct TRansacResults
{
        typedef std::set < TFittingLine <Point>, sortFittingLinesByHash <TFittingLine <Point> > > Type;
};


//Fit curves through the dataset using RANSAC algorithm
class Ransac
{
        public:

                template <typename Point>
                static bool ransacFitLine ( const Container <Point> &input, TFittingLine <Point> &acceptable_solution, const typename Point::Type acceptable_error, const bool find_best, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point>
                static void ransacFitAllLines ( const Container <Point> &input, typename TRansacResults <Point> ::Type & ransac_results, const typename Point::Type acceptable_error, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static void ransacFitMeridiansAndParallels ( const Container <Point3DGeographic <T> *> &pl_geographic,  typename TMeridiansList <T> ::Type &meridians,  typename TParallelsList <T> ::Type &parallels, const T acceptable_error = 0.5, const T angle_tolerance = 1.0,
                                const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );
};

#include "Ransac.hpp"

#endif
