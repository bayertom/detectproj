// Description: Functor, compute matrix cost of results for cartometric analysis

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


#ifndef FAnalyzeProjC_H
#define FAnalyzeProjC_H


#include "libalgo/source/algorithms/matrixoperations/MatrixOperations.h"


using namespace MatrixOperations;


template <typename T>
class FAnalyzeProjC
{
        public:

                FAnalyzeProjC () {}

                T operator () ( const Matrix <T> &X, Matrix <T> &C )
                {
                        //Compute cost of the solution using the standard deviation
                        return sqrt ( sum2 ( C ) / ( 0.5 * C.rows() - 2 ) );
                }

};


#endif
