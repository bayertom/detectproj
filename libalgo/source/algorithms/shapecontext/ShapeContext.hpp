// Description: Compute Shape Context used for Inner Distance Analysis

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


#ifndef ShapeContext_HPP
#define ShapeContext_HPP


#include "libalgo/source/exceptions/MathMatrixNotSquareException.h"


template <typename T>
void ShapeContext::computeShapeContext ( const Matrix <T> &D, const Matrix <unsigned short> &TT, Matrix <unsigned short> &SC, unsigned int point_index )
{
        //Compute shape context of the point (index of the matrix) from given distance and theta matrices of the Face
        const unsigned int m = D.rows();

        //Compute shape context for matrices D (length) and T (theta)
        for ( unsigned int j = 0; j < m; j++ )
        {
                //Get row and col index: use only valid non-zero items
                if ( ( D ( point_index, j ) > 0 ) && ( D ( point_index, j ) <= 1 ) && ( TT ( point_index, j ) > 0 ) )
                {
                        unsigned short row = ( unsigned short ) ( D ( point_index, j ) / 0.2 );
                        unsigned short col = ( TT ( point_index, j ) / 30 );

                        //Throw exception
                        if ( row > 5 )
                        {
                                throw IndexOutOfBoundException ( "IndexOutOfBoundException, ", "can not compute the shape context: invalid dimension of the matrix (row > rows_count)" );
                        }

                        //Throw exception
                        if ( col > 11 )
                        {
                                throw IndexOutOfBoundException ( "IndexOutOfBoundException, ", "can not compute the shape context: invalid dimension of the matrix (row > rows_count)" );
                        }

                        //Increment value
                        SC ( row, col ) += 1;
                }
        }
}


template <typename T>
T ShapeContext::compare2ShapeContexts ( const Matrix <unsigned short> &SC1, const Matrix <unsigned short> &SC2 )
{
        //Compare 2 shape contexts using shape context (X - quadrat comparision)
        T sum = 0;

        for ( unsigned int i = 0; i < 5; i++ )
        {
                for ( unsigned int j = 0; j < 12; j++ )
                {
                        if ( ( SC1 ( i, j ) + SC2 ( i, j ) ) != 0 )
                        {
                                sum += ( SC1 ( i, j ) - SC2 ( i, j ) ) * ( SC1 ( i, j ) - SC2 ( i, j ) ) / ( T ) ( SC1 ( i, j ) + SC2 ( i, j ) );
                        }
                }
        }

        return  0.5 * sum;
}

#endif
