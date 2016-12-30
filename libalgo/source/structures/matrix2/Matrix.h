// Description: Matrix and basic operator

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


#ifndef Matrix_H
#define Matrix_H

#include <vector>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <array>

#include "libalgo/source/types/TVector2D.h"


//Matrix definition and basic operators
template <typename T>
class Matrix
{
        private:
                TVector2D <T> items;				//Matrix items
                const unsigned int rows_count;			//Total rows_count
                const unsigned int columns_count;		//Total columns_count

        public:

                Matrix ( const unsigned int rows_count_, const unsigned int columns_count_, const T item = 0 ) : rows_count ( rows_count_ ), columns_count ( columns_count_ ), items ( rows_count_, std::vector <T> ( columns_count_, item ) ) {}
                Matrix ( const unsigned int rows_count_, const unsigned int columns_count_, const T non_diag_item_val, const T diag_val );
		
		template <std::size_t rows_count_, std::size_t columns_count_>
		Matrix (const std::array<std::array<T, columns_count_>, rows_count_> &data);

                template <typename U>
                Matrix ( const Matrix <U> &M ) ;

                ~Matrix();

        public:
	
		//Get rows and columns count
                inline unsigned int rows() const {return rows_count;}
		inline unsigned int cols() const {return columns_count;}

                //Get Matrix content
		TVector2D <T> const & getItems () const {return items;}
		TVector2D <T> & getItems ()  {return items;}

                //Get row, col
                Matrix <T> row ( const unsigned int r ) const;
                Matrix <T> col ( const unsigned int c ) const;

                //Set row, col, submatrix
                void row ( const Matrix <T> & R, const unsigned int r );
                void col ( const Matrix <T> & C, const unsigned int c );
                void resize ( const int r, const int c, const T val = 0.0 );

                //Other methods
                void print ( std::ostream * output = &std::cout ) const;

        public:

                //Matrix operators =, ==, !=
                Matrix <T> & operator = ( const Matrix <T> &M );
                bool operator == ( const Matrix <T> &M ) const;
                bool operator != ( const Matrix <T> &M ) const;

                //Matrix operators + : Matrix + Matrix
                template <typename U>
                Matrix <T> operator + ( Matrix <U> const &M ) const;

		//Matrix operators + : Matrix + scalar
		template <typename U>
		Matrix <T> operator + (const U val) const;

                //Matrix operators - : Matrix - Matrix
                template <typename U>
                Matrix <T> operator - ( Matrix <U> const &M ) const;

		//Matrix operators + : Matrix - scalar
		template <typename U>
		Matrix <T> operator - (const U val) const;

                //Matrix operators += : Matrix += Matrix
                template <typename U>
                Matrix <T> & operator += ( Matrix <U> const &M );

                //Matrix operators -= : Matrix -= Matrix
                template <typename U>
                Matrix <T> & operator -= ( Matrix <U> const &M );

                //Matrix operators * : Matrix * Matrix
                template <typename U>
                Matrix <T> operator * ( Matrix <U> const &M ) const;

                //Matrix operators * : Matrix * Scalar
                template <typename U>
                Matrix <T> operator * ( const U & val ) const;

                //Matrix operators *= : Matrix *= Scalar
                template <typename U>
                Matrix <T> & operator *= ( const U & val );

                //Matrix operators / : Matrix / Scalar
                template <typename U>
                Matrix <T> operator / ( const U & val ) const;

                //Matrix operators /= : Matrix /= Scalar
                template <typename U>
                Matrix <T> & operator /= ( const U & val );

                //Matrix operators /= : Matrix /= Matrix
                template <typename U>
                Matrix <T> & operator /= ( Matrix <U> const &M );

                //Matrix operators .* : Matrix .* Matrix (Hadamard product)
                template <typename U>
                Matrix <T> operator % ( Matrix <U> const &M ) const;

                //Matrix operator (row, col)
                inline T & operator() ( const unsigned int row, const unsigned int col );

                inline T const & operator() ( const unsigned int row, const unsigned int col ) const;

                //Matrix operator (r1, r2, c1, c2): get submatrix of the matrix
                inline Matrix <T> operator () ( const unsigned int r1, const unsigned int r2,
                                          const unsigned int c1, const unsigned int c2 ) const;

		//Matrix operator (M, r, c): replace part of the matrix A with M
		void operator () (const Matrix <T> &M, const unsigned int row, const unsigned int col);

                //Overriden operator <<
                friend void operator << ( std::ostream & output, const Matrix <T> &m ) { m.print ( &output ); }
};

//Multiply matrix with a scalar, non-member function
template <typename T>
Matrix<T> operator * (const T &val, const Matrix <T> & M) 
{
    return M * val;
}


#include "Matrix.hpp"

#endif
