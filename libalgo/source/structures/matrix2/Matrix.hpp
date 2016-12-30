// Description: Matrix and basic operators

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


#ifndef Matrix_HPP
#define Matrix_HPP

#include "libalgo/source/const2/Const.h"

#include "libalgo/source/exceptions/BadDataException.h"
#include "libalgo/source/exceptions/IndexOutOfBoundException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"
#include "libalgo/source/exceptions/MathMatrixDifferentSizeException.h"


//Constructor of the matrix
template <typename T>
Matrix <T> :: Matrix ( const unsigned int rows_count_, const unsigned int columns_count_, const T non_diag_val, const T diag_val ) :
        rows_count ( rows_count_ ), columns_count ( columns_count_ ), items ( rows_count_, std::vector <T> ( columns_count_, non_diag_val ) )
{
        //Create matrix with a specific value on the main diagonal
        for ( unsigned int i = 0; i < std::min ( rows_count, columns_count ) ; i++ )
        {
                items[i][i] = diag_val;
        }
}


//Create matrix from 2D array
template <typename T>
template <std::size_t rows_count_, std::size_t columns_count_>
Matrix <T> :: Matrix(const std::array<std::array<T, columns_count_>, rows_count_> &data) : 
	rows_count(rows_count_), columns_count(columns_count_), items(rows_count_, std::vector <T>(columns_count_, 0))
{
	//Process all items
	for (int i = 0; i < rows_count; i++)
	{
		for (int j = 0; j < columns_count; j++)
		{
			items[i][j] = data[i][j];
		}
	}
}


//Copy constructor
template <typename T>
template <typename U>
Matrix <T> :: Matrix ( const Matrix <U> &M )
        : rows_count ( M.rows() ), columns_count ( M.cols() ), items ()
{
	const unsigned int n = M.rows();

        //Create copy of the matrix
        for ( unsigned int i = 0; i < n; i++ )
        {
                items.push_back ( std::vector <T> ( M.getItems() [i].begin(), M.getItems() [i].end() ) );
        }
}


//Destructor
template <typename T>
Matrix <T> :: ~Matrix() {}


//Assignment operator =
template <typename T>
Matrix <T> & Matrix <T> ::operator = ( const Matrix <T> &M )
{
        if ( this != &M )
        {

                //Matrix dimension are invalid, throw exception
                if ( rows_count !=  M.rows_count )
                {
                        throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " different rows count in operator =. Can not assign matrices  ", *this, M );
                }

                //Matrix dimension are invalid, throw exception
                if ( columns_count != M.columns_count )
                {
                        throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " different columns count in operator =. Can not assign matrices.  ", *this, M );
                }

                //Process all lines
                for ( unsigned int i = 0; i < rows_count; i++ )
                {
                        //Process columns_count of the actual row
                        for ( unsigned int j = 0 ; j < columns_count; j++ )
                        {
                                //Add to vector
                                items[i][j] = M.items[i][j];
                        }
                }
        }

        return *this;
}



//Matrix operators + : Matrix + Matrix
template <typename T>
template <typename U>
Matrix <T> Matrix <T> :: operator + (  Matrix <U> const &M ) const
{
        return Matrix <T> ( *this ).operator += ( M );

}

//Matrix operators + : Matrix + scalar
template <typename T>
template <typename U>
Matrix <T> Matrix <T> :: operator + (const U val) const
{
	const unsigned int m = this->rows_count, n = this->columns_count;
	const Matrix <T> A(m, n, 1, 1);

	return Matrix <T>(*this).operator += (A * val);

}


//Matrix operators - : Matrix - Matrix
template <typename T>
template <typename U>
Matrix<T> Matrix <T> :: operator - ( Matrix < U > const &M ) const
{
        return Matrix < T > ( *this ).operator -= ( M );
}


//Matrix operators - : Matrix - scalar
template <typename T>
template <typename U>
Matrix <T> Matrix <T> :: operator - (const U val) const
{
	const unsigned int m = this->rows_count, n = this->columns_count;
	const Matrix <T> A(m, n, 1, 1);

	return Matrix <T>(*this).operator -= (A * val);
}


//Matrix operators * : Matrix * Matrix
template <typename T>
template <typename U>
Matrix<T> Matrix <T> :: operator * ( Matrix < U > const &M ) const
{
        const unsigned int m2 = M.rows(), n2 = M.cols();

        //Matrix dimension invalid, throw exception
        if ( columns_count !=  m2 )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows and columns count. Cannot compute A *= B. ", *this, M );
        }

        //Create temporary matrix for results
        Matrix <T> C ( rows_count, n2 );
	
	//Multiplication of matrices
	for (unsigned int i = 0; i < rows_count; i++)
	{
		for (unsigned int k = 0; k < columns_count; k++)
		{
			for (unsigned int j = 0; j < n2; j++)
			{
				C(i, j) += items[i][k] * M(k, j);
			}
		}
	}
	
        return C;
}


//Matrix operators * : Matrix * Scalar
template <typename T>
template <typename U>
Matrix<T> Matrix <T> :: operator * ( const U & val ) const
{
        return Matrix < T > ( *this ).operator *= ( val );
}


//Matrix operators / : Matrix / Scalar
template <typename T>
template <typename U>
Matrix<T> Matrix <T> :: operator / ( const U & val ) const
{
        return Matrix < T > ( *this ).operator /= ( val );
}


//Matrix operators += : Matrix += Matrix
template <typename T>
template <typename U>
Matrix <T> & Matrix <T> :: operator += ( Matrix <U> const &M )
{
        //Matrix dimension are invalid, throw exception
        if ( this->rows_count !=  M.rows() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows count.  Cannot compute A += B. " , *this, M );
        }

        //Matrix dimension are invalid, throw exception
        if ( this->columns_count != M.cols() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different columns count.  Cannot compute A += B. ", *this, M );
        }

        //Process all lines
        for ( unsigned int i = 0;  i < rows_count; i++ )
        {
                //Process columns_count of the actual row
                for ( unsigned int j = 0 ; j < columns_count; j++ )
                {
                        //Add to vector
                        items[i][j] = items[i][j]  + M ( i, j );
                }
        }

        return *this;
}


//Matrix operators -= : Matrix -= Matrix
template <typename T>
template <typename U>
Matrix <T> & Matrix <T> :: operator -= ( Matrix <U> const &M )
{
        //Matrix dimension invalid, throw exception
        if ( rows_count !=  M.rows() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows count. Cannot compute A -= B.", *this, M );
        }

        //Matrix dimension invalid, throw exception
        if ( columns_count != M.cols() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different columns count.  Cannot compute A -= B. ", *this, M );
        }

        //Process all lines
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                //Process columns_count of the actual row
                for ( unsigned int j = 0 ; j < columns_count; j++ )
                {
                        //Add to vector
                        items[i][j] = items[i][j] - M ( i, j );
                }
        }

        return *this;
}


//Matrix operators *= : Matrix *= Scalar
template <typename T>
template <typename U>
Matrix <T> & Matrix <T> ::operator *= ( const U & val )
{
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                for ( unsigned int j = 0; j < columns_count; j++ )
                {
                        items[i][j] = val * items[i][j];
                }
        }

        return *this;
}


//Matrix operators /= : Matrix /= Scalar
template <typename T>
template <typename U>
Matrix <T> & Matrix <T> ::operator /= ( const U & val )
{
        //Divider is zero
        if ( fabs ( ( T ) val ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: divider is zero. ",  "Can not divide matrix by scalar.", val );
        }

        //Multiplication of matrices
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                for ( unsigned int j = 0; j < columns_count; j++ )
                {
                        items[i][j] = items[i][j] / val;
                }
        }

        return *this;
}



//Matrix operators .* : Matrix .* Matrix (Hadamard product)
template <typename T>
template <typename U>
Matrix <T> Matrix <T> ::operator % ( Matrix <U> const &M ) const
{
        const unsigned int m2 = M.rows(), n2 = M.cols();

        //Matrix dimension invalid, throw exception
        if ( rows_count !=  m2 )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows_count count. Cannot compute A .* B.", *this, M );
        }

        //Matrix dimension invalid, throw exception
        if ( columns_count != n2 )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different columns_count count.  Cannot compute A .* B. ", *this, M );
        }

	//Create temporary matrix for results
	Matrix <T> C(m2, n2);

        //Hadamard product
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                for ( unsigned int j = 0; j < columns_count; j++ )
                {
                        C(i, j) = items[i][j] * M ( i, j );
                }
        }

        return C;
}



//Get row of the matrix
template <typename T>
Matrix <T> Matrix <T> ::row ( const unsigned int r ) const
{
        //Matrix dimension invalid, throw exception
        if ( r > rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " row index exceeds rows_count.  " );
        }

	//Call operator A(r1, r2, c1, c2)
	return (*this)(r, r, 0, columns_count - 1);
}


//Get column of the matrix
template <typename T>
Matrix < T > Matrix <T> :: col ( const unsigned int c ) const
{
        //Matrix dimension invalid, throw exception
        if ( c > columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " col index exceeds columns_count.  " );
        }

	//Call operator A(r1, r2, c1, c2)
	return (*this)(0, rows_count - 1, c, c);
}


//Set row / col
template <typename T>
void Matrix <T> ::row ( const Matrix <T> &M, const unsigned int r )
{
        //Matrix dimension invalid, throw exception
        if ( r > rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " row index exceeds rows_count.  " );
        }

        //Matrix dimension invalid, throw exception
        if ( M.cols() != columns_count )
        {
                throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " invalid dimension of the matrices (columns_count count).  ", *this, M );
        }

	//Call operator A(M, r, c);
	(*this)(M, r, 0);
}


//Set submatrix
template <typename T>
void Matrix <T> ::col ( const Matrix <T> &M, const unsigned int c )
{
        //Matrix dimension invalid, throw exception
        if ( c > columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " col index exceeds columns_count.  " );
        }

        //Matrix dimension are invalid, throw exception
        if ( M.rows() != rows_count )
        {
                throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " invalid dimension of the matrices (rows_count count).  ", *this, M );
        }

	//Call operator A(M, r, c);
	(*this)(M, 0, c);
}


//Resize matrix
template <typename T>
void Matrix <T> ::resize ( const int r, const int c, const T val )
{
        //Not enough rows
        if ( r < 1 )
        {
                throw BadDataException ( "BadDataException: bad rows count", "matrix can not be resized " );
        }

        //Not enough cols
        if ( c < 1 )
        {
                throw BadDataException ( "BadDataException: bad columns count", "matrix can not be resized " );
        }

        //Matrix dimensions have been changed
        if ( ( r  != rows_count ) || ( c != columns_count ) )
        {
                //Resize rows
                items.resize ( r );

                //Resize columns
                for ( unsigned int i = 0; i < r ; i++ ) items[i].resize ( c, val ) ;

                //Set new dimensions
                rows_count = r;
                columns_count = c;
        }
}


//Matrix operator ()()
template <typename T>
T & Matrix <T> :: operator() ( const unsigned int row, unsigned int col )
{
        //Matrix dimension invalid, throw exception
        if ( row >= rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", "invalid dimension of the matrix (row > rows_count), row = " );
        }

        //Matrix dimension are invalid, throw exception
        if ( col >= columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", "invalid dimension of the matrix (col > columns_count), col = " );
        }

        return items[row][col];
}


//Matrix operator ()()
template <typename T>
T const & Matrix <T> :: operator() ( const unsigned int row, unsigned int col ) const
{
        //Matrix dimension are invalid, throw exception
        if ( row >= rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", "invalid dimension of the matrix (row > rows_count), row = " );
        }

        //Matrix dimension are invalid, throw exception
        if ( col >= columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", "invalid dimension of the matrix (col > columns_count), col = " );
        }

        return items[row][col];
}


//Matrix operator (r1, r2, c1, c2): get submatrix of the matrix
template <typename T>
Matrix <T>  Matrix <T> ::operator () ( const unsigned int r1, const unsigned int r2, const unsigned int c1, const unsigned int c2 ) const
{
        //Bad row index
        if ( r2 >  rows_count )
                throw BadDataException ( "BadDataException: row index r2 must not be greater than A.rows_count. ", " Can not create the submatrix. " );

        //Bad row index
        if ( c2 > columns_count )
                throw BadDataException ( "BadDataException: col index c2 must not be greater than A.columns_count. ", " Can not create the submatrix. " );

        //Bad row index interval
        if ( r1 > r2 )
                throw BadDataException ( "BadDataException: row index r2 must not be smaller then r1. ", " Can not create the submatrix. " );

        //Bad col index interval
        if ( c1 > c2 )
                throw BadDataException ( "BadDataException: col index c2 must not be smaller then c1. ", " Can not create the submatrix. " );

        //Create sub-matrix
        Matrix <T> M ( r2 - r1 + 1, c2 - c1 + 1 ) ;

        for ( unsigned int i = r1; i <= r2; i++ )
        {
                for ( unsigned int j = c1; j <= c2; j++ )
                {
                        M ( i - r1, j - c1 ) = items[i][j];
                }
        }

        return M;
}


template <typename T>
void Matrix <T> ::operator () (const Matrix <T> &M, const unsigned int row, const unsigned int col)
{
	//Matrix operator (M, r, c): replace part of the matrix A at the position [row, col] with M
	const unsigned int m = M.rows(), n = M.cols();

	if (m + row > rows_count)
	{
		throw IndexOutOfBoundException("IndexOutOfBoundException: a submatrix does not fit at the specified row position, ", "can not append a submatrix to the matrix. ");
	}

	if (n + col > columns_count)
	{
		throw IndexOutOfBoundException("IndexOutOfBoundException: a submatrix does not fit at the specified col position, ", "can not append a submatrix to the matrix.");
	}

	//Copy submatrix
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < n; j++)
		{
			items[i + row][j + col] = M(i, j);
		}
	}
}

template <typename T>
void Matrix <T> ::print ( std::ostream * output ) const
{
        //Print matrix
        *output << std::showpoint << std::fixed << std::right;
        *output << '\n';

        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                *output << "| ";

                for ( unsigned int j = 0; j < columns_count; j++ )
                {
                        *output <<  std::setw ( 18 ) << std::setprecision ( 9 );
                        items[i][j] < MAX_FLOAT ? *output << items[i][j] : *output << "---";
                }

                *output << " |" << '\n';
        }
}


#endif
