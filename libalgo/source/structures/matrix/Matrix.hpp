// Description: Matrix and basic operators

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


#ifndef Matrix_HPP
#define Matrix_HPP

//#include "Matrix.h"

#include "libalgo/source/const/Const.h"

#include "libalgo/source/exceptions/BadDataException.h"
#include "libalgo/source/exceptions/IndexOutOfBoundException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"
#include "libalgo/source/exceptions/MathMatrixDifferentSizeException.h"

//Forward declaration
template <typename T>
class MathMatrixDifferentSizeException;


//Constructor of the matrix
template <typename T>
Matrix <T> :: Matrix ( const unsigned int rows_count_, const unsigned int columns_count_, const T non_diag_item_val, const T diag_val ) :
        rows_count ( rows_count_ ), columns_count ( columns_count_ ), items ( rows_count_, std::vector <T> ( columns_count_, non_diag_item_val ) )
{
        //Create matrix with eye value on main diagonal
        for ( unsigned int i = 0; i < std::min ( rows_count, columns_count ) ; i++ )
        {
                items[i][i] = diag_val;
        }
}


//Copy constructor
template <typename T>
template <typename U>
Matrix <T> :: Matrix ( const Matrix <U> &M )
        : rows_count ( M.rows() ), columns_count ( M.cols() ), items ()
{
        //Create copy of the matrix
        for ( unsigned int i = 0; i < M.getItems().size(); i++ )
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
                        throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " different rows_count count in operator =. Can not assign matrices  ", *this, M );
                }

                //Matrix dimension are invalid, throw exception
                if ( columns_count != M.columns_count )
                {
                        throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " different columns_count count in operator =. Can not assign matrices.  ", *this, M );
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
Matrix <T> Matrix <T> :: operator + ( const Matrix <U> &M ) const
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
Matrix<T> Matrix <T> :: operator - ( const Matrix < U > &M ) const
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
Matrix<T> Matrix <T> :: operator * ( const Matrix < U > &M ) const
{
        const unsigned int m2 = M.rows(), n2 = M.cols();

        //Matrix dimension invalid, throw exception
        if ( columns_count !=  m2 )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows_count count. Cannot compute A *= B. ", *this, M );
        }

        //Create temporary matrix for results
        Matrix <T> C ( rows_count, n2 );

        //Multiplication of matrices
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                for ( unsigned int j = 0; j < n2; j++ )
                {
                        T sum = 0;

                        for ( unsigned int k = 0; k < columns_count; k++ )
                        {
                                sum += items[i][k] * M ( k, j );
                        }

                        C ( i, j ) = sum;
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
Matrix <T> & Matrix <T> :: operator += ( const Matrix <U> &M )
{
        //Matrix dimension are invalid, throw exception
        if ( this->rows_count !=  M.rows() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows_count count.  Cannot compute A += B. " , *this, M );
        }

        //Matrix dimension are invalid, throw exception
        if ( this->columns_count != M.cols() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different columns_count count.  Cannot compute A += B. ", *this, M );
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
Matrix <T> & Matrix <T> :: operator -= ( const Matrix <U> &M )
{
        //Matrix dimension invalid, throw exception
        if ( rows_count !=  M.rows() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different rows_count count. Cannot compute A -= B.", *this, M );
        }

        //Matrix dimension invalid, throw exception
        if ( columns_count != M.cols() )
        {
                throw MathMatrixDifferentSizeException <Matrix <U> > ( "MathMatrixDifferentSizeException: ", " different columns_count count.  Cannot compute A -= B. ", *this, M );
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
Matrix <T> Matrix <T> ::operator % ( const Matrix <U> &M ) const
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



//Get row / col
template <typename T>
Matrix <T> Matrix <T> ::row ( const unsigned int r ) const
{
        //Matrix dimension invalid, throw exception
        if ( r > rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " row index exceeds rows_count.  " );
        }

        //Create empty row matrix
        Matrix <T> m_temp ( 1, columns_count );

        //Copy items to the matrix
        for ( unsigned int i = 0; i < columns_count; i++ )
        {
                m_temp ( 0, i ) =  items[r][i];
        }

        //Get matrix
        return m_temp;
}


template <typename T>
Matrix < T > Matrix <T> :: col ( const unsigned int c ) const
{
        //Matrix dimension invalid, throw exception
        if ( c > columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: ", " col index exceeds columns_count.  " );
        }

        //Create empty column matrix
        Matrix <T> m_temp ( rows_count, 1 );

        //Copy items to the matrix
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                m_temp ( i, 0 ) = items[i][c];
        }

        //Get matrix
        return m_temp;
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

        //Copy row
        for ( unsigned int i = 0; i < columns_count; i++ )
        {
                items[r][i] = M ( 0, i );
        }
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

        //Copy col
        for ( unsigned int i = 0; i < rows_count; i++ )
        {
                items[i][c] = M ( i, 0 );
        }
}


//Replace part of the matrix with a matrix
template <typename T>
void Matrix <T> ::replace ( const Matrix <T> & M, const unsigned int row, const unsigned int col )
{
        const unsigned int m = M.rows(), n = M.cols();

        if ( m + row > rows_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: a submatrix does not fit at the specified row position, ", "can not append a submatrix to the matrix. " );
        }

        if ( n + col > columns_count )
        {
                throw IndexOutOfBoundException ( "IndexOutOfBoundException: a submatrix does not fit at the specified col position, ", "can not append a submatrix to the matrix." );
        }

        //Copy submatrix
        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        items [i + row][j + col] = M ( i, j );
                }
        }
}

template <typename T>
Matrix <T>  Matrix <T> ::operator () (const Matrix <T> &M, const unsigned int row, const unsigned int col)
{
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

        for ( unsigned int i = 0; i < ( r2 - r1 + 1 ); i++ )
        {
                for ( unsigned int j = 0; j < ( c2 - c1 + 1 ); j++ )
                {
                        M ( i, j ) = items[i + r1 ][j + c1];
                }
        }

        return M;
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
