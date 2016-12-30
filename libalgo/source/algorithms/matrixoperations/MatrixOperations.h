// Description: Basic matrices operation

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


#ifndef MatrixOperations_H
#define MatrixOperations_H

#include <cmath>

#include "libalgo/source/structures/matrix/Matrix.h"


//Result of the Cholesky factorization
typedef enum
{
	upper = 0,
	lower
} TCholeskyFactorization;


//Define namespace: advanced matrix operations
namespace MatrixOperations
{
        //New user defined type: items of vector
        template <typename T>
        struct TVector
        {
                typedef std::vector <T> Type;
        };
               
	template <typename T>
	Matrix <T> eye(const unsigned int r, const unsigned int c, const T mult);
	
	template <typename T>
	Matrix <T> ones(const unsigned int r, const unsigned int c, const T mult);
	
	template <typename T> 
	int sign(const T val);
	
        template <typename T>
	Matrix <T> trans(Matrix <T> A);
        
	template <typename T, typename U>
	bool isequal(const Matrix <T> &A, const Matrix <U> &B);

        template <typename T>
	T min(const Matrix <T> &A, unsigned int & row, unsigned int & col);
        
        template <typename T>
	T min(const Matrix <T> &A);
       
        template <typename T>
	T max(const Matrix <T> &A, unsigned int & row, unsigned int & col);
        
        template <typename T>
	T max(const Matrix <T> &A);
       
        template <typename T>
	Matrix <T> abs(Matrix <T> A);
        
	template <typename T>
	T trace(const Matrix <T> &A);
	
        template <typename T>
	T sumRow(const Matrix <T> &A, const unsigned int row);
       
        template <typename T>
	T sumCol(const Matrix <T> &A, const unsigned int col);
        
        template <typename T>
	Matrix <T> sumRows(const Matrix <T> &A);
        
        template <typename T>
	Matrix <T> sumCols(const Matrix <T> &A);
        
        template <typename T>
	T sum(const Matrix <T> &A);
        
        template <typename T>
	T sum2(const Matrix <T> &A);
        
        template <typename T>
	T norm(const Matrix <T> &A);
        
        template <typename T>
	Matrix <T> sqrtm(const Matrix <T> &A);
        
	template <typename T>
	T median(const Matrix <T>  &A, const unsigned int col = 0);
	
	template <typename T>
	T mad(const Matrix <T> &A, const unsigned int c = 0);
	
        template <typename T>
	Matrix <T> diag(const Matrix <T> &A, const int offset = 0);
       
	template <typename T>
	Matrix <T> tridiag(const Matrix <T> &M, Matrix <T> &l, Matrix <T> &r);
	
	template <typename T>
	Matrix <T> tridiagh(Matrix <T> A);
	
        template <typename T>
	void sort(Matrix <T> &A, Matrix <unsigned int> &IX, const int c = -1);
        
	template <typename T>
	void sortrows(Matrix <T> &A, Matrix <unsigned int> &IX, const int c = -1);
	
        template <typename T>
	Matrix <T> gem(Matrix <T> A, Matrix <T> B);
        
        template <typename T>
	void gemF(Matrix <T> &A_triangle, Matrix <T> *B = NULL);
        
        template <typename T>
	void gemFSign(Matrix <T> &A_triangle, Matrix <T> *B = NULL);
        
        template <typename T>
	void gemB(Matrix <T> &A, const Matrix <T> &B, Matrix <T> &X);
        
        template <typename T>
	T det(Matrix <T> A);
        
        template <typename T>
	void lu(Matrix <T> A, Matrix <T> &L, Matrix <T> &U, Matrix <T> &P);
       
        template <typename T>
	void qr(const Matrix <T> &A, Matrix <T> &Q, Matrix <T> &R);
       
        template <typename T>
	void qr(const Matrix <T> &A, Matrix <T> &Q, Matrix <T> &R, Matrix <unsigned int> &P);
        
	template <typename T>
	bool posdef(const Matrix <T> &A);
	
	template <typename T>
	Matrix <T> chol(Matrix <T> A, const TCholeskyFactorization f);
	
	template <typename T>
	void gill(const Matrix <T> &A, Matrix <T> &R, Matrix <T> &E, bool &indefinite, const T max_error = 2.2204e-16);
	
	template <typename T>
	void gill(const Matrix <T> &A, Matrix <T> &L, Matrix <T> &d, Matrix <T> &e, Matrix <T> &v, const T max_error = 2.2204e-16);
	
        template <typename T>
	Matrix <T> inv(Matrix <T> A);
        
        template <typename T>
	Matrix <T> pinv(const Matrix <T> &A);
       
        template <typename T>
	Matrix <T> mlsqr(const Matrix <T> &A, const Matrix <T> &W, const Matrix <T> &L);
       
	template <typename T>
	void eig(const Matrix <T> &A, Matrix<T> &V, Matrix<T> &L, const unsigned int max_iter = 1000, const T tolerance = 1.0e-15);
	
        template <typename T>
	Matrix <T> pinvs(const Matrix <T> &A, const T tolerance = MAX_FLOAT_OPER_ERROR);
       
        template <typename T>
	Matrix <T> pinv1(const Matrix <T> &A);
        
        template <typename T>
	Matrix <T> pinv2(const Matrix <T> &A);
        
        template <typename T>
	void svd(const Matrix <T> &A, Matrix <T> &U, Matrix <T> &B, Matrix <T> &V, unsigned int max_iterations);
      
        template <typename T>
	void bidiag(const Matrix <T> &A, Matrix <T> &U, Matrix <T> &B, Matrix <T> &V);
        
        template <typename T>
	void qr0Shift(Matrix <T> &D, Matrix <T> &E, Matrix <T> &U, Matrix <T> &V);
        
        template <typename T>
	Matrix <T> hous(const Matrix <T> &A, const unsigned int i, const unsigned int j);
        
        template <typename T>
	Matrix <T> hous(const Matrix <T> &A, const unsigned int i);
        
        template <typename T>
	void givens(const T f, const T g, T &c, T &s, T &r);
        
        template <typename T>
	void updgivens(const T cs, const T sn, Matrix <T> & V1, Matrix <T> &V2);
       
        template <typename T>
	Matrix <T> load(const char *file, const bool print_exception = true, std::ostream * output = &std::cout);
       
};

#include "MatrixOperations.hpp"

#endif
