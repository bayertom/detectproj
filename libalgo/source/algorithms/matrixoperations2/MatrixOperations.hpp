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


#ifndef MatrixOperations_HPP
#define MatrixOperations_HPP

#include <cmath>

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TVector2D.h"

#include "libalgo/source/comparators/indexComp.h"

#include "libalgo/source/io/File.h"

#include "libalgo/source/exceptions/MathOverflowException.h"
#include "libalgo/source/exceptions/MathMatrixDifferentSizeException.h"
#include "libalgo/source/exceptions/MathMatrixSingularException.h"
#include "libalgo/source/exceptions/MathMatrixNotSquareException.h"
#include "libalgo/source/exceptions/MathMatrixNotPositiveDefiniteException.h"

namespace MatrixOperations
{
	//Identity matrix
	template <typename T>
	Matrix <T> eye(const unsigned int r, const unsigned int c, const T mult)
	{
		Matrix <T> A(r, c, 0, 1);
		return mult * A;
	}


	//Matrix of ones
	template <typename T>
	Matrix <T> ones(const unsigned int r, const unsigned int c, const T mult)
	{
		Matrix <T> A(r, c, 1, 1);
		return mult * A;
	}


	//Signum
	template <typename T>
	int sign(const T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	//Transpose matrix
	template <typename T>
	Matrix <T> trans(const Matrix <T> &A)
	{
		//Create trans matrix
		const unsigned int m = A.rows(), n = A.cols();
		Matrix <T> A_trans(n, m);

		//Copy m to m_trans
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				A_trans(j, i) = A(i, j);
			}
		}

		return A_trans;
	}


	//Are two matrices equal ?
	template <typename T, typename U>
	bool isequal(const Matrix <T> &A, const Matrix <U> &B)
	{

		const unsigned int m1 = A.rows(), n1 = A.cols(), m2 = B.rows(), n2 = B.cols();

		//Rectangular matrix
		if (m1 != n1 || m2 != n2)
		{
			return false;
		}

		//Process matrix element by element
		for (unsigned int i = 0; i < m1; i++)
		{
			for (unsigned int j = 0; j < n1; j++)
			{
				if (fabs(A(i, j) - B(i, j)) > MIN_FLOAT)
					return false;
			}
		}

		return true;
	}


	//Find min value, row and col position
	template <typename T>
	T min(const Matrix <T> &A, unsigned int & row, unsigned int & col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Initialize minimum
		T min_value = A(0, 0);

		//Process all items
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				if (A(i, j) < min_value)
				{
					min_value = A(i, j);
					row = i; col = j;
				}
			}
		}

		return min_value;
	}


	//Find min value
	template <typename T>
	T min(const Matrix <T> &A)
	{
		unsigned int row_index, column_index;
		return min(A, row_index, column_index);
	}


	//Find max value, row and col index
	template <typename T>
	T max(const Matrix <T> &A, unsigned int & row, unsigned int & col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Initialize maximum
		T max_value = A(0, 0);

		//Process all items
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				if (A(i, j) > max_value)
				{
					max_value = A(i, j);
					row = i; col = j;
				}
			}
		}

		return max_value;
	}


	//Find max value
	template <typename T>
	T max(const Matrix <T> &A)
	{
		unsigned int row_index, column_index;
		return max(A, row_index, column_index);
	}


	//Abolute value of a matrix
	template <typename T>
	Matrix <T> abs(const Matrix <T> &A)
	{
		//Create trans matrix
		const unsigned int m = A.rows(), n = A.cols();
		Matrix <T> A_a(m, n);

		//Copy m to m_trans
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				A_a(i, j) = std::fabs(A(i, j));
			}
		}

		return A_a;
	}


	//Trace of the matrix
	template <typename T>
	T trace(const Matrix <T> &A)
	{
		return sum(diag(A));
	}

	//Sum selected row
	template <typename T>
	T sumRow(const Matrix <T> &A, const unsigned int row)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Matrix row index is invalid, throw exception
		if (row < 0 || row > m)
		{
			throw IndexOutOfBoundException("IndexOutOfBoundException: ", "can not compute row sum, invalid row (row < 0 || row > rows_count)");
		}

		T sum = 0;

		//Sum all items in a row
		for (unsigned int j = 0; j < n; j++)
		{
			sum += A(row, j);
		}

		return sum;
	}


	//Sum selected col
	template <typename T>
	T sumCol(const Matrix <T> &A, const unsigned int col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Matrix  col index is invalid, throw exception
		if (col < 0 || col > n)
		{
			throw IndexOutOfBoundException("IndexOutOfBoundException: ", "can not compute col sum, invalid col (col < 0 || col > columns_count)");
		}

		T sum = 0;

		//Sum all items in a column
		for (unsigned int i = 0; i < m; i++)
		{
			sum += A(i, col);
		}

		return sum;
	}


	//Sum of rows
	template <typename T>
	Matrix <T> sumRows(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		Matrix <T> SR(m, 1);

		//Sum all rows
		for (unsigned int i = 0; i < m; i++)
		{
			SR(i, 0) = sumRow(A, i);
		}

		return SR;
	}


	//Sum of columns
	template <typename T>
	Matrix <T> sumCols(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		Matrix <T> SC(1, n);

		//Sum all columns
		for (unsigned int i = 0; i < n; i++)
		{
			SC(0, i) = sumCol(A, i);
		}

		return SC;
	}



	//Sum of the matrix items
	template <typename T>
	T sum(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		T sum = 0;

		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				sum += A(i, j);
			}
		}

		return sum;
	}


	//Sum2 of the matrix items
	template <typename T>
	T sum2(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		T sum = 0;

		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				sum += A(i, j) * A(i, j);
			}
		}

		return sum;
	}


	//Norm of the matrix
	template <typename T>
	T norm(const Matrix <T> &A)
	{
		return sqrt(sum2(A));
	}


	//Sqrt of the matrix
	template <typename T>
	Matrix <T> sqrtm(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		Matrix <T> B(m, n);

		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				B(i, j) = sqrt(A(i, j));
			}
		}

		return B;
	}


	//Median of a column vector
	template <typename T>
	T median(const Matrix <T>  &A, const unsigned int col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Matrix  col index is invalid, throw exception
		if (col < 0 || col > n)
		{
			throw IndexOutOfBoundException("IndexOutOfBoundException: ", "can not compute median, invalid col (col < 0 || col > columns_count)");
		}

		//Matrix contains one row
		if (m == 1)
		{
			return A(0, col);
		}

		//Sort specified column of matrix
		Matrix <T> AS = A;
		Matrix <unsigned int> IX(m, n);
		sort(AS, IX, col);

		//Even case
		if (m % 2 == 0)
		{
			return 0.5 * (AS(m / 2, col) + AS(m / 2 - 1, col));
		}

		//Odd case
		else
		{
			return AS(m / 2, col);
		}
	}


	//Median absolute deviation
	template <typename T>
	T mad(const Matrix <T> &A, const unsigned int c)
	{
		const unsigned int m = A.rows(), n = A.cols();

		return median(abs(A(0, m - 1, c, c) - median(A, c)));
	}



	//Get matrix diagonal and convert to column vector
	template <typename T>
	Matrix <T> diag(const Matrix <T> &A, const int offset)
	{
		const int m = A.rows(), n = A.cols();

		//A is a matrix
		if ((n > 1) || (m == 1) && (n == 1))
		{
			//Bad offset: to large
			if (offset > n - 1)
			{
				throw BadDataException("BadDataException: offset value to large, offset > n. ", "Can not get diagonal of the matrix.");
			}

			//Bad offset: to low
			if (offset + m < 0)
			{
				throw BadDataException("BadDataException: offset value to small, offset + n < 0. ", "Can not get diagonal of the matrix.");
			}

			//Set size of the diagonal
			const int mn = (offset > 0 ? m - offset : n + offset);

			//Create output diagonal matrix (m, 1)
			Matrix <T> A_diag(mn, 1);

			//Copy A elements to the diagonal matrix A_diag
			for (unsigned int i = std::abs(offset); (i < m) && (i < n); i++)
			{
				const T item = (offset >= 0 ? A(i - offset, i) : A(i, i + offset));
				A_diag(i - std::abs(offset), 0) = item;
			}

			return A_diag;
		}

		//A is a column vector
		const int mn = m + std::abs(offset);

		//Create new diagonal matrix
		Matrix <T> A_diag(mn, mn);

		//Copy A elements to the diagonal matrix A_diag
		for (unsigned int i = 0; i < m; i++)
		{
			if (offset >= 0) A_diag(i, i + offset) = A(i, 0);
			else A_diag(i + std::abs(offset), i) = A(i, 0);
		}

		return A_diag;
	}


	//Compute tridiagonal matrix using the Lanczos bidiagonalization
	//Works for a general sqaure matrix
	template <typename T>
	Matrix <T> tridiag(const Matrix <T> &M, Matrix <T> &l, Matrix <T> &r)
	{
		//Slightly modified algorithm by Z. Bai, 1997 
		//L, R are left and right initial vectors
		const unsigned int m = M.rows(), n = M.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute tridiagonal matrix; (rows_count, columns_count):  ", M);
		}

		const unsigned int m2 = l.rows(), m3 = r.rows();

		//Different amount of rows
		if (m != m2)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different amount of rows of M, L, can not compute tridiagonal matrix:  ", M, l);
		}

		//Different amount of rows
		if (m != m3)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different amount of rows of M, R, can not compute tridiagonal matrix:  ", M, r);
		}

		//Left and right Lanczos vector
		unsigned int i = 0;
		Matrix <T> W(m, m), V(m, m);

		//Parameters of the tridiagonal matrix
		Matrix <T> A(m, 1), B(m, 1), G(m, 1), D(m, 1), E(m, 1), R(m, 1);

		//Initialize parameters A - R
		R(i, 0) = norm(r);
		E(i, 0) = norm(l);
		V(r * (1.0 / R(i, 0)), 0, i);
		W(l * (1.0 / E(i, 0)), 0, i);
		D(i, 0) = norm(trans(W(0, m - 1, i, i)) * V(0, m - 1, i, i));

		Matrix <T> MV = M * V(0, m - 1, i, i);
		A(i, 0) = norm(trans(W(0, m - 1, i, i)) * MV * (1.0 / D(i, 0)));

		Matrix <T> v = MV - V(0, m - 1, i, i)* A(i, 0);
		R(i + 1, 0) = norm(v);
		V(v * (1 / R(i + 1, 0)), 0, i + 1);

		Matrix <T> MW = trans(M) * W(0, m - 1, i, i);
		Matrix <T> w = MW - W(0, m - 1, i, i) * A(i, 0);
		E(i + 1, 0) = norm(w);
		W(w * (1.0 / E(i + 1, 0)), 0, i + 1);

		//Perform the Lanczos method
		for (i = 1; i < m; i++)
		{
			//Compute parameters A-G
			D(i, 0) = norm(trans(W(0, m - 1, i, i)) * V(0, m - 1, i, i));
			MV = trans(M) * V(0, m - 1, i, i);
			A(i, 0) = norm(trans(W(0, m - 1, i, i)) * MV * (1.0 / D(i, 0)));
			B(i, 0) = D(i, 0) / D(i - 1, 0) * E(i, 0);
			G(i, 0) = D(i, 0) / D(i - 1, 0) * R(i, 0);

			//construct left and right Lanczos vectors
			if (i < n - 1)
			{
				v = MV - V(0, m - 1, i, i) * A(i, 0) - V(0, m - 1, i - 1, i - 1) * B(i, 0);

				//Perform re-orthogonalization
				for (unsigned int j = 0; j <= i; j++)
				{
					v = v - V(0, m - 1, j, j) * trans(W(0, m - 1, j, j)) * v * (1.0 / D(i, 0));
				}

				//Right Lanzos vector
				R(i + 1, 0) = norm(v);
				V(v * (1.0 / R(i + 1, 0)), 0, i + 1);

				MW = trans(M) * W(0, m - 1, i, i);
				w = MW - W(0, m - 1, i, i) * A(i, 0) - W(0, m - 1, i - 1, i - 1) * G(i, 0);

				//Perform re-orthogonalization
				for (unsigned int j = 0; j <= i; j++)
				{
					w = w - W(0, m - 1, j, j) * trans(V(0, m - 1, j, j)) * w * (1.0 / D(i, 0));
				}

				//Left Lanzos  vector
				E(i + 1, 0) = norm(w);
				W(w * (1.0 / E(i + 1, 0)), 0, i + 1);
			}
		}

		//Construct tridiagonal matrix
		Matrix <T> TR(m, n);
		for (unsigned int i = 1; i < m; i++)
		{
			TR(i - 1, i - 1) = A(i - 1, 0);
			TR(i, i - 1) = R(i, 0);
			TR(i - 1, i) = B(i, 0);
		}

		TR(m - 1, m - 1) = A(m - 1, 0);

		//const Matrix<T> TR = diag(A) + diag(R(1, m - 1, 0, m - 1), -1) + diag(B(1, m - 1, 0, m - 1), 1);

		return TR;
	}


	//Compute tridiagonal matrix using the Householder rotations
	//Works for positive deinite matrices
	template <typename T>
	Matrix <T> tridiagh(const Matrix <T> &M)
	{
		//Slightly modified algorithm by Matt Fig, 2008
		Matrix <T> A = M;
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute tridiagonal matrix; (rows_count, columns_count):  ", A);
		}

		//Symmetric matrix
		if (!isequal(A, trans(A)))
		{
			throw BadDataException("BadDataException: Matrix A is not symmetric. ", "Can not compute tridiagonal matrix.");
		}

		//Auxilliary matrices
		Matrix <T> V(m, 1);
		Matrix <T> I(m, m, 0, 1);

		//Perform tridiagonalization
		for (unsigned int i = 0; i < m - 2; i++)
		{
			//Reset elements till the procesed index
			for (unsigned int j = 0; j <= i; j++)
				V(j, 0) = 0;

			//Compute sum
			const T sumA2 = sqrt(norm(trans(A(i + 1, m - 1, i, i)) * A(i + 1, m - 1, i, i)));

			//Perform Householder rotations
			V(i + 1, 0) = sqrt(.5 * (1 + fabs(A(i + 1, i)) / (sumA2 + 2 * MIN_FLOAT)));

			//Create temporary matrix
			const Matrix <T> VT = A(i + 2, m - 1, i, i) * sign(A(i + 1, i)) / (2 * V(i + 1, 0) * sumA2 + 2 * MIN_FLOAT);

			//Update V 
			V(VT, i + 2, 0);

			//Compute rotation
			const Matrix <T> P = I - V * trans(V) * (2);

			//Update matrix
			Matrix <T> AN = P * A * P;
			A = AN;
		}

		return A;
	}


	//Sort matrix by columns and return array of indices
	template <typename T>
	void sort(Matrix <T> &A, Matrix <unsigned int> &IX, const int c)
	{
		//Sort matrix by a column c, if c = -1 sort all columns
		const int m1 = A.rows(), n1 = A.cols();
		const int m2 = IX.rows(), n2 = IX.cols();

		//Test number of rows A and IX
		if (m1 != m2)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different rows count for A, IX, can not sort a matrix:  ", A, IX);
		}

		//Test number of columns A and IX
		if (n1 != n2)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different columns count for A, IX, can not sort a matrix:  ", A, IX);
		}

		//Test, if c is correct
		if ((c > n1 - 1) || (c < -1))
		{
			throw IndexOutOfBoundException("IndexOutOfBoundException: ", "can not sort matrix by a column c, invalid c (c < -1 || col > columns_count)");
		}

		//Initialize indices: i1 = start, i2 = end
		int i1 = (c >= 0 ? c : 0);
		int i2 = (c >= 0 ? c + 1 : n1);

		//Sort columns
		for (; i1 < i2; i1++)
		{
			//Get a column of a matrix
			Matrix <T> col = trans(A(0, m1 - 1, i1, i1));

			//Convert to 2D std::vector
			TVector2D <T> col2D = col.getItems();

			//Convert to 1D vector
			TVector <T>  col1D;
			col1D.insert(col1D.end(), col2D[0].begin(), col2D[0].end());

			//Create 1D vector of indices
			std::vector <unsigned int> ix(m1);
			for (unsigned int j = 0; j < m1; j++)  ix[j] = j;

			//Sort 1D vector and change indices
			std::sort(ix.begin(), ix.end(), indexComp <typename TVector <T>::const_iterator>(col1D.begin(), col1D.end()));

			//Add sorted column to A and sorted indices to IX
			for (unsigned int j = 0; j < m1; j++)
			{
				A(j, i1) = col1D[ix[j]];
				IX(j, i1) = ix[j];
			}
		}
	}


	//Sort rows of the matrix by a specific column
	template <typename T>
	void sortrows(Matrix <T> &A, Matrix <unsigned int> &IX, const int c)
	{
		//Sort matrix by a column c, if c = -1 sort all columns
		const int m1 = A.rows(), n1 = A.cols();
		const int m2 = IX.rows(), n2 = IX.cols();

		//Test number of rows A and IX
		if (m1 != m2)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different rows count for A, IX, can not sort rows of a matrix:  ", A, IX);
		}

		//Test number of columns A and IX
		if (n2 != 1)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " IX does not have 1 column, can not sort rows of a matrix:  ", A, IX);
		}

		//Test, if c is correct
		if ((c > n1 - 1) || (c < -1))
		{
			throw IndexOutOfBoundException("IndexOutOfBoundException: ", "can not sort rows of the matrix by a column c, invalid c (c < -1 || col > columns_count)");
		}

		//Get a column of a matrix
		Matrix <T> col = trans(A(0, m1 - 1, c, c));

		//Convert to 2D std::vector
		TVector2D <T> col2D = col.getItems();

		//Convert to 1D vector
		TVector <T>  col1D;
		col1D.insert(col1D.end(), col2D[0].begin(), col2D[0].end());

		//Create 1D vector of indices
		std::vector <unsigned int> ix(m1);

		for (unsigned int j = 0; j < m1; j++)  ix[j] = j;

		//Sort 1D vector
		std::sort(ix.begin(), ix.end(), indexComp <typename TVector <T>::const_iterator>(col1D.begin(), col1D.end()));

		//Sort rows of the matrix
		Matrix <T> AP = A;
		for (unsigned int i = 0; i < m1; i++)
		{
			//Store permutation
			IX(i, 0) = ix[i];

			//Copy row to permutated matrix
			AP(A(ix[i], ix[i], 0, n1 - 1), i, 0);
		}

		//Assign permutated matrix
		A = AP;
	}


	//Gaussian elimination algorithm with the partial pivotation
	template <typename T>
	Matrix <T> gem(const Matrix <T> &A, Matrix <T> &B)
	{
		const unsigned int m1 = A.rows(), n1 = A.cols();
		const unsigned int m2 = B.rows(), n2 = B.cols();

		//Different rows and cols count
		if (m1 != m2)
		{
			throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " invalid dimension of the matrices, can not perform GEM; (rows_count != columns_count):  ", A, B);
		}

		//Rectangle matrix
		if (m1 != n1)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not perform GEM; (rows_count, columns_count):  ", A, B);
		}

		//Create matrix X
		Matrix <T> X(m1, 1);

		//Forward step of the Gauss elimination: convert A to triangle matrix
		Matrix <T> A_triangle(A);
		gemF(A_triangle, B);

		//Backward step of the Gaussian elimination: compute A * x = B
		gemB(A_triangle, B, X);

		//Return solution
		return X;
	}


	//Gaussian elimination algorithm (Forward step for matrices A, B, store sign in A(m-1,0) item)
	template <typename T>
	void gemF(Matrix <T> &A_triangle, Matrix <T> &B)
	{
		const unsigned int m = A_triangle.rows(), n = A_triangle.cols();

		//Create triangle matrix A_triangle from A
		for (int k = 0; k < n; k++)
		{
			//Initialize index maximum item
			unsigned int index_max = k;

			//Find max index
			for (unsigned j = k + 1; j < n; j++)
			{
				if (fabs(A_triangle(j, k)) > fabs(A_triangle(index_max, k)))
				{
					//Remember index of the maximum item
					index_max = j;
				}
			}

			//Test, if maximum index has been changed
			if (k != index_max)
			{
				//Swap rows in matrix A
				for (unsigned int j = k; j < n; j++)
				{
					const T pom = A_triangle(k, j);
					A_triangle(k, j) = A_triangle(index_max, j);
					A_triangle(index_max, j) = pom;
				}

				//Swap rows in matrix B
				const T pom = B(k, 0);
				B(k, 0) = B(index_max, 0);
				B(index_max, 0) = pom;
			}

			for (unsigned int i = k + 1; i < n; i++)
			{
				//Singularity test
				if (fabs(A_triangle(k, k)) > MIN_FLOAT * fabs(A_triangle(i, k)))
				{
					//Get  fraction
					const T item = A_triangle(i, k) / A_triangle(k, k);

					//Process row in matrix B
					B(i, 0) = B(i, 0) - item * B(k, 0);

					for (int j = n - 1; j >= k; j--)
					{
						A_triangle(i, j) = A_triangle(i, j) - item * A_triangle(k, j);
					}
				}

				//Throw exception
				else
				{
					throw MathMatrixSingularException <Matrix <T> >("MathMatrixSingularException: ", " singular matrix, cannot perform GEM (forward); (rows_count, columns_count):  ", A_triangle);
				}
			}
		}
	}


	//Solving of the equation (Backward step)
	template <typename T>
	void gemB(const Matrix <T> &A, const Matrix <T> &B, Matrix <T> &X)
	{
		const unsigned int n = A.cols();

		for (int i = n - 1; i >= 0; i--)
		{
			T sum = 0;

			for (unsigned int j = i + 1; j < n; j++)
			{
				sum += A(i, j) * X(j, 0);
			}

			sum = B(i, 0) - sum;

			//Throw exception, singular matrix
			if (fabs(A(i, i)) < fabs(sum) * MIN_FLOAT)
			{
				throw MathMatrixSingularException <Matrix <T> >("MathMatrixSingularException: ", " singular matrix, cannot perform GEM (backward); (rows_count, columns_count):  ", A);
			}

			X(i, 0) = sum / A(i, i);
		}
	}


	//Determinant of the matrix using LU decomposition
	template <typename T>
	T det(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute determinant; (rows_count, columns_count):  ", A);
		}

		//Create LU decomposition
		Matrix <T> L(m, m);
		Matrix <T> U(m, m);
		Matrix <T> P(m, m, 0, 1);
		short sign = 1;
		lu(A, L, U, P, sign);

		//Compute determinant of the triangle matrix
		//det(A) = det(L) * det(P) * det(U) = 1 * (+-1) * det(U)
		T determ = 1.0;

		for (unsigned int i = 0; i < m; i++)
		{
			determ *= U(i, i);
		}

		//Return determinant
		return sign * determ;
	}


	//LU decomposition of the matrix
	template <typename T>
	void lu(const Matrix <T> &A, Matrix <T> &L, Matrix <T> &U, Matrix <T> &P, short &sign)
	{
		//LU decomposition of A, L = lower triangular matrix, U = upper triangular matrix, P = permutation matrix
		const unsigned int m = A.rows(), n = A.cols();

		//Set the determinant det(U) sign to 1
		sign = 1;

		//Is A rectangular matrix ?
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not perform LU decomposition; (rows_count, columns_count):  ", A);
		}

		//Create row permutation vector
		Matrix <unsigned int> PR(1, n);

		//Create scale vector
		Matrix <T> S(1, n);

		//Set diagonal items of L to 1, otherwise to 0
		//Set items of the row permutation matrix to <0; n-1>
		for (unsigned int i = 0; i < m; i++)
		{
			L(i, i) = 1.0;
			PR(0, i) = i;

			for (unsigned int j = 0; j < m; j++)
			{
				if (j != i)  L(i, j) = 0;

				P(i, j) = 0;
			}
		}

		//Initialize U = A
		U = A;

		//Find max item in each row to compute the scale vector
		for (unsigned int i = 0; i < n; i++)
		{
			T max_val = 0.0;

			for (unsigned int j = 0; j < n; j++)
			{
				if (fabs(U(i, j)) > max_val)
					max_val = fabs(U(i, j));
			}

			//Actualize scale vector
			if (max_val > MIN_FLOAT)
				S(0, i) = 1.0 / max_val;
		}

		//Start LU decomposition
		for (unsigned int j = 0; j < n; j++)
		{
			for (unsigned int i = 0; i < j; i++)
			{
				T sum = U(i, j);

				//Compute new U ( i, j ) item: multiply ith row and j-th column
				for (unsigned int k = 0; k < i; k++) sum -= U(i, k) * U(k, j);

				U(i, j) = sum;
			}

			//Initialize max_val and pivot index
			T max_val = 0.0;
			unsigned int i_pivot = n;

			//Find row that will be swapped and actualize row index
			for (unsigned int i = j; i < n; i++)
			{
				T sum = U(i, j);

				//Compute new U ( i, j ) item: multiply ith row and j-th column
				for (unsigned int k = 0; k < j; k++) sum -= U(i, k) * U(k, j);

				//Compute new U (i, j)
				U(i, j) = sum;

				//Compute index of the pivot
				const T val = S(0, i) * fabs(sum);

				if (val >= max_val)
				{
					max_val = val;
					i_pivot = i;
				}
			}

			//Perform row swaps in U,PR: j <-> i_pivot
			if ((j != i_pivot) && (i_pivot < n))
			{
				//Perform swap in U matrix
				const Matrix <T> U_temp = U(i_pivot, i_pivot, 0, n - 1);
				U(U(j, j, 0, n - 1), i_pivot, 0);
				U(U_temp, j, 0);

				//Perform swap in the row permutation matrix
				const unsigned int perm_temp = PR(0, i_pivot);
				PR(0, i_pivot) = PR(0, j);
				PR(0, j) = perm_temp;

				//Actualize also the scale vector
				S(0, i_pivot) = S(0, j);

				//Actualize sign of the determinant det(U)
				sign *= -1;
			}

			//Change diagonal item U ( j, j ) = 0 to "small" value before the devision
			if (U(j, j) == 0.0)
				U(j, j) = MIN_FLOAT;

			//Actualize U (i, j) from diagonal items
			if (j != n - 1)
			{
				const T val = 1.0 / U(j, j);

				for (unsigned int i = j + 1; i < n; i++)
					U(i, j) *= val;
			}
		}

		//Process L matrix together with U matrix
		for (unsigned int i = 0; i < n; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				L(i, j) = U(i, j);
				U(i, j) = 0.0;
			}

			//Actualize permutation matrix from the row permutation matrix
			P(i, PR(0, i)) = 1.0;
		}
	}


	//QR decomposition of the matrix
	//Based on modified Gram-Schmidt with reorthogonalization ( high accuracy ), Gander algorithm
	//Will work with rank deficient matrices
	template <typename T>
	void qr(const Matrix <T> &A, Matrix <T> &Q, Matrix <T> &R)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//First iteration : Q = A
		Q(A, 0, 0);

		//Perform QR decomposition
		for (int i = 0; i < n; i++)
		{
			//Compute norm
			T tt = 0, t = norm(Q(0, m - 1, i, i));

			bool orthogonalize = true;

			//Perform reortogonalization
			while (orthogonalize)
			{
				for (int j = 0; j <= i - 1; j++)
				{
					//Compute norm
					const Matrix <T> S = trans(Q(0, m - 1, j, j)) * Q(0, m - 1, i, i);

					//Actualize R item
					R(j, i) = R(j, i) + S(0, 0);

					Q(Q(0, m - 1, i, i) - Q(0, m - 1, j, j) * S(0, 0), 0, i);
				}

				//Compute norm
				tt = norm(Q(0, m - 1, i, i));

				//Perform reorthogonalization
				if ((tt > 10.0 * MIN_FLOAT * t) && (tt < t / 10.0))
				{
					orthogonalize = true;
					t = tt;
				}

				//Stop orthogonalization process
				else
				{
					orthogonalize = false;

					//If column linear dependent, stop orthogonalization
					if (tt < 10.0 * MIN_FLOAT  * t)
						tt = 0.0;
				}
			}

			//Compute diagonal item in R
			R(i, i) = tt;

			if (tt * MIN_FLOAT != 0.0)
				tt = 1.0 / tt;
			else
				tt = 0.0;

			//Update Q matrix
			Q(Q(0, m - 1, i, i)  * tt, 0, i);
		}
	}


	//QR decomposition of the matrix, return also a permutation matrix P
	//Based on Businger and Golub algorithm
	//Will not work with rank deficient matrices
	template <typename T>
	void qr(const Matrix <T> &A, Matrix <T> &Q, Matrix <T> &R, Matrix <unsigned int> &P)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Bad dimensions of A: rank defficient matrix
		if (m < n)
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix A, m < n.", "Can not compute QR decomposition with columns pivoting.");
		}

		//Bad dimensions of Q
		if ((Q.rows() != m) || (Q.cols() != m))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix Q(m, m).", "Can not compute QR decomposition with columns pivoting.");
		}

		//Bad dimensions of R
		if ((R.rows() != m) || (R.cols() != n))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix R (m, n).", "Can not compute QR decomposition with columns pivoting.");
		}

		//Bad dimensions of P
		if ((P.rows() != n) || (P.cols() != n))
		{
			throw BadDataException("BadDataException: Matrix A is rank deficient (m <n).", "Can not compute QR decomposition with columns pivoting, must be of m >=n.");
		}

		//First iteration: R = A, Q = E, P = E
		R = A;

		for (unsigned int i = 0; i < m; i++) 
			Q(i, i) = 1.0;

		for (unsigned int i = 0; i < n; i++)
			P(i, i) = 1.0;

		//Compute the norms
		Matrix <T> col_norms(1, n);

		for (unsigned int i = 0; i < n; i++)
		{
			const Matrix <T> col_norm = trans(R(0, m - 1, i, i)) * R(0, m - 1, i, i);
			col_norms(0, i) = col_norm(0, 0);
		}

		//Perform QR decomposition with columns pivoting
		for (unsigned int i = 0; i < n - 1; i++)
		{
			//Find max column norm
			T max_col_norm = col_norms(0, i);
			unsigned int permut_index = i;

			for (unsigned int j = i + 1; j < n; j++)
			{
				//Remeber new max and permutation index
				if (col_norms(0, j) > max_col_norm)
				{
					max_col_norm = col_norms(0, j);
					permut_index = j;
				}
			}

			//Stop, can not compute QR decomposition
			if (col_norms(0, permut_index) == 0.0)
				break;

			//Performs pivoting
			if (permut_index != i)
			{
				//Swap permutaion matrix
				const Matrix <unsigned int> TMP = P(0, n - 1, i, i);
				P(P(0, n - 1, permut_index, permut_index), 0, i);
				P(TMP, 0, permut_index);

				//Swap R matrix
				const Matrix <T> TMP2 = R(0, m - 1, i, i);
				R(R(0, m - 1, permut_index, permut_index), 0, i);
				R(TMP2, 0, permut_index);

				//Swap column norms matrix
				const T tmp = col_norms(0, i);
				col_norms(0, i) = col_norms(0, permut_index);
				col_norms(0, permut_index) = tmp;
			}

			//Compute Householder vector
			Matrix <T> H = hous(R(0, m - 1, i, i), i, m - 1);

			//Apply left Householder transformation
			R = R - H * (trans(H) * R);

			//Apply right Householder transformation
			Q = Q - (Q * H) * trans(H);

			//Norm downdate using Hadamard product
			const Matrix <T> col_norm = R(i, i, i + 1, n - 1) % R(i, i, i + 1, n - 1);
			col_norms(col_norms(0, 0, i + 1, n - 1) - col_norm, 0, i + 1);
		}
	}


	//Is matrix positive definite?
	template <typename T>
	bool posdef(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute inverse matrix; (rows_count, columns_count):  ", A);
		}

		//Test of the submatrix
		for (unsigned int i = 0; i < m; i++)
		{
			//Create submatrix
			Matrix <T> AS = A(0, i, 0, i);

			//Test determinant
			if (det(AS) <= 0)
				return false;
		}

		return true;
	}


	//Perform Cholesky decomposition of the matrix
	//Use chol (A, 'upper') or chol (A,'lower')
	template <typename T>
	Matrix <T> chol(const Matrix <T> &A, const TCholeskyFactorization f)
	{
		//Is matrix positive definite?
		const bool positive_definite = posdef(A);

		//Throw exception
		if (!positive_definite)
		{
			throw MathMatrixNotPositiveDefiniteException <Matrix <T> >("MathMatrixNotPositiveDefiniteException: ", " invalid dimension of the matrix (rectangle matrix), can not compute Cholesky decomposition:  ", A);
		}

		const unsigned int m = A.rows(), n = A.cols();

		//Create empty matrix
		Matrix <T> F(m, m);

		//Perform decomposition
		for (unsigned int i = 0; i < m; i++)
		{
			//Process parts of the columns till the actual row
			for (unsigned int k = 0; k < i + 1; k++)
			{
				//Compute sum ( F(k,j)^2)
				T sum_f2 = 0;
				for (unsigned int j = 0; j < k; j++)
				{
					sum_f2 += F(i, j) * F(k, j);
				}

				//Compute new diagonal element F(k, k )
				if (i == k)
				{
					F(i, k) = sqrt(A(i, i) - sum_f2);
				}

				//Compute element F(i, k ) bellow the diagonal
				else
				{
					F(i, k) = (A(i, k) - sum_f2) / F(k, k);
				}
			}
		}

		//Return lower triangular matrix
		if (f == lower)
		{
			return F;
		}

		//Return upper triangular matrix
		else
		{
			return trans(F);
		}
	}


	//Perform Gill-Murray modified Cholesky decomposition
	//Modified algorithm by  M. Overton, 2005
	//Decomposition: R'*R + E = A
	template <typename T>
	void gill(const Matrix <T> &A, Matrix <T> &R, Matrix <T> &E, bool &indefinite, const T max_error)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute tridiagonal matrix; (rows_count, columns_count):  ", A);
		}

		//Symmetric matrix
		if (!isequal(A, trans(A)))
		{
			throw BadDataException("BadDataException: Matrix A is not symmetric. ", "Can not compute tridiagonal matrix.");
		}

		//Matrix is apriori positive definite
		indefinite = false;

		//Get diagonal and its maximum element
		const Matrix <T> DA = diag(A);
		const T max_diag = max(abs(DA));

		//Get maximum off-diagonal element
		Matrix <T> AD = diag(diag(A));
		const T max_off_diag = max(abs(A) - AD);

		//Get threshold
		const T delta = max_error * std::max(max_diag + max_off_diag, 1.0);
		const T beta = sqrt(std::max(max_diag, std::max(max_off_diag / m, max_error)));

		//Initialize matrices
		Matrix <T> D(m, 1), L(m, m, 0, 1);

		//Perform the decomposition
		for (unsigned int i = 0; i < m; i++)
		{
			//Get diagonal element cjj
			const T dc = (i > 0 ? norm(L(i, i, 0, i - 1) * (D(0, i - 1, 0, 0) % trans(L(i, i, 0, i - 1)))) : 0);
			const T cjj = A(i, i) - dc;

			//Do not process last row of A
			if (i < m - 1)
			{
				//Prevent D(i,1) to be too small and L(:, i) too high
				Matrix <T> DC(m - i - 1, 1);

				//Compute update
				if (i > 0)
				{
					DC = L(i + 1, m - 1, 0, i - 1) * (D(0, i - 1, 0, 0) % trans(L(i, i, 0, i - 1)));
				}

				Matrix <T> Cij = A(i + 1, m - 1, i, i) - DC;
				const T theta = max(abs(Cij));

				//Compute D(j,1)
				D(i, 0) = std::max(fabs(cjj), std::max(theta * theta / (beta * beta), delta));

				//Compute L(:, i)
				L(Cij * (1.0 / D(i, 0)), i + 1, i);
			}

			//Last row
			else
			{
				D(i, 0) = std::max(fabs(cjj), delta);
			}

			//Test, if matrix is positive definite
			if (D(i, 0) > cjj)
			{
				//Matrix is indefinite
				indefinite = true;
			}
		}

		//Create output matrices
		for (unsigned int i = 0; i < m; i++)
		{
			L(L(0, m - 1, i, i) * sqrt(D(i, 0)), 0, i);
		}

		//Assign factorized matrix
		R = trans(L);

		//Assign E matrix
		E = A - trans(R) * R;
	}


	//Perform Gill-Murray modified Cholesky decomposition
	//Algorithm by  Brian Borchers, modified by Michael Zibulevsky
	//Decomposition: A + diag(e) = L*diag(d)*L'
	//Returns: v_neg'*A*v_neg < 0
	template <typename T>
	void gill(const Matrix <T> &A, Matrix <T> &L, Matrix <T> &d, Matrix <T> &e, Matrix <T> &v, const T max_error)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute tridiagonal matrix; (rows_count, columns_count):  ", A);
		}

		//Symmetric matrix
		if (!isequal(A, trans(A)))
		{
			throw BadDataException("BadDataException: Matrix A is not symmetric. ", "Can not compute tridiagonal matrix.");
		}

		//Get diagonal and its maximum element
		const T max_diag = max(diag(A));

		//Get maximum off-diagonal element
		Matrix <T> C = diag(diag(A));
		const T max_off_diag = max(A - C);

		//Get threshold
		const T nu = std::max(1.0, sqrt(m * m - 1));
		const T beta = std::max(max_diag, std::max(max_off_diag / nu, max_error));

		//Perform the decomposition
		Matrix <T> theta(m, 1);
		for (unsigned int i = 0; i < m; i++)
		{
			//Compute i-th row of L
			if (i > 0)
			{
				for (unsigned int j = 0; j <= i - 1; j++)
				{
					L(i, j) = C(i, j) / d(j, 0);
				}
			}

			//Update i-th column of C
			if (i > 0)
			{
				if (i < n - 1)
				{
					C(A(i + 1, n - 1, i, i) - C(i + 1, n - 1, 0, i - 1)  * trans(L(i, i, 0, i - 1)), i + 1, i);
				}
			}

			else
			{
				C(A(i + 1, n - 1, i, i), i + 1, i);
			}

			//Compute theta
			if (i == m - 1)
			{
				//The last row
				theta(i, 0) = 0;
			}
			else
			{
				//The common row
				theta(i, 0) = max(abs(C(i + 1, m - 1, i, i)));
			}

			//Compute d
			d(i, 0) = std::max(1.0e-16, std::max(fabs(C(i, i)), theta(i, 0) * theta(i, 0) / beta));

			//Compute e,  a "difference" from positive definite matrix
			e(i, 0) = d(i, 0) - C(i, i);

			//Compute new deigonal elements of C 
			for (unsigned int j = i + 1; j < m; j++)
			{
				C(j, j) += -C(j, i) * C(j, i) / d(i, 0);
			}
		}

		//Set diagonal elements of L to 1
		for (unsigned int i = 0; i < m; i++)
			L(i, i) = 1.0;

		//Find a descent direction v'Bv < 0
		Matrix <T> l(m, 1);

		//Find minimum
		unsigned int row = 0, col = 0;
		const T cmin = min(diag(C), row, col);

		//Create vector l, |l| = 1
		l(row, 0) = 1.0;

		//Solve L'v=l
		v = inv(trans(L)) * l;
	}

	//Inverse matrix calculation using LU decomposition
	template <typename T>
	Matrix <T> inv(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute inverse matrix; (rows_count, columns_count):  ", A);
		}

		//Find maximum
		const T max_val = max(A);

		if (max_val > MAX_DOUBLE)
			throw MathOverflowException <T>("MathOverflowException: bad scaled matrix, can not compute inverse matrix. ", "Max item > MAX_FLOAT.", max_val);

		//Create LU decomposition
		Matrix <T> L(m, m);
		Matrix <T> U(m, m);
		Matrix <T> P(m, m, 0, 1);
		short sign = 1;
		lu(A, L, U, P, sign);

		//Compute X = L^-1 (lower triangular matrix)
		Matrix <T> X(m, m);

		for (unsigned int j = 0; j < m; j++)
		{
			X(j, j) = 1.0;

			for (unsigned int i = j + 1; i < m; i++)
			{
				T sum = 0;

				for (unsigned int k = j; k <= i - 1; k++)
				{
					sum -= L(i, k) * X(k, j);
				}

				X(i, j) = sum;
			}
		}

		//Compute Y = U^-1 (upper triangular matrix)
		Matrix <T> Y(m, m);

		for (unsigned int j = 0; j < m; j++)
		{
			Y(j, j) = 1 / U(j, j);

			for (int i = j - 1; i >= 0; i--)
			{
				T sum = 0.0;

				for (unsigned int k = i + 1; k <= j; k++)
				{
					sum -= U(i, k) * Y(k, j) / U(i, i);
				}

				Y(i, j) = sum;
			}
		}

		//Compute inverse matrix A^-1 = U^-1 * L^-1 = X * Y * P
		return Y * X * P;
	}


	//Pseudo-inverse matrix calculation using double QR decomposition
	//Use  new algorithm  based on modified Goodall, Golub & Loan solution,
	//Works for rank deficient matrix
	template <typename T>
	Matrix <T> pinv(const Matrix <T> &A)
	{
		unsigned int m = A.rows(), n = A.cols();
		Matrix <T> L(std::max(m, n), std::max(m, n), 0.0, 1.0), W(m, m, 0.0, 1.0);

		//Compute pseudoinverse using MLS with double QR factorization
		return mlsqr(A, W, L);
	}


	//Solve MLS using QR factorization, weighted variant
	//Use  new algorithm  based on modified Goodall, Golub & Loan solution
	template <typename T>
	Matrix <T> mlsqr(const Matrix <T> &A, const Matrix <T> &W, const Matrix <T> &L)
	{
		unsigned int m = A.rows(), n = A.cols();
		bool transpose = false;
		Matrix <T> AT(n, m);

		//Rand deficient matrix
		if (m < n)
		{
			//Transpose matrix
			AT = trans(A);

			//Swap m, n
			const unsigned int temp = m;
			m = n; n = temp;

			//Set flag as true
			transpose = true;
		}

		//Initialize matrix M
		Matrix <T> AA = (transpose ? AT : A);

		//QR decomposition of A with a permutation matrix P
		Matrix <T> Q(m, m), R(m, n), B();
		Matrix <unsigned int> P(n, n);
		qr(AA, Q, R, P);

		//Create matrix of floats
		Matrix <T> PF = P * 1.0;

		//Set tolerance
		const T eps = std::max(m, n) * norm(A) * MAX_FLOAT_OPER_ERROR;

		//Find index of the first diagonal element < eps
		int k = 0;

		for (; (k < n) && (fabs(R(k, k)) >= eps); k++);

		//Correct matrix (not rank deficient matrix)
		if ((k == n) && (fabs(R(k - 1, k - 1)) >= eps))
		{
			//Compute MLS problem (L = E, pseudoinverse)
			Matrix <T> AA_I = PF * inv(trans(R) * trans(Q) * W * Q * R) * trans(R) * trans(Q) * W * L;

			return AA_I;
		}

		//Rank deficient matrix
		else
		{
			//Get R submatrix containing only "large" elements
			Matrix <T> RS = R(0, k - 1, 0, n - 1);

			//Get Q submatrix
			Matrix <T> QS = Q(0, m - 1, 0, k - 1);

			//Compute second QR factorization from remaining elements
			Matrix <T> Q2(n, k), R2(k, k);
			qr(trans(RS), Q2, R2);

			//Trim Q to submatrix
			Matrix <T> Q_S = Q(0, m - 1, 0, k - 1);

			//Compute MLS problem (L = E, pseudoinverse)
			Matrix <T> AA_I = PF * Q2 * inv(trans(QS) * W * QS * trans(R2)) * trans(Q_S) * W * L;

			//Return A_I or transposed A_I
			return (transpose ? trans(AA_I) : AA_I);
		}
	}


	//Compute eigevalues L and eigenvectors V of symmetric tridiagonal matrix
	//Computation based on the QR algorithm with the Wilkinson shift
	//Eigenvalues and eigenvecors are sorted in the ascending order
	template <typename T>
	void eig(const Matrix <T> &A, Matrix<T> &V, Matrix<T> &L, const unsigned int max_iter, const T tolerance)
	{
		//Modified algorithm by Brian Brodie
		unsigned int m = A.rows(), n = A.cols();

		//Not rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute eigenvalues; (rows_count, columns_count):  ", A);
		}

		/*
		//Not positive definite matrix
		const bool positive_definite = posdef(A);

		//Throw exception
		if (!positive_definite)
		{
		throw MathMatrixNotPositiveDefiniteException <Matrix <T> >("MathMatrixNotPositiveDefiniteException: ", " invalid dimension of the matrix (rectangle matrix), can not compute eigenvalue decomposition:  ", A);
		}
		*/

		//Create tridiagonal matrix
		Matrix <T> LV(m, 1, 1), RV(m, 1, 1);
		const Matrix <T> AT = tridiag(A, LV, RV);

		//Get diagonal and upper diagonal
		Matrix <T> B = diag(AT);
		Matrix <T> CT = diag(AT, -1);

		//Get size of the diagonal
		unsigned int mb = B.rows(), nb = B.cols();

		//Modify C matrix: first element twice
		Matrix <T> C(mb, 1);
		C(0, 0) = CT(0, 0);
		C(CT(0, mb - 2, 0, 0), 1, 0);

		//Create additional variables
		short index = mb - 1;
		T w_shift = 0;
		Matrix <T> S(mb, 1), D(mb, 1);

		//Perform iterations
		for (unsigned int iter = 0; (iter < max_iter) && (index != 0); iter++)
		{
			//Compute trace and determinant
			const T trace = B(index - 1, 0) + B(index, 0);
			const T deter = B(index - 1, 0) * B(index, 0) - C(index, 0) * C(index, 0);

			//Solve quadratic equation
			const T discriminant = sqrt(trace * trace - 4 * deter);
			const T m1 = (trace + discriminant) / 2;
			const T m2 = (trace - discriminant) / 2;

			//Test roots to diagonal elements and actualize Wilkinson shift
			if (fabs(m1 - B(index, 0)) < fabs(m2 - B(index, 0)))
				S(0, 0) = m1;
			else
				S(0, 0) = m2;

			w_shift += S(0, 0);

			//Actualize diagonal elements
			for (unsigned int i = 0; i <= index; i++)
				B(i, 0) -= S(0, 0);

			T c_old = C(1, 0);

			//QR algorithm
			for (unsigned int i = 1; i <= index; i++)
			{
				int j = i - 1;

				//Compute norm
				const T r = sqrt(B(j, 0) * B(j, 0) + c_old * c_old);
				D(i, 0) = B(j, 0) / r;
				S(i, 0) = c_old / r;
				B(j, 0) = r;

				//Actualize B, C: diagonal and off diagonal
				const T t1 = D(i, 0) * C(i, 0) + S(i, 0) * B(i, 0);
				const T t2 = -S(i, 0) * C(i, 0) + D(i, 0) * B(i, 0);
				C(i, 0) = t1;
				B(i, 0) = t2;

				if (i != index)
				{
					c_old = C(i + 1, 0);
					C(i + 1, 0) *= D(i, 0);
				}
			}

			//Actualize first diagonal and second off diagonal elements
			B(0, 0) = D(1, 0) * B(0, 0) + S(1, 0) * C(1, 0);
			C(1, 0) = S(1, 0) * B(1, 0);

			//Actualize remaining diagonal and off diagonal elements
			for (unsigned int i = 1; i <= index - 1; i++)
			{
				B(i, 0) = S(i + 1, 0) * C(i + 1, 0) + D(i, 0) * D(i + 1, 0) * B(i, 0);
				C(i + 1, 0) = S(i + 1, 0) * B(i + 1, 0);
			}

			B(index, 0) *= D(index, 0);

			//B.print();
			//C.print();
			//D.print();
			//S.print();

			//Compute eigenvectors
			for (unsigned int i = 1; i <= index; i++)
			{
				Matrix <T> CO = V(0, mb - 1, i - 1, i - 1) * D(i, 0) + V(0, mb - 1, i, i) * S(i, 0);
				V(V(0, mb - 1, i - 1, i - 1) * (-S(i, 0)) + V(0, mb - 1, i, i) * D(i, 0), 0, i);
				V(CO, 0, i - 1);
			}

			//Sufficient accuracy, compute the eigenvalue
			if (fabs(C(index, 0)) < tolerance)
			{
				L(index, 0) = B(index, 0) + w_shift;
				index--;
			}
		}

		//Actualize first eigenvalue
		L(0, 0) = B(0, 0) + w_shift;

		//Sort eigenvalues
		Matrix <unsigned int> IX(m, 1);
		sortrows(L, IX, 0);

		//Sort eigenvectors accoridng to the eigenvalues
		Matrix <T> VP = V;
		for (unsigned int i = 0; i < m; i++)
		{
			//Copy row to permutated matrix
			VP(V(IX(i, 0), IX(i, 0), 0, m - 1), i, 0);
		}

		//Assign permutated matrix
		V = VP;
	}


	//Pseudo-inverse matrix calculation using SVD algorithm
	//Works also for rank deficient matrices
	template <typename T>
	Matrix <T> pinvs(const Matrix <T> &A, const T tolerance)
	{
		unsigned int m = A.rows(), n = A.cols();
		bool transpose = false;
		const T MAX_SVD_ITER = 1000;
		Matrix <T> AT(n, m);

		//Rank deficient matrix
		if (m < n)
		{
			//Transpose matrix
			AT = trans(A);

			//Swap m, n
			const unsigned int temp = m;
			m = n; n = temp;

			//Set flag as true
			transpose = true;
		}

		//Initialize matrix M
		Matrix <T> M = (transpose ? AT : A);

		//Create matrices
		Matrix <T> U(m, m), B(m, n), V(n, n);

		//Compute SVD decomposition
		svd(M, U, B, V, MAX_SVD_ITER);

		//Get B diagonal
		Matrix <T> B_D = diag(B);

		//Compute B inverse
		Matrix <T> B_I(n, n);

		//Compute threshold
		const T treshold = std::max(m, n) * norm(M) * MAX_FLOAT_OPER_ERROR;

		//Compute inverse matrix
		for (unsigned int i = 0; i < n; i++)
		{
			//Inverse diagonal elements
			if (fabs(B_D(i, 0)) > treshold)
				B_I(i, i) = 1.0 / B_D(i, 0);

			//Reset small diagonal elements to zero
			else B_I(i, i) = 0.0;
		}

		//Compute pseudoinverse
		Matrix <T> A_I = V *  B_I * trans(U(0, m - 1, 0, n - 1));

		//Return A_I or transposed A_I
		return (transpose ? trans(A_I) : A_I);
	}


	//Pseudo-inverse matrix calculation using double QR decomposition
	//More numerically stable than Moore-Penrose algorithm
	//Use Goodall, Golub & Loan algorithm, works for rank deficient matrix
	template <typename T>
	Matrix <T> pinv1(const Matrix <T> &A)
	{
		unsigned int m = A.rows(), n = A.cols();
		bool transpose = false;
		Matrix <T> AT(n, m);

		//Rand deficient matrix
		if (m < n)
		{
			//Transpose matrix
			AT = trans(A);

			//Swap m, n
			const unsigned int temp = m;
			m = n; n = temp;

			//Set flag as true
			transpose = true;
		}

		//Initialize matrix M
		Matrix <T> AA = (transpose ? AT : A);

		//QR decomposition of A with a permutation matrix E
		Matrix <T> Q(m, m), R(m, n), B();
		Matrix <unsigned int> P(n, n);
		qr(AA, Q, R, P);

		//Permutation matrix represented by the row matrix
		Matrix <T> PR(1, n);
		Matrix <unsigned int> M(1, n);

		for (unsigned int i = 0; i < n; i++)
			M(0, i) = i;

		//Convert E to row matrix PR
		PR = M  * P;

		//Reversed row permutation matrix
		Matrix <unsigned int> PR_R(1, n);

		for (unsigned int i = 0; i < n; i++)
			PR_R(0, PR(0, i)) = i;

		//Set tolerance
		const T eps = std::max(m, n) * norm(A) * MAX_FLOAT_OPER_ERROR;

		//Find index of the first diagonal element < eps
		int k = 0;

		for (; k < n; k++)
		{
			if (fabs(R(k, k)) < eps)
				break;
		}

		//Get R submatrix
		Matrix <T> RS = R(0, k - 1, 0, n - 1);

		//Get transposed Q submatrix
		Matrix <T> QTS = trans(Q(0, m - 1, 0, k - 1));

		//Correct matrix (not rank deficient matrix)
		if (k == n)
		{
			//Compute reversed permutation matrix ER (n,n) from reversed row permutation matrix
			Matrix <T> P_R(n, n);

			for (unsigned int i = 0; i < n; i++)
			{
				P_R(i, PR_R(0, i)) = 1.0;
			}

			//Compute inverse matrix
			Matrix <T> RI = inv(RS) * QTS;

			//Compute inverse matrix using the reversed permutation
			return P_R * RI;
		}

		//Rank deficient matrix
		else
		{
			//Compute second QR factorization
			Matrix <T> Q2(n, k), R2(k, k);

			qr(trans(RS), Q2, R2);

			//Compute reversed permutation matrix PR (n,n) from reversed row permutation matrix
			Matrix <T> P_R(n, k);

			for (unsigned int i = 0; i < n; i++)
			{
				for (int j = 0; j < k; j++)
				{
					P_R(i, j) = Q2(PR_R(0, i), j);
				}
			}

			//Transpose R2 matrix
			Matrix <T> R2T = trans(R2);

			//Compute inverse: from squared matrices
			Matrix <T> R2I = inv(R2T) * QTS;

			//Compute pseudoinverse
			Matrix <T> AA_I = P_R * R2I;

			//Return A_I or transposed A_I
			return (transpose ? trans(AA_I) : AA_I);
		}
	}


	//Pseudo-inverse matrix calculation using Moore - Penrose inverse
	//Algorithm by Pierre Courrieu, 2005 with the Cholesky factorization
	//Fast but inappropriate for large numbers
	template <typename T>
	Matrix <T> pinv2(const Matrix <T> &A)
	{
		unsigned int m = A.rows(), n = A.cols();
		bool transpose = false;

		//Compute A*A
		Matrix <T> A2(std::max(m, n), std::min(m, n));

		//Input matrix must satisfy m > n, otherwise transpose
		if (m < n)
		{
			A2 = A * trans(A);
			transpose = true;
			n = m;
		}

		else A2 = trans(A) * A;

		//Find minimum positive value in diagonal
		T min_val = MAX_FLOAT;

		for (unsigned int i = 0; i < n; i++)
			if ((A2(i, i) > 0) && (A2(i, i) < min_val)) min_val = A2(i, i);

		//Compute the tolerance from this value
		Matrix <T> L(n, n);
		const T eps = min_val * 1.0e-20;

		//Compute L using the Cholesky factorization
		int j = -1;

		for (unsigned int i = 0; i < n; i++)
		{
			//Increment index
			j++;

			//Set new matrix L
			L((j == 0 ? A2(i, n - 1, i, i) : A2(i, n - 1, i, i) - L(i, n - 1, 0, j - 1) * trans(L(i, i, 0, j - 1))), i, j);

			//Iterative computation of L matrix using Cholesky algorithm
			if (L(i, j) > eps)
			{
				L(i, j) = sqrt(L(i, j));

				//Actualize L submatrix
				if (i < n - 1)
				{
					L(L(i + 1, n - 1, j, j) / L(i, j), i + 1, j);
				}
			}

			//Decrement j
			else j--;
		}

		//Get submatrix: first j columns of L
		const Matrix <T> L2 = L(0, n - 1, 0, j);

		//Compute inverse matrix M = inv ( L2'*L2 )
		const Matrix <T> M = inv(trans(L2) * L2);

		//Compute pseudoinverse matrix using Moore - Penrose algorithm: I = A'*L2*M*M*L2'
		if (transpose)
			return trans(A) * L2 * M * M * trans(L2);

		//Or I = L2*M*M*L2'A'
		return L2 * M * M * trans(L2) * trans(A);
	}


	//Compute SVD using Demel-Kahan algorithm,
	//Accurate singular values of bidiagonal matrices, 1990, pp. 11
	//Bidigonal matrix A is iteratively processed using the QR zero shift algorithm with Givens rotations
	template <typename T>
	void svd(const Matrix <T> &A, Matrix <T> &U, Matrix <T> &B, Matrix <T> &V, unsigned int max_iterations)
	{
		unsigned int m = A.rows(), n = A.cols();

		//Bad dimensions of U
		if ((U.rows() != m) || (U.cols() != m))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix A.", "Can not compute SVD.");
		}

		//Bad dimensions of B
		if ((B.rows() != m) || (B.cols() != n))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix B.", "Can not compute SVD.");
		}

		//Bad dimensions of V
		if ((V.rows() != n) || (V.cols() != n))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix V.", "Can not compute SVD.");
		}

		bool transpose = false;
		Matrix <T> AT(n, m);

		//Rank deficient matrix
		if (m < n)
		{
			//Transpose matrix
			AT = trans(A);

			//Swap m, n
			const unsigned int temp = m;
			m = n; n = temp;

			//Set flag as true
			transpose = true;

			//Set new size for U, B, V matrices
			U.res(m, m); B.res(m, n); V.res(n, n);
		}

		//Initialize matrix M
		Matrix <T> M = (transpose ? AT : A);

		//Perform A reduction to the bidiagonal form using Householder rotations
		bidiag(M, U, B, V);

		//Get diagonal and super diagonal of the B matrix
		Matrix <T> D = diag(B);
		Matrix <T> E = diag(B, 1);

		//Set convergence parameters
		const T tolerance = 100.0 * MAX_FLOAT_OPER_ERROR;
		max_iterations = std::max(max_iterations, 500 * n * n);

		//Process diagonal
		Matrix <T> LA(n, 1);
		LA(n - 1, 0) = fabs(D(n - 1, 0));

		for (int i = n - 2; i >= 0; i--)
		{
			const T denom = LA(i + 1, 0) + fabs(E(i, 0));
			LA(i, 0) = (denom != 0.0 ? fabs(D(i, 0)) * LA(i + 1, 0) / denom : 0.0);
		}

		//Process superdiagonal
		Matrix <T> MU(n, 1);
		MU(0, 0) = fabs(D(0, 0));

		for (unsigned int i = 0; i < n - 1; i++)
		{
			const T denom = MU(i, 0) + fabs(E(i, 0));
			MU(i + 1, 0) = (denom != 0.0 ? fabs(D(i + 1, 0)) * MU(i, 0) / denom : 0.0);
		}

		//Get smallest singular value
		const T sigma_min = std::min(min(LA), min(MU));

		//Compute treshold
		const T treshold = 0.01 * std::max(sigma_min * tolerance, MIN_FLOAT * max_iterations);

		//Perform iterations using QR zero shift algorithm
		unsigned int iterations = 0, i1 = 0, i2 = n - 2;

		while (iterations < max_iterations)
		{
			//Find first super diagonal item abs(E) > treshold, start from right bottom
			while ((fabs(E(i1, 0)) < treshold) && (i1 < i2)) i1++;

			//Find first super diagonal item abs(E) > treshold, start from left top
			while ((fabs(E(i2, 0)) < treshold) && (i2 > i1)) i2--;

			//All E elements on diagonal are < treshold
			if ((i1 == i2) && (fabs(E(i1, 0)) <= treshold))
			{
				break;
			}

			//Get part of matrix need to be changed
			Matrix <T> D_S = D(i1, i2 + 1, 0, 0);
			Matrix <T> E_S = E(i1, i2, 0, 0);
			Matrix <T> U_S = U(0, m - 1, i1, i2 + 1);
			Matrix <T> V_S = V(0, n - 1, i1, i2 + 1);

			//Perform next iteration
			qr0Shift(D_S, E_S, U_S, V_S);

			//Actualize matrices
			D(D_S, i1, 0);
			E(E_S, i1, 0);
			U(U_S, 0, i1);
			V(V_S, 0, i1);

			//Increment iterations
			iterations++;
		}

		//Not transposed matrix A
		if (!transpose)
		{
			//Copy D elements to the B
			for (unsigned int i = 0; i < n; i++) B(i, i) = D(i, 0);

			return;
		}

		//Transposed matrix A: change dimensions
		B.res(n, m);

		for (unsigned int i = 0; i < n; i++) B(i, i) = D(i, 0);

		//Remember U,V
		Matrix <T> TMP1 = U; Matrix <T> TMP2 = V;

		//Resize U, V
		U.res(n, n); V.res(m, m);

		//Switch U <-> V
		U = TMP2; V = TMP1;
	}


	//Convert matrix to bidiagonal form using Householder rotations
	template <typename T>
	void bidiag(const Matrix <T> &A, Matrix <T> &U, Matrix <T> &B, Matrix <T> &V)
	{
		unsigned int m = A.rows(), n = A.cols(), m_old = m;

		//Bad dimensions of U
		if ((U.rows() != m) || (U.cols() != m))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix A.", "Can not compute bidiagonalization.");
		}

		//Bad dimensions of B
		if ((B.rows() != m) || (B.cols() != n))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix B.", "Can not compute bidiagonalization.");
		}

		//Bad dimensions of V
		if ((V.rows() != n) || (V.cols() != n))
		{
			throw BadDataException("BadDataException: invalid dimension of the matrix V.", "Can not compute bidiagonalization.");
		}

		//Matrix is not rank deficient
		bool rank_def = false;
		Matrix <T> AE = A;

		//Rank deficient matrix
		if (m < n)
		{
			//Set flag as true
			rank_def = true;

			//Set dimensions
			m = n;

			//Set new size for A, U, B, matrices
			AE.res(n, n); U.res(n, n); B.res(n, n);
		}

		//Initialize matrices
		B = AE;
		U = Matrix <T>(m, m, 0, 1);
		V = Matrix <T>(n, n, 0, 1);

		//Set amount of rotations: n-1 for squared matrix, n for rectangular matrix
		unsigned int mn = (m > n ? n : n - 1);

		//Perform bidiagonalization using Householder transformation
		for (int i = 0; i < mn; i++)
		{
			//Eliminate all non-zero elements below the diagonal using Householder transformation
			Matrix <T> H1 = hous(B(0, m - 1, i, i), i);

			//Left multiplication
			B = H1 * B;
			U = U * H1;

			//Eliminate all non-zero elements to the right of the superdiagonal Householder transformation
			if (i < n - 2)
			{
				Matrix <T> H2 = hous(trans(B(i, i, 0, n - 1)), i + 1);

				//Right multiplication
				B = B * trans(H2);
				V = H2 * V;
			}
		}

		//Transpose V
		V = trans(V);

		//Rank deficient matrix, trim U, B matrices to m_old rows
		if (rank_def)
		{
			U.res(m_old, m_old);
			B.res(m_old, n);
		}
	}


	//Perform QR zero shift using givens transformation, Demel & Kahan,
	//Accurate singular values of bidiagonal matrices, 1990, pp. 14
	template <typename T>
	void qr0Shift(Matrix <T> &D, Matrix <T> &E, Matrix <T> &U, Matrix <T> &V)
	{
	
		const unsigned int m = D.rows(), n1 = U.rows(), n2 = V.rows();
		T cs = 1.0, oldcs = 1.0, r = 0.0, sn = 0.0, oldsn = 0.0;

		for (unsigned int i = 0; i < m - 1; i++)
		{
			//Compute givens rotation
			givens(cs * D(i, 0), E(i, 0), cs, sn, r);

			//Get submatrices
			Matrix <T> V_S = V(0, n2 - 1, i, i);
			Matrix <T> V_S1 = V(0, n2 - 1, i + 1, i + 1);

			//Update items of V matrix applying Givens rotations
			updgivens(cs, sn, V_S, V_S1);

			//Update U matrices
			V(V_S, 0, i);
			V(V_S1, 0, i + 1);

			if (i > 0)
				E(i - 1, 0) = r * oldsn;

			//Compute givens rotation
			givens(oldcs * r, D(i + 1, 0) * sn, oldcs, oldsn, D(i, 0));

			//Get submatrices
			Matrix <T> U_S = U(0, n1 - 1, i, i);
			Matrix <T> U_S1 = U(0, n1 - 1, i + 1, i + 1);

			//Update items of V matrix applying Givens rotations
			updgivens(oldcs, oldsn, U_S, U_S1);

			//Update U matrices
			U(U_S, 0, i);
			U(U_S1, 0, i + 1);
		}

		//Actualize diagonal and superdiagonal
		const T h = cs * D(m - 1, 0);
		E(m - 2, 0) = h * oldsn;
		D(m - 1, 0) = h * oldcs;
	}


	//Compute Householder vector for column vector given start and end indices i, j
	template <typename T>
	Matrix <T> hous(const Matrix <T> &A, const unsigned int i, const unsigned int j)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Not a vector, n > 1
		if (n > 1)
		{
			throw BadDataException("BadDataException: Matrix A is not a column vector. ", "Can not compute the Householder vector.");
		}

		//Bad index condition: i > j
		if (i > j)
		{
			throw BadDataException("BadDataException: index i > j. ", "Can bot compute the Householder vector.");

		}

		//Bad index number: j > m - 1
		if (j > m - 1)
		{
			throw BadDataException("BadDataException: index i > n - 1 . ", "Can bot compute the Householder vector.");
		}

		//Initialize Householder vector
		Matrix <T> H(m, n);
		H(A(i, j, 0, 0), i, 0);

		//Compute Householder vector
		H(i, 0) = H(i, 0) + (A(i, 0) >= 0 ? norm(A(i, j, 0, 0)) : -norm(A(i, j, 0, 0)));

		//Normalize Householder vector
		const Matrix <T>  norm2 = trans(H) * H;

		if (norm2(0, 0) > 0)
			H = H * sqrt(2.0 / (norm2(0, 0)));

		return H;
	}


	//Compute Householder vector H from column vector A so as all elements of H * b bellow k-index are set to zero
	template <typename T>
	Matrix <T> hous(const Matrix <T> &A, const unsigned int i)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Not a vector, n > 1
		if (n > 1)
		{
			throw BadDataException("BadDataException: Matrix A is not a column vector. ", "Can not compute the Householder vector.");
		}

		//Bad index number: j > m - 1
		if (i > m - 1)
		{
			throw BadDataException("BadDataException: index i > n - 1 . ", "Can not compute the Householder vector.");
		}

		//Initialize Householder vector
		Matrix <T> H(m, m);

		//Matrix D
		Matrix <T> D = A(i, m - 1, 0, 0);

		//Compute norm
		const T d_norm = (D(0, 0) >= 0 ? -norm(D) : norm(D));

		//Zero vector in input: stop computation
		if (fabs(d_norm) < MIN_FLOAT)
		{
			H = Matrix <T>(m, m, 0.0, 1.0);
			return H;
		}

		//Create V matrix and initialize
		Matrix <T> V(m - i, 1);
		V(0, 0) = sqrt(0.5 * (1 - D(0, 0) / d_norm));

		//Actualize V
		const T  p = -d_norm * V(0, 0);
		V(D(1, m - i - 1, 0, 0) / (2 * p), 1, 0);

		//Create W matrix
		Matrix <T> W(m, 1);
		W(V, i, 0);

		//Compute H
		Matrix <T> E(m, m, 0.0, 1.0);
		H = E - W * trans(W) * 2.0;

		return H;
	}


	//Compute Givens rotation
	template <typename T>
	void givens(const T f, const T g, T &c, T &s, T &r)
	{
		//   |  c s | | a | = | r |
		//   | -s c | | b |   | 0 |

		if (f == 0.0)
		{
			c = 0.0;
			s = 1.0;
			r = g;
		}

		else if (fabs(f) > fabs(g))
		{
			const T t = g / f;
			const T t1 = sqrt(1 + t * t);
			c = 1.0 / t1;
			s = t * c;
			r = f * t1;
		}

		else
		{
			const T t = f / g;
			const T t1 = sqrt(1 + t * t);
			s = 1.0 / t1;
			c = t * s;
			r = g * t1;
		}
	}


	//Update U, V matrix using Givens rotation results Demel & Kahan,
	//Accurate singular values of bidiagonal matrices, 1990, pp. 14
	template <typename T>
	void updgivens(const T cs, const T sn, Matrix <T> & V1, Matrix <T> &V2)
	{
		const unsigned int m = V1.rows(), n = V1.cols();

		for (unsigned int i = 0; i < m; i++)
		{
			const T t = V1(i, 0);
			V1(i, 0) = cs * t + sn * V2(i, 0);
			V2(i, 0) = -sn * t + cs * V2(i, 0);
		}
	}


	//Load matrix from file
	template <typename T>
	Matrix <T> load(const char *file, const bool print_exception, std::ostream * output)
	{
		try
		{
			//Load file
			const TVector2D <std::string> file_content = File::loadFileToWords(file);

			//Get rows count
			const unsigned int m = file_content.size();

			//At least one row was load from file
			if (m > 0)
			{
				//Get columns count
				const unsigned n = file_content[0].size();

				//Create matrix
				Matrix <T> A(m, n);

				//Add items into matrix
				for (unsigned int i = 0; i < m; i++)
				{
					for (unsigned int j = 0; j < n; j++)
					{
						A(i, j) = atof(file_content[i][j].c_str());
					}
				}

				return A;
			}

			//Throw exception
			else throw BadDataException("BadDataException: no items in matrix. ", "Can not load matrix from a file.");
		}

		//Some error during points processing has appeared
		catch (Exception & error)
		{
			//Print error and do not add projection to the list
			if (print_exception)
			{
				error.printException();
			}

			//Throw exception
			throw;
		}
	}

}

#endif
