// Description: 1D optimization based on the bisection

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


#ifndef BISECTION_HPP
#define BISECTION_HPP

template <typename Function, typename T>
void Bisection::bisection(Function function, Matrix <T> &A, Matrix <T> &B, const T eps, const T max_diff, Matrix <T> &XMIN, T &fmin, unsigned short &iterations, const unsigned short max_iterations )
{
	//Performs 1D optimization using the bisection
	T fc = 0;
	iterations = 0;

	//Half of the interval
	Matrix <T> C = 0.5 * (A + B);

	//Perform iterations
	while (norm(B - A) > eps)
	{
		//Half of the subinterval
		Matrix <T>  D = 0.5 * (A + C);
		Matrix <T>  E = 0.5 * (B + C);

		fc = function(C);
		const T fd = function(D);
		const T fe = function(E);

		//Mid point of the left sub interval has better cost than midpoint of the interval
		if (fd < fc)
		{
			B = C;
			C = D;
		}

		//Midpoint of the interval has better cost than bost mid points of sub intervals
		else if ((fc < fd) && (fc < fe))
		{
			A = D;
			B = E;
		}

		//Mid point of the right sub interval has better cost than midpoint of the interval
		else if (fe < fc)
		{
			A = C;
			C = E;
		}

		//Too many iterations or the same function value inside the interval
		if ((iterations > max_iterations)/*  (fabs(function(A) - function(B)) < max_diff) && (fabs(function(c) - function(b)) < max_diff )*/)
		break;

		//std::cout << fc << "  " << A(4, 0) << "  " << B(4, 0) << "  " << C(4, 0) << "  " << D(4, 0) << "  " << E(4, 0) << '\n';

		iterations++;
	}

	//Assign minimum
	XMIN = C;
	fmin = fc;
	
	//std::cout << " xmin = " << C(4, 0) << " fmin = " << fc << '\n';
}

#endif

