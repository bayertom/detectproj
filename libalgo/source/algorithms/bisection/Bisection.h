// Description: 1D optimization based on the bisection

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

#ifndef BISECTION_H
#define BISECTION_H

#include "libalgo/source/structures/matrix/Matrix.h"

// 1D optimization based on the bisection
class Bisection
{
	public:

		template <typename Function, typename T>
		static void bisection(Function function, Matrix <T> &A ,Matrix <T> &B, const T eps, const T max_diff, Matrix <T> &XMIN, T &fmin, unsigned short &iterations, const unsigned short max_iterations );

};

#include "Bisection.hpp"

#endif