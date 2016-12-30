// Description: Create random permutation of n-elements od vector

// Copyright (c) 2010 - 2014
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

#ifndef RandomPermutation_H
#define RandomPermutation_H


#include "libalgo/source/structures/matrix/Matrix.h"


class RandomPermutation
{
        public:
		static Matrix <unsigned int> randperm(const unsigned int n, const unsigned int k);

};

#include "RandomPermutation.hpp"

#endif

