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


#include "RandomPermutation.h"

#include <vector>  
#include <algorithm>


Matrix <unsigned int> RandomPermutation::randperm(const unsigned int n, const unsigned int k)
{
	//Create random permutation of indices
	std::vector <unsigned int> indices;

	//Store indices
	for (unsigned int i = 1; i < n; i++)
	{
		indices.push_back(i);
	}

	//Use random number generator
	std::random_shuffle(indices.begin(), indices.end());

	//Copy selected k indices
	Matrix <unsigned int> I(1, k);
	for (unsigned int i = 0; i < k; i++)
	{
		I(0, i) = indices[i];
	}

	return I;
}