// Description: Sort fitting lines computed by RANSAC using hash

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

#ifndef sortFittingLinesByHash_H
#define sortFittingLinesByHash_H

#include "libalgo/source/algorithms/leastsquaresfitting/LeastSquaresFitting.h"

//Sort lines by the hash
template <typename Line>
class sortFittingLinesByHash
{
        public:
                bool operator() ( const Line & l1, const Line & l2 ) const
                {
                        return l1.hash_val <  l2.hash_val;
                }
};

#endif
