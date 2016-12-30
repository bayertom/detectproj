// Description: Test if meridian / parallel point is present in TDevIndexPairs
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


#ifndef isEqualMeridianParallelPointIndex_H
#define isEqualMeridianParallelPointIndex_H

#include "libalgo/source/algorithms/transformation/Transformation2D.h"


//Test if meridian / parallel point is present in TDevIndexPairs
template <typename T>
class isEqualMeridianParallelPointIndex
{
        private:
                unsigned int index;

        public:
                isEqualMeridianParallelPointIndex ( const unsigned int index_ ) : index ( index_ ) {}

                bool operator () ( const typename TDevIndexPair <T>::Type &pair )
                {
                        return index == pair.second;
                }
};

#endif
