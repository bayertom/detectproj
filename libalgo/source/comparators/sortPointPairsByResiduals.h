// Description: Sort point pairs by standard deviation

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


#ifndef sortPointPairsByDeviation_H
#define sortPointPairsByDeviation_H

//#include "libalgo/source/algorithms/transformation/Transformation2D.h"

template <typename T>
struct TDevIndexPair;

//Sort point pairs by standard deviation
template <typename T>
class sortPointPairsByResiduals
{
        public:

                bool operator() ( const typename TDevIndexPair <T>::Type & p1, const typename TDevIndexPair <T>::Type & p2 ) const
                {
                        return  p1.first < p2.first;
                }
};

#endif
