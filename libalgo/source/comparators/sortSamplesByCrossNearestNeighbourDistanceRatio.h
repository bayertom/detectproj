// Description: Sort cartographic samples by cross nearest neighbour distance ratio
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

#ifndef sortSamplesByCrossNearestNeighbourDistanceRatio_H
#define sortSamplesByCrossNearestNeighbourDistanceRatio_H

//Forward declaration
template <class T>
class Sample;


//Sorter by cross nearest neighbour distance ratio
class sortSamplesByCrossNearestNeighbourDistanceRatio
{
        public:
                template <typename T>
                bool operator() ( const Sample <T> &s1, const Sample <T> &s2 ) const
                {
                        return s1.getCrossNearestNeighbourDistanceRatio() < s2.getCrossNearestNeighbourDistanceRatio();
                }
};


#endif

