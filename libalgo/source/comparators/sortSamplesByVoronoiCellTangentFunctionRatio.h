// Description: Sort cartographic samples by Voronoi cell tangent function differences

// Copyright (c) 2010 - 2011
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

#ifndef sortSamplesByVoronoiCellTangentFunctionRatio_H
#define sortSamplesByVoronoiCellTangentFunctionRatio_H

//Forward declaration
template <class T>
class Sample;

//Sorter by Voronoi cell tangent function difference ratio
class sortSamplesByVoronoiCellTangentFunctionRatio
{
        public:
                template <typename T>
                bool operator() ( const Sample <T> &s1, const Sample <T> &s2 ) const
                {
                        return s1.getVoronoiCellTangentFunctionRatio() < s2.getVoronoiCellTangentFunctionRatio();
                }
};

#endif

