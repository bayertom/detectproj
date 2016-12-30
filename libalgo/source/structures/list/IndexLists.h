// Description: Lists of indices

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


#ifndef IndexLists_H
#define IndexLists_H

#include <vector>
#include <set>

//Forward declaration
template <typename T>
class sortPointsIndicesByX;


//Unsorted list of indices
class TIndexList: public std::vector <unsigned int>
{
        public:
                TIndexList () : std::vector <unsigned int> ( 0 ) {}
                TIndexList ( const unsigned int n ) : std::vector <unsigned int> ( n ) {}
};


//Set of indices (points indices) sorted by X coordinates and subsequently by Y coordinates
template <typename Point>
struct TIndexSet
{
        typedef std::set <unsigned int, sortPointsIndicesByX <Point> >  Type;
};


#endif
