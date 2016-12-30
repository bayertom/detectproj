// Description: Length of the polyline

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


#ifndef PolyLineLength_HPP
#define PolyLineLength_HPP

#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"


template <typename Point>
typename Point::Type PolyLineLength::getPolyLineLength ( const Container <Point> &points )
{
        //Get length of the polyline
        const unsigned int points_size = points.size();
        typename Point::Type length = 0;

        for ( unsigned int i = 0; i < points_size - 1; i++ )
        {
                length += EuclDistance::getEuclDistance2D ( points [i].getX(), points [i].getY(),
                                points [i + 1].getX(), points [i + 1].getY() );
        }

        return length;
}


template <typename T>
T PolyLineLength::getPolyLineLength ( const Container <Node3DCartesian <T> *> &points, const TIndexList & il )
{
        //Get length of the polyline
        const unsigned int n = il.size();
        const unsigned int points_size = points.size();
        T length = 0;

        for ( unsigned int i = 0; i < n - 1; i++ )
        {
                length += EuclDistance::getEuclDistance2D ( points [il[i]], points [il[i + 1]] );
        }

        return length;
}

#endif
