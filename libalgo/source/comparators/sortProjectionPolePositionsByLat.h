// Description: Sort all generated map projection positions by latitude
// Only appropriate positions will be further analyzed
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


#ifndef sortProjectionPolePositionsByLat_H
#define sortProjectionPolePositionsByLat_H


#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


class sortProjectionPolePositionsByLat
{
        public:
                template <typename T>
                bool operator() ( const TProjectionPolePosition <T> & pa1, const TProjectionPolePosition <T> & pa2 ) const
                {
                        return ( ( pa1.cart_pole.getLat() < pa2.cart_pole.getLat() ) || ( pa1.cart_pole.getLat() == pa2.cart_pole.getLat() ) && ( pa1.cart_pole.getLon() < pa2.cart_pole.getLon() ) );
                }
};


#endif
