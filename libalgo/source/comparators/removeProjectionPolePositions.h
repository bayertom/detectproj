// Description: Remove projection pole positions by complex criteria values
// Only appropriate map projection positions will be further analyzed
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


#ifndef removeProjectionPolePositions_H
#define removeProjectionPolePositions_H


#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

template <typename T>
class removeProjectionPolePositions
{
        private:
                const T complex_crit;

        public:
                removeProjectionPolePositions ( const T complex_crit_ ) : complex_crit ( complex_crit_ ) { }

                bool operator () ( const TProjectionPolePosition <T> & p ) const
                {
                        return ! ( p.complex_crit < complex_crit );
                }

};


#endif
