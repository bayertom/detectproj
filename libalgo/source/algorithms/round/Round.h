// Description: Number rounding

// Copyright (c) 2010 - 2016
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


#ifndef Round_H
#define Round_H

#include <cmath>


//Round number
class Round
{
        public:
                template <typename T>
                static T roundNumber ( const T num, const unsigned short dec_places );

		template <typename T>
		static T roundToMultipleCeil( const T value, const T multiple );

		template <typename T>
		static T roundToMultipleFloor( const T value, const T multiple );

};

#include "Round.hpp"

#endif
