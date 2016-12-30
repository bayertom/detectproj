// Description: Some has algorithms

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


#ifndef Hash_H
#define Hash_H

//Hashing algorithms
class Hash
{
        public:
                template <typename T>
                static  unsigned int hashXY ( const T x, const T y, const unsigned int x_ratio, const unsigned int y_ratio );

                template <typename T>
                static  unsigned int hashXYZ ( const T x, const T  y, const T z, const unsigned int x_ratio, const T y_ratio, const T z_ratio );

                template <typename T>
                static int hashLatLon ( const T lat, const T lon );

};


#include "Hash.hpp"

#endif
