// Description: Some hash algorithms

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

#ifndef Hash_HPP
#define Hash_HPP


template <typename T>
unsigned int Hash::hashXY ( const T x, const T y, const unsigned int x_ratio, const unsigned int y_ratio )
{
        //Hash coordinates x,y (for cm accuracy ratio = 100
        const unsigned int x_hash = ( ( unsigned int ) ( x_ratio * x ) ) & 0xFFFF;
        const unsigned int y_hash = ( ( unsigned int ) ( y_ratio * y ) ) & 0xFFFF;

        return ( x_hash << 16 ) | y_hash;
}


template <typename T>
unsigned int Hash::hashXYZ ( const T x, const T  y, const T z, const unsigned int x_ratio, const T y_ratio, const T z_ratio )
{
        //Hash coordinates x,y,z (for cm accuracy ratio = 100)
        const unsigned int x_hash = ( ( 0xcb1ab31f * ( int ( x_ratio * x ) ) ) >> 16 ) ^ ( int ( x_ratio * x ) );
        const unsigned int y_hash = ( ( 0xd8163841 * ( int ( y_ratio * y ) ) ) >> 16 ) ^ ( int ( y_ratio * y ) );
        const unsigned int z_hash = ( ( 0x6b4a8be7 * ( int ( z_ratio * z ) ) ) >> 16 ) ^ ( int ( z_ratio * z ) );

        return ( x_hash + y_hash + z_hash ) % 16;
}


template <typename T>
int Hash::hashLatLon ( const T lat, const T  lon )
{
        //Hash geographic coordinates lat, lon
        return ( ( int ) ( lat * 1000 ) * 10000 ) + ( int ) ( lon * 1000 );
}

#endif
