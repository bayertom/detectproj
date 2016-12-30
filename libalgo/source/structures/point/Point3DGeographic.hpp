// Description: 3D geographic point

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


#ifndef Point3DGeographic_HPP
#define Point3DGeographic_HPP

#include <string.h>

//Include correct version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/const/Const.h"
#else
	#include "libalgo/source/const2/Const.h"
#endif


//Set point ID  to 0, initialize static variable
template <typename T>
unsigned int Point3DGeographic <T>::points_geo_id_counter = 0;

template <typename T>
Point3DGeographic <T> & Point3DGeographic <T> ::operator = ( const Point3DGeographic <T> &p )
{
        if ( this != &p )
        {
                point_id = p.point_id;
		point_label = p.point_label;
                lat = p.lat;
                lon = p.lon;
                H = p.H;
        }

        return *this;
}


template <typename T>
bool Point3DGeographic <T> ::operator == ( const Point3DGeographic <T> &p ) const
{
        return ( lat - p.lat ) * ( lat - p.lat ) + ( lon - p.lon ) * ( lon - p.lon ) < MIN_POSITION_DIFF * MIN_POSITION_DIFF;
}


template <typename T>
void Point3DGeographic <T>::print ( std::ostream * output ) const
{
        //Print point
        *output << std::fixed << std::setprecision ( 3 );
        *output << point_id << "   "  << lat << "   " << lon  << "   " << H << '\n';
}

#endif
