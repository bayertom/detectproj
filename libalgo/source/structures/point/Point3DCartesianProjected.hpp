// Description: 3D cartesian cartographic point

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


#ifndef Point3DCartesianProjected_HPP
#define Point3DCartesianProjected_HPP

//Include C++98/C++11 version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/structures/projection/Projection.h"
#else
	#include "libalgo/source/structures/projection2/Projection.h"
#endif


template <typename T>
bool Point3DCartesianProjected <T>:: operator == ( const Point3DCartesianProjected <T> &p ) const
{
        return ( static_cast <const Point3DCartesian <T> &> ( *this ) == static_cast <const Point3DCartesian <T> &> ( p ) );
}


template <typename T>
void Point3DCartesianProjected <T>::print ( std::ostream * output ) const
{
        //Print point
        *output << std::fixed << std::setprecision ( 3 );
        Point3DCartesian <T>::print ( output );
        *output << "h = " << h << " k = " << k << " s = " << s  << " theta = " << theta << " a_tiss = " << tiss.a_tiss << '\n'
                << " b_tiss = " << tiss.b_tiss << " Ae = " << tiss.Ae << " b_mer = " << tiss.b_mer << " w = " << w << '\n';
}

#endif
