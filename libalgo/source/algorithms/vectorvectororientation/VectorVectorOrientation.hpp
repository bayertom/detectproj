// Description: Compute 2 vector orientations

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


#ifndef VectorVectorOrientation_HPP
#define VectorVectorOrientation_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/algorithms/pointlinedistance/PointLineDistance.h"

template <typename T>
short VectorVectorOrientation::getVectorVectorOrientation2D ( const T dx1, const T dy1, const T dx2, const T dy2 )
{
        //Return orientation of two vectors
        //		1 = Vectors are positively oriented,
        //		0 = Vectors are negatively oriented,
        // 		2 = Vectors are colinear.

        //Compute signed distance (point, line)
        const T dist_signed = PointLineDistance::getPointLineDistance2DSigned ( dx2, dy2, 0.0, 0.0, dx1, dy1 );

        //Point is too close from edge or lies on the edgeresult 
        if ( fabs ( dist_signed ) < POSITION_ROUND_ERROR )
        {
                return 2;
        }

        //Vectors are positively oriented
        if ( dist_signed > 0.0f )
        {
                return 1;
        }

        //Vectors are negatively oriented
        return 0;
}

#endif
