// Description: Compute bearing for line given by 3 points

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


#ifndef Bearing_HPP
#define Bearing_HPP

#include <cmath>


//#include "libalgo/source/const/Const.h"

#include "libalgo/source/exceptions/MathInvalidArgumentException.h"


template <typename Point>
typename Point::Type Bearing::getBearing ( const Point * p1, const Point * p2 )
{
        //Compute bearing
        const typename Point::Type dy = p2->getY() - p1->getY();
        const typename Point::Type dx = p2->getX() - p1->getX();

        //Points not identical
        if ( ( dx != 0.0f ) && ( dy != 0.0 ) )
        {
                //Compute phi
                const typename Point::Type phi = atan ( fabs ( dy / dx ) ) * 180 / M_PI ;

                //First quadrant: (0; 90) deg
                if ( ( dx > 0.0 ) && ( dy > 0.0 ) )
                {
                        return phi;
                }

                //Second quadrant: (90; 180) deg
                if ( ( dx < 0.0 ) && ( dy > 0.0 ) )
                {
                        return 180.0 - phi;
                }

                //Third quadrant: (180; 270) deg
                if ( ( dx < 0.0 ) && ( dy < 0.0 ) )
                {
                        return 180.0 + phi;
                }

                //Fourth quadrant: (270; 360) deg
                if ( ( dx > 0.0 ) && ( dy < 0.0 ) )
                {
                        return 360.0 - phi;
                }
        }

        //90 or 180 deg
        else if ( dx == 0.0 )
        {
                //90 deg
                if ( dy > 0.0 )
                {
                        return 90.0;
                }

                //270 deg
                if ( dy < 0.0 )
                {
                        return 270.0;
                }

                //Identical points, throw exception
                if ( dy == 0.0 )
                {
                        throw MathInvalidArgumentException <typename Point::Type> ( "EErorrMathInvalidArgument: ", "can not compute bearing(p1, p2), dist(p1,p2) = ", sqrt ( dx * dx + dy * dy ) );
                }
        }

        // 0 or 90 deg
        else if ( dy == 0.0 )
        {
                //0 deg
                if ( dx > 0.0 )
                {
                        return 0.0;
                }

                //180 deg
                if ( dx < 0.0 )
                {
                        return 180.0;
                }
        }
}

#endif
