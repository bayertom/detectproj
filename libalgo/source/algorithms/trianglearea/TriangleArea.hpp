// Description: Compute area of the triangle

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


#ifndef TriangleArea_HPP
#define TriangleArea_HPP


template <typename Point>
typename Point::Type TriangleArea::getTriangleArea ( const Point * p1, const Point * p2, const Point * p3 )
{
        //Caculate unsigned triangle area
        return fabs ( getTriangleAreaSigned ( p1, p2, p3 ) );
}


template <typename Point>
typename Point::Type TriangleArea::getTriangleAreaSigned ( const Point * p1, const Point * p2, const Point * p3 )
{
        //Compute signed triangle area
        return 0.5 * ( ( p2->getX() - p1->getX() ) * ( p3->getY() - p1->getY() ) -
                       ( p3->getX() - p1->getX() ) * ( p2->getY() - p1->getY() ) );
}

#endif
