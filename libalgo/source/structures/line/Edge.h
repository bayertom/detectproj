// Description: Edge having 2 end points

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

#ifndef Edge_H
#define Edge_H

#include <stdio.h>

//Edge having 2 start points
template <typename Point>
class Edge
{
        protected:
                Point * p1;		//Start point of the edge
                Point * p2;		//End point of the edge

        public:
                Edge() : p1 ( NULL ), p2 ( NULL ) {}
                Edge ( const Point * p1_, const Point * p2_ ) : p1 ( p1_ ), p2 ( p2_ ) {}
                ~Edge() {}

        public:
                Point * getP1() const {return p1;}
                Point * getP2() const {return p2;}
                void setP1 ( const Point * p1_ ) {p1 = p1_;}
                void setP2 ( const Point * p2_ ) {p2 = p2_;}
};

#endif
