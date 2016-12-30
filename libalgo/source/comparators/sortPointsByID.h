// Description: Sort cartesian point by ID

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

#ifndef sortPointsByID_H
#define sortPointsByID_H


//Sorter by ID (Generic sorter)
template <typename Point>
class sortPointsByID
{
        public:

                bool operator() ( const Point & p1, const Point & p2 ) const
                {
                        return ( p1.getPointID() < p2.getPointID() );
                }

};


//Sorter by ID (Partial specialization for Point *)
template <typename Point>
class sortPointsByID <Point *>
{
        public:
                bool operator() ( const Point * p1, const Point * p2 ) const
                {
                        return ( p1->getPointID() < p2->getPointID() );
                }
};

#endif
