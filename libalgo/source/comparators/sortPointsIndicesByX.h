// Description: Sort indices of cartesian points stored in list by X coordinate

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


#ifndef sortPointsIndicesByX_H
#define sortPointsIndicesByX_H

#include "libalgo/source/structures/list/Container.h"

//Sorter of points defined using indices by x coordinate and subsequently by y coordinate (Generic sorter)
template <typename Point>
class sortPointsIndicesByX
{
        private:
                Container <Point> pl;

        public:
                sortPointsIndicesByX ( const Container <Point> &pl_ ) : pl ( pl_ ) {}

                bool operator() ( const unsigned int & i_p1, const unsigned int & i_p2 ) const
                {
                        return ( pl [i_p1].getX() < pl [i_p2].getX() ) || ( ( pl [i_p1].getX() == pl [i_p2].getX() ) && ( pl [i_p1].getY() < pl [i_p2].getY() ) );
                }
};


//Sorter of points defined using indices by x coordinate and subsequently by y coordinate (Partial specialization for Point *)
template <typename Point >
class sortPointsIndicesByX <Point *>
{

        private:
                const Container <Point *> pl;

        public:
                sortPointsIndicesByX ( const Container <Point *> &pl_ ) : pl ( pl_ ) {}

                bool operator() ( const unsigned int & i_p1, const unsigned int & i_p2 ) const
                {
                        return ( pl [i_p1]->getX() < pl [i_p2]->getX() ) || ( ( pl [i_p1]->getX() == pl [i_p2]->getX() ) && ( pl [i_p1]->getY() < pl [i_p2]->getY() ) );
                }

};

#endif
