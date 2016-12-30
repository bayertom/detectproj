// Description: Sort points by distances from the center of gravity

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


#ifndef sortPointsByGravityCenterDistance_H
#define sortPointsByGravityCenterDistance_H


#include "libalgo/source/structures/list/Container.h"
#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"


//Sort points according to distances from the center of gravity (Generic sorter)
template <typename Point>
class sortPointsByGravityCenterDistance
{
        private:
                double xt, yt;
                const Container <Point> *l;

        public:
                sortPointsByGravityCenterDistance ( const Container <Point> *l_ ) : xt ( 0 ), yt ( 0 ), l ( l_ )
                {
                        unsigned int n = this->items.size();

                        //Compute coordinates of the center of the gravity
                        for ( unsigned int i = 0; i < n; i++ )
                        {
                                xt += ( *l ) [i].getX();
                                yt += ( *l ) [i].getY();
                        }

                        xt = 1.0 / n * xt;
                        yt = 1.0 / n * yt;
                }


                bool operator() ( const Point & p1, const Point & p2 ) const
                {
                        //Sort points according to distances from the center of gravity
                        return 	EuclDistance::getEuclDistance2D ( p1.getX(), xt, p1.getY(), yt ) <
                                EuclDistance::getEuclDistance2D ( p2.getX(), xt, p2.getY(), yt );
                }
};


//Sort points according to distances from the center of gravity (Partial specialization for Point*)
template <typename Point>
class sortPointsByGravityCenterDistance <Point*>
{
        private:
                double xt, yt;
                const Container <Point> *l;

        public:
                sortPointsByGravityCenterDistance ( const Container <Point> *l_ ) : xt ( 0 ), yt ( 0 ), l ( l_ )
                {
                        unsigned int n = this->items.size();

                        //Compute coordinates of the center of the gravity
                        for ( unsigned int i = 0; i < n; i++ )
                        {
                                xt += ( *l ) [i]->getX();
                                yt += ( *l ) [i]->getY();
                        }

                        xt = 1.0 / n * xt;
                        yt = 1.0 / n * yt;
                }

                bool operator() ( const Point * p1, const Point * p2 ) const
                {
                        //Sort points according to distances from the center of gravity
                        return 	EuclDistance::getEuclDistance2D ( p1->getX(), xt, p1->getY(), yt ) <
                                EuclDistance::getEuclDistance2D ( p2->getX(), xt, p2->getY(), yt );
                }
};

#endif
