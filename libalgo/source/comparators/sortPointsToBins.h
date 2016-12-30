// Description: Sort points to bins (Sloan, 1992)

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

#ifndef sortPointsToBins_H
#define sortPointsToBins_H

#include <cmath>

#include <algorithm>

//Sort points to bins: heuristic for incremental DT 2D algorithm (Generic sorter)
template <typename Point>
class sortPointsToBins
{
        private:
                unsigned int n;
                double xmin, ymin, xmax, ymax;

        public:

                sortPointsToBins ( const double xmin_, const double ymin_, const double xmax_, const double ymax_, const int n_ ) :
                        xmin ( xmin_ ), ymin ( ymin_ ), xmax ( xmax_ ), ymax ( ymax_ ), n ( n_ ) {}

                bool operator() ( const Point & p1, const Point & p2 ) const
                {
                        //Sort points to bins (Sloan, 1993)
                        double dmax = ( std::max ) ( xmax - xmin, ymax - ymin );
                        double  x_max = ( xmax - xmin ) / dmax;
                        double  y_max = ( ymax - ymin ) / dmax;

                        //Reduce coordinates
                        double  x_1 = ( p1.getX() - xmin ) / dmax;
                        double  y_1 = ( p1.getY() - ymin ) / dmax;
                        double  x_2 = ( p2.getX() - xmin ) / dmax;
                        double  y_2 = ( p2.getY() - ymin ) / dmax;

                        //Indices
                        unsigned int i1 = 0, j1 = 0, i2 = 0, j2 = 0;

                        if ( x_max != 0 && y_max != 0 )
                        {
                                //First point index
                                i1 = ( unsigned int ) ( 0.99 * n * y_1 / y_max );
                                j1 = ( unsigned int ) ( 0.99 * n * x_1 / x_max );

                                //Second point index
                                i2 = ( unsigned int ) ( 0.99 * n * y_2 / y_max );
                                j2 = ( unsigned int ) ( 0.99 * n * x_2 / x_max );
                        }

                        //Is index "i1" odd or even?
                        double b1, b2;

                        if ( i1 % 2 == 0 )
                        {
                                b1 = i1 * n + j1 + 1;
                        }
                        else
                        {
                                b1 = ( i1 + 1 ) * n - j1;
                        }

                        //Is index "i2" odd or even?
                        if ( i2 % 2 == 0 )
                        {
                                b2 = i2 * n + j2 + 1;
                        }
                        else
                        {
                                b2 = ( i2 + 1 ) * n - j2;
                        }

                        return b1 < b2;
                }
};



//Sort points to bins: heuristic for incremental DT 2D algorithm (Partial specialization for Point *)
template <typename Point>
class sortPointsToBins <Point *>
{
        private:
                unsigned int n;
                double xmin, ymin, xmax, ymax;

        public:

                sortPointsToBins ( const double xmin_, const double ymin_, const double xmax_, const double ymax_, const unsigned int n_ ) :
                        xmin ( xmin_ ), ymin ( ymin_ ), xmax ( xmax_ ), ymax ( ymax_ ), n ( n_ ) {}

                bool operator() ( const Point * p1, const Point * p2 ) const
                {
                        //Sort points to bins (Sloan, 1993)
                        double dmax = ( std::max ) ( xmax - xmin, ymax - ymin );
                        double  x_max = ( xmax - xmin ) / dmax;
                        double  y_max = ( ymax - ymin ) / dmax;

                        //Reduce coordinates
                        double  x_1 = ( p1->getX() - xmin ) / dmax;
                        double  y_1 = ( p1->getY() - ymin ) / dmax;
                        double  x_2 = ( p2->getX() - xmin ) / dmax;
                        double  y_2 = ( p2->getY() - ymin ) / dmax;

                        //Indices
                        unsigned int i1 = 0, j1 = 0, i2 = 0, j2 = 0;

                        if ( x_max != 0 && y_max != 0 )
                        {
                                //First point index
                                i1 = ( unsigned int ) ( 0.99 * n * y_1 / y_max );
                                j1 = ( unsigned int ) ( 0.99 * n * x_1 / x_max );

                                //Second point index
                                i2 = ( unsigned int ) ( 0.99 * n * y_2 / y_max );
                                j2 = ( unsigned int ) ( 0.99 * n * x_2 / x_max );
                        }

                        //Is index "i1" odd or even?
                        double b1, b2;

                        if ( i1 % 2 == 0 )
                        {
                                b1 = i1 * n + j1 + 1;
                        }
                        else
                        {
                                b1 = ( i1 + 1 ) * n - j1;
                        }

                        //Is index "i2" odd or even?
                        if ( i2 % 2 == 0 )
                        {
                                b2 = i2 * n + j2 + 1;
                        }
                        else
                        {
                                b2 = ( i2 + 1 ) * n - j2;
                        }

                        return b1 < b2;
                }
};




#endif

