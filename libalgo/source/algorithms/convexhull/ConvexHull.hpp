// Description: 2D convex hull using Q-hull algorithm

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


#ifndef ConvexHull_HPP
#define ConvexHull_HPP

#include <algorithm>

#include "libalgo/source/algorithms/pointlineposition/PointLinePosition.h"
#include "libalgo/source/algorithms/pointlinedistance/PointLineDistance.h"

#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByY.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename Point>
void ConvexHull::getConvexHull ( const Container <Point> &nl, Container <Point, NonDestructable > &hull, const bool print_exception, std::ostream * output )
{
        //Construct Convex Hull using Q-HULL algorithm
        try
        {
                //More than 2 points, it could be possible to create convex hull
                if ( nl.size() > 2 )
                {
                        //Find MBR points
                        const unsigned int pos_p_x_max = std::max_element ( nl.begin(), nl.end(), sortPointsByX () ) - nl.begin();
                        const unsigned int pos_p_x_min = std::min_element ( nl.begin(), nl.end(), sortPointsByX () ) - nl.begin();
                        const unsigned int pos_p_y_max = std::max_element ( nl.begin(), nl.end(), sortPointsByY () ) - nl.begin();
                        const unsigned int pos_p_y_min = std::min_element ( nl.begin(), nl.end(), sortPointsByY () ) - nl.begin();

                        //Create four indices lists
                        TIndexList il1, il2, il3, il4;
                        unsigned int n1 = 0, n2 = 0, n3 = 0, n4 = 0;

                        //Split nl to 4 nl inices lists an add points
                        for ( unsigned int i = 0; i < nl.size(); i++ )
                        {
                                //Point is right from existing edge (x_min, y_min)
                                if ( ( PointLinePosition::getPointLinePosition2D ( nl [i], nl [pos_p_x_min], nl [pos_p_y_min], 0 ) == 0 ) )
                                {
                                        //Add to the nl1
                                        il1.push_back ( i );
                                        n1++;
                                }

                                //Point is right from existing edge  (y_min, x_max)
                                if ( ( PointLinePosition::getPointLinePosition2D ( nl [i], nl [pos_p_y_min], nl [pos_p_x_max], 0 ) == 0 ) )
                                {
                                        //Add to the nl2
                                        il2.push_back ( i );
                                        n2++;
                                }

                                //Point is right from existing edge  (x_max, y_max)
                                if ( ( PointLinePosition::getPointLinePosition2D ( nl [i], nl [pos_p_x_max], nl [pos_p_y_max], 0 ) == 0 ) )
                                {
                                        //Add to the nl3
                                        il3.push_back ( i );
                                        n3++;
                                }

                                //Point is right from existing edge (y_max, x_mmin)
                                if ( ( PointLinePosition::getPointLinePosition2D ( nl [i], nl [pos_p_y_max], nl [pos_p_x_min], 0 ) == 0 ) )
                                {
                                        //Add to the nl4
                                        il4.push_back ( i );
                                        n4++;
                                }
                        }


                        //Throw error, only two extremal points and no other point right from their straight line (i.e. all points lies on the line) has been found
                        if ( ( n1 + n2 + n3 + n4 < 1 ) && ( ( pos_p_x_min == pos_p_y_min ) && ( pos_p_x_max  == pos_p_y_max )  ||
                                                            ( pos_p_x_min == pos_p_y_max ) && ( pos_p_x_max  == pos_p_y_min ) ) )
                        {
                                throw BadDataException ( "BadDataException, can not create convex hull, ", "all points lies on the straight line." );
                        }

                        //Approximation of the Convex Hull using quadrilateral connecting MBR (call Q-HULL for the first edge)
                        if ( pos_p_x_min != pos_p_y_min )
                        {
                                //Add point to the CH
                                hull.push_back ( nl [pos_p_x_min] ) ;

                                //Call Q-HULL for the first edge
                                if ( n1 > 0 )
                                {
                                        QHULL ( pos_p_x_min, pos_p_y_min, nl, il1, hull, n1 );
                                }
                        }

                        //Approximation of the Convex Hull using quadrilateral connecting MBR (call Q-HULL for the second edge)
                        if ( pos_p_y_min != pos_p_x_max )
                        {
                                //Add point to the CH
                                hull->push_back ( nl [pos_p_y_min] );

                                //Call Q-HULL for the second edge
                                if ( n2 > 0 )
                                {
                                        QHULL ( pos_p_y_min, pos_p_x_max, nl, il2, hull, n2 );
                                }
                        }

                        //Approximation of the Convex Hull using quadrilateral connecting MBR (call Q-HULL for the third edge)
                        if ( pos_p_x_max != pos_p_y_max )
                        {
                                //Add point to the CH
                                hull.push_back ( nl [pos_p_x_max] );

                                //Call Q-HULL for the third edge
                                if ( n3 > 0 )
                                {
                                        QHULL ( pos_p_x_max, pos_p_y_max, nl, il3, hull, n3 );
                                }
                        }

                        //Approximation of the Convex Hull using quadrilateral connecting MBR (call Q-HULL for the fourth edge)
                        if ( pos_p_y_max != pos_p_x_min )
                        {
                                //Add point to the CH
                                hull.push_back ( nl [pos_p_y_max] );

                                //Call Q-HULL for the fourth edge
                                if ( n4 > 0 )
                                {
                                        QHULL ( pos_p_y_max, pos_p_x_min, nl, il4, hull, n4 );
                                }
                        }
                }

                //Throw exception, not enough points
                else
                {
                        throw BadDataException ( "BadDataException: can not create Convex Hull,", "not enough points: points < 3." );
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                //Print exception
                if ( print_exception )
                {
                        error.printException();
                        *output << "Can not create convex hull ... '\n' ";
                }

                //Clear hu;;
                hull->clear();

                //Throw exception
                throw;
        }
}


template <typename Point>
void ConvexHull::QHULL ( const unsigned int p1, const unsigned int p2, const Container <Point> &nl, const TIndexList & il, Container <Point, NonDestructable > &hull, const unsigned int n )
{

        //Construct Convex hull using modified Q-HULL
        typename Point::Type max_dist_point_line = 0;

        //No point between p1 and p2, stop recursion, no further split
        if ( n == 0 )
        {
                return;
        }

        //1 point between p1 and p2, stop recursion, no further split
        else if ( n == 1 )
        {
                //Add point to the hull
                hull.push_back ( nl [ il [0]] );

                return;
        }

        //More than 1 point between p1 and p2, split set
        else
        {
                //Initialize point index right from the edge (p1, p2)
                int p3 = -1;

                //Find the furthest point on the right side of the edge (p1 - p2)
                for ( unsigned int i = 0; i < n; i++ )
                {
                        //Point is on the right side of the edge (p1, p2)
                        if ( ( PointLinePosition::getPointLinePosition2D ( nl [ il [i]], nl [p1], nl [p2], 0 ) == 0 ) )
                        {
                                //Get point - line distance
                                typename Point::Type dist_point_line = PointLineDistance::getPointLineDistance2D ( nl [ il [i] ], nl [p1], nl [p2] );

                                //Find points with min distance
                                if ( dist_point_line > max_dist_point_line )
                                {
                                        max_dist_point_line = dist_point_line;
                                        p3 = il [i];
                                }
                        }
                }

                //No point on the right to the edge (p1, p2) has been found, no further split
                if ( p3 == -1 )
                {
                        return;
                }

                //Add found point p3 to the Convex hull
                hull.push_back ( nl [p3] );

                //Split set into two subsets located on the right sides of edges
                TIndexList il_first, il_second;
                unsigned int n_first = 0, n_second = 0;

                //Find points right from (p1, p3) or (p3, p2)
                for ( unsigned int i = 0; i < n; i++ )
                {
                        //Point is right from the edge (p1, p3)
                        if ( PointLinePosition::getPointLinePosition2D ( nl [ il [i]], nl [p1], nl [p3], 0 ) == 0 )
                        {
                                il_first.push_back ( il [i] );
                                n_first ++;
                        }

                        //Point is right from the edge (p3, p2)
                        if ( PointLinePosition::getPointLinePosition2D ( nl [ il [i]], nl [p3], nl [p2], 0 ) == 0 )
                        {
                                il_second.push_back ( il [i] );
                                n_second ++;
                        }
                }

                //Recursive calling of the Q-HULL for both intervals
                if ( n_first > 0 ) //Process left subset
                {
                        QHULL ( p1, p3, nl, il_first, hull, n_first );
                }

                if ( n_second > 0 ) //Process right subset
                {

                        QHULL ( p3, p2, nl, il_second, hull, n_second );
                }
        }
}

#endif
