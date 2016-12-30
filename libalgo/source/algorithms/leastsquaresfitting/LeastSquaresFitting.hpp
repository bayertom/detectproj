// Description: Least squares fitting, 2D line

// Copyright (c) 2010 - 2011
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


#ifndef LeastSquaresFitting_HPP
#define LeastSquaresFitting_HPP

#include <cmath>


#include "libalgo/source/exceptions/BadDataException.h"
#include "libalgo/source/exceptions/IndexOutOfBoundException.h"


template <typename Point>
void LeastSquaresFitting::fitLine ( const Container <Point> &nl, TFittingLine <Point> &regression_line, const bool print_message, const bool print_exception, std::ostream * output )
{
        //Find best fitting line using least squares approximation
        const unsigned int regression_line_points = regression_line.points_indices.size();
        const unsigned int total_points = nl.size();

        try
        {
                //Are there enough points?
                if ( ( regression_line_points >= 2 ) && ( total_points >= 2 ) )
                {
                        typename Point::Type sumx = 0.0, sumy = 0.0;
                        typename TIndexSet <Point> ::Type ::iterator i_points_indices = regression_line.points_indices.begin();

                        //Calculate coordinates of the center of gravity
                        for ( i_points_indices = regression_line.points_indices.begin(); i_points_indices != regression_line.points_indices.end(); ++i_points_indices )
                        {
                                sumx += nl [*i_points_indices].getX();
                                sumy += nl [*i_points_indices].getY();
                        }

                        //Center of the gravity
                        regression_line.xt = sumx / regression_line_points;
                        regression_line.yt = sumy / regression_line_points;

                        //Calculate sums
                        typename Point::Type sumxx = 0.0, sumxy = 0.0, sumyy = 0.0;

                        for ( i_points_indices = regression_line.points_indices.begin(); i_points_indices != regression_line.points_indices.end(); ++i_points_indices )
                        {
                                //Get coordinaets of the point
                                const typename Point::Type x = nl [*i_points_indices].getX(),
                                                           y = nl [*i_points_indices].getY();

                                //Sums
                                sumxx += ( x - regression_line.xt ) * ( x - regression_line.xt );
                                sumxy += ( x - regression_line.xt ) * ( y - regression_line.yt );
                                sumyy += ( y - regression_line.yt ) * ( y - regression_line.yt );
                        }

                        //Calculate numerator and denominator
                        const typename Point::Type denom = sumxx - sumyy;

                        //Compute alpha2
                        const typename Point::Type alpha2 = atan2 ( 2.0 * sumxy, denom ) * 180 / M_PI;

                        //Compute alpha
                        regression_line.alpha = ( alpha2 < 0 ? 0.5 * alpha2 + 360.0 : 0.5 * alpha2 );

                        //Horizontal line
                        if ( ( regression_line.alpha == 0.0 ) || ( regression_line.alpha == 180.0 ) )
                        {
                                regression_line.b = regression_line.yt;
                                regression_line.c = 1;
                        }

                        //Vertical line
                        else if ( ( regression_line.alpha == 90.0 ) || ( regression_line.alpha == 270.0 ) )
                        {
                                regression_line.b = regression_line.xt;
                                regression_line.c = 0;
                        }

                        //Common line
                        else
                        {
                                regression_line.b = regression_line.yt - regression_line.xt * tan ( regression_line.alpha * M_PI / 180 );
                                regression_line.c = 1;
                        }

                        //Compute regression error: sum of roots of distances
                        typename Point::Type root_dist = 0.0;

                        for ( i_points_indices = regression_line.points_indices.begin(); i_points_indices != regression_line.points_indices.end(); ++i_points_indices )
                        {

                                //Get coordinates of the point
                                const typename Point::Type x = nl [*i_points_indices].getX(),
                                                           y = nl [*i_points_indices].getY();

                                //Distance of the point to the regression line
                                const typename Point::Type dist = ( x - regression_line.xt ) *  sin ( regression_line.alpha * M_PI / 180 ) - ( y - regression_line.yt ) * cos ( regression_line.alpha * M_PI / 180 );

                                //Error: root of the distance
                                root_dist += dist * dist;
                        }

                        //Compute error
                        regression_line.error = sqrt ( root_dist / regression_line_points );
                }

                //Throw exception
                else
                {
                        throw BadDataException ( "ErorrBadData: ", "Least square fitting, total points n < 3." );
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                //Print exception
                if ( print_exception )
                {
                        error.printException ( output );
                }

                //Throw exception
                throw;
        }
}

#endif
