// Description: Cartographic parallel given by start / end pointers

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


#ifndef Parallel_HPP
#define Parallel_HPP

#include <string.h>
#include <algorithm>
#include <iostream>

//#include "point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Node3DCartesianProjected.h"


template <typename T>
Parallel <T>::Parallel ( const Parallel <T> &p )
{
        this->lat = p.lat;
        this->points_indices = p.points_indices;
}


template <typename T>
template <typename Point>
Parallel <T> ::Parallel ( const TFittingLine <Point> &line )
{
        //Create new parallel from RANSAC detected line
        lat = line.yt;

        //Copy indices from set to vector
        std::copy ( line.points_indices.begin(), line.points_indices.end(), std::back_inserter ( points_indices ) );
}



template <typename T>
void Parallel <T> ::print ( std::ostream * output ) const
{
        //Print all points of the parallel
        *output << "PARALLEL" << std::endl ;
        *output << "lat = " << lat << std::endl;
        *output << "Points: " << std::endl << std::endl;

        //Print all points
        for ( unsigned int i = 0; i < points_indices.size(); i++ )
        {
                //Print point inedx
                *output << points_indices[i] << '\t';
        }

        *output << std::endl;
}


template <typename T>
template <typename Point>
void Parallel <T> ::print ( const Container <Point> *points, std::ostream * output ) const
{
        //Print all points of the parallel
        *output << "PARALLEL" << std::endl ;
        *output << "lat = " << lat << std::endl;
        *output << "Points: " << std::endl << std::endl;

        //Print all points
        for ( unsigned int i = 0; i < points_indices.size(); i++ )
        {
                //Print point
                *output << ( *points ) [points_indices[i]]->getLat() << '\n';
                *output << ( *points ) [points_indices[i]]->getLon() << "\n\n";
        }
}


#endif
