// Description: 3D cartesian node storing cartographic parameters, derived class

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


#ifndef Node3DCartesianProjected_HPP
#define Node3DCartesianProjected_HPP

#include <algorithm>

#include "Point3DCartesianProjected.h"
#include "Node3DCartesian.h"


template <typename T>
Node3DCartesianProjected <T> ::Node3DCartesianProjected ( const Node3DCartesianProjected <T> *n ) : Point3DCartesian <T> ( *n ), Node3DCartesian <T> ( *n ), Point3DCartesianProjected <T> ( *n )
{
        parallel_point_index_prev = n->parallel_point_index_prev ;
        parallel_point_index_next = n->parallel_point_index_next;
        meridian_point_index_next = n->meridian_point_index_next;
        meridian_point_index_prev = n->meridian_point_index_prev;
}


template <typename T>
Node3DCartesianProjected <T> ::~Node3DCartesianProjected()
{
        parallel_point_index_prev = -1;
        parallel_point_index_next = -1;
        meridian_point_index_next = -1;
        meridian_point_index_prev = -1;
}


template <typename T>
bool Node3DCartesianProjected <T>::operator == ( const Node3DCartesianProjected <T> &n ) const
{
        return ( static_cast <const Point3DCartesian <T> &> ( *this ) == static_cast <const Point3DCartesian <T> &> ( n ) );
}


template <typename T>
void Node3DCartesianProjected <T> ::print ( std::ostream * output ) const
{
        //Print node
        Node3DCartesian <T>::print ( output );
        /*
        *output << "Parallel point index prev: " << parallel_point_index_prev << std::endl;
        *output << "Parallel point index next: " << parallel_point_index_next << std::endl;
        *output << "Meridian point index prev: " << meridian_point_index_prev << std::endl;
        *output << "Meridian point index next: " << meridian_point_index_next << std::endl;*/
}

#endif
