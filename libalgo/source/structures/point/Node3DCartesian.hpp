// Description: 3D cartesian node storing pointers to face / half edge used for topological operations

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


#ifndef Node3DCartesian_HPP
#define Node3DCartesian_HPP


#include <iostream>
#include <iomanip>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/line/HalfEdge.h"

#include "libalgo/source/structures/face/Face.h"



template <typename T>
bool Node3DCartesian <T> ::operator == ( const Node3DCartesian <T> &p ) const
{
        return ( static_cast <Point3DCartesian <T> > ( *this ) ) == ( static_cast <Point3DCartesian <T> > ( p ) );
}


template <typename T>
void Node3DCartesian <T>::print ( std::ostream * file ) const
{
        //Print point into file
        Point3DCartesian <T>::print ( file );
}

#endif
