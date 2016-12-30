// Description: 2D half edge class

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


#ifndef HalfEdge_HPP
#define HalfEdge_HPP

#include <stdio.h>

#include <iomanip>

#include "libalgo/source/structures/face/Face.h"
#include "libalgo/source/structures/point/Node3DCartesian.h"



template <typename T>
HalfEdge <T> ::~HalfEdge()
{
        //Set pointers
        start = NULL;
        previous = NULL;
        next = NULL;
        twin = NULL;
        dual = NULL;
        face = NULL;
        simplex = false;
}


template <typename T>
void HalfEdge <T> ::print ( std::ostream * output ) const
{
        //Print half edge
        *output << std::fixed << std::setprecision ( 5 );

        *output << "Start: " ;
        start != NULL ?  *output << start-> getX() << "  " << start->getY() << std::endl : *output << "NULL" << std::endl;

        *output << "Previous: ";
        previous != NULL ? *output << previous->getPoint()-> getX() << "  " << previous->getPoint()->getY() << std::endl : *output << "NULL" << std::endl;

        *output << "Next: ";
        next != NULL ? *output << next->getPoint()-> getX() << "  " << next->getPoint()->getY() << std::endl : *output << "NULL" << std::endl;

        *output << "Twin: ";
        twin != NULL ? *output << twin->getPoint()-> getX() << "  " << twin->getPoint()->getY() << std::endl : *output << "NULL" << std::endl;

        *output << "Dual: ";
        dual != NULL ? *output << dual->getPoint()-> getX() << "  " << dual->getPoint()->getY()  << std::endl : *output << "NULL" << std::endl;

        *output << "Simplex: ";
        *output << simplex << std::endl;
}


#endif
