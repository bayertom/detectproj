// Description: VArious swapping criteria for data depending triangulations

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

#ifndef SwappingCriteria_H
#define SwappingCriteria_H

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class Point3DCartesian;

//Calculate swapping criteria for DT or DDT
class SwappingCriteria
{
        public:
                template <typename T>
                static T getAbn ( const Point3DCartesian <T> *n1, const Point3DCartesian <T> *n2, const Point3DCartesian <T> *n3, const Point3DCartesian <T> *n4, const Point3DCartesian <T> *n5, const Point3DCartesian <T> *n6 );

                template <typename T>
                static T getSco ( const Point3DCartesian <T> *n1, const Point3DCartesian <T> *n2, const Point3DCartesian <T> *n3, const Point3DCartesian <T> *n4, const Point3DCartesian <T> *n5, const Point3DCartesian <T> *n6 );

                template <typename T>
                static bool getClineRenka ( const Node3DCartesian <T> *p, const Node3DCartesian <T> *p1, const Node3DCartesian <T> *p2, const Node3DCartesian <T> *p3 );

                template <typename T>
                static T getEmptyCircleTest ( const Node3DCartesian <T> *p, const Node3DCartesian <T> *p1, const Node3DCartesian <T> *p2, const Node3DCartesian <T> *p3 );
};

#include "SwappingCriteria.hpp"

#endif
