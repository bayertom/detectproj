// Description: Compute TAR criterion

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

#ifndef TARCriterion_H
#define TARCriterion_H

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class Face;

//Computer Triangulate Area Criterion for the Face
class TARCriterion
{
        public:
                template <typename T>
                static T compare2FacesUsingTARCriterion ( const Face <T> *p1, const Face <T> *p2 );

                template <typename T>
                static T getTARCriterion ( HalfEdge <T> *e, unsigned int vertices_step );

};

#include "TARCriterion.hpp"

#endif

