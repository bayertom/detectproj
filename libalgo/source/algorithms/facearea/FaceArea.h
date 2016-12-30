// Description: Compute area of the face

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

#ifndef  FaceArea_H
#define  FaceArea_H

#include "libalgo/source/structures/list/Container.h"

//Forward declaration
template <typename T>
class Point3DCartesian;

template <typename T>
class Face;


//Compute area of the Face
class FaceArea
{
        public:
                template <typename T>
                static T getFaceArea ( const Face <T> *face );

                template <typename T>
                static T getFaceAreaSigned ( const Face <T> *face );

                template <typename Point, TDestructable destructable>
                static typename Point::Type getFaceArea ( const Container <Point *, destructable> *points );

                template <typename Point, TDestructable destructable>
                static typename Point::Type getFaceAreaSigned ( const Container <Point *, destructable> *points );

};

#include "FaceArea.hpp"

#endif

