// Description: Compute perimeter of the face

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


#ifndef FacePerimeter_H
#define FacePerimeter_H

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class Face;


//Compute perimeter of the Face
class FacePerimeter
{
        public:
                template <typename T>
                static T getFacePerimeter ( const Face <T> *f );
        
                template <typename T>
                static T getFacePerimeter ( const HalfEdge <T> *e_start );
};

#include "FacePerimeter.hpp"

#endif

