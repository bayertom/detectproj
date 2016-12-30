// Description: Compute 2 vector orientations

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

#ifndef VectorVectorOrientation_H
#define VectorVectorOrientation_H


//Compute orientation of two vectors
class VectorVectorOrientation
{
        public:
                template <typename T>
                static short getVectorVectorOrientation2D ( const T dx1, const T dy1, const T dx2, const T dy2 );
};

#include "VectorVectorOrientation.hpp"

#endif
