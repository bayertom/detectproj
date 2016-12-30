// Description: Compute Shape Context used for Inner Distance Analysis

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

#ifndef ShapeContext_H
#define ShapeContext_H


#include "libalgo/source/structures/matrix/Matrix.h"


//Compute shape context
class ShapeContext
{
        public:
                template <typename T>
                static void computeShapeContext ( const Matrix <T> &D, const Matrix <unsigned short> &TT, Matrix <unsigned short> &SC, unsigned int point_index );

                template <typename T>
                static T compare2ShapeContexts ( const Matrix <unsigned short> &SC1, const Matrix <unsigned short> &SC2 );
};

#include "ShapeContext.hpp"

#endif
