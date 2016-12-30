// Description: Compute inner distance of the face

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


#ifndef InnerDistance_H
#define InnerDistance_H

#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class Face;

template <typename T>
class Matrix;


//Create matrix of  inner distances and performs analysis
class InnerDistance
{
        public:
                template <typename T>
                static T compare2FacesUsingInnerDistances ( const Face <T> *f1, const Face <T> *f2 );
		
        private:
                template <typename T>
                static void computeTheta ( const Face <T> *f, const Matrix <T> &D, const Matrix <int> &P, Matrix <unsigned short> &TH );
};


#include "InnerDistance.hpp"

#endif
