// Description: Compute tangent function

// Copyright (c) 2010 - 2011
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

#ifndef TangentFunction_H
#define TangentFunction_H

#include <map>


#include "libalgo/source/structures/list/Container.h"


//Forward declarations
template <typename T>
class Face;

template <typename T>
class HalfEdge;

class TIndexList;

//New user type definition
template <typename T>
struct TTangentFunction
{
        typedef std::map <T, T> Type;
};

//Rotation dependent/invariant
typedef enum
{
        RotationInvariant,
        RotationDependent,
} TTangentFunctionRotationMethod;


//Scale dependent/invariant
typedef enum
{
        ScaleInvariant, ScaleDependent
} TTangentFunctionScaleMethod;



//Compute tangent function for the Face
class TangentFunction
{
        public:
                template <typename T>
                static T compare2FacesUsingTangentFunction ( const Face <T> *p1, const Face <T> *p2, const TTangentFunctionRotationMethod & rotation_method = RotationInvariant,
			const TTangentFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename Point>
                static typename Point::Type compare2PolyLinesUsingTangentFunction ( const Container <Point> *pl1, const Container <Point> *pl2, const TTangentFunctionRotationMethod & rotation_method = RotationInvariant,
			const TTangentFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename T>
                static void computeTangentFunctionFace ( const Face <T> *face, typename TTangentFunction <T>::Type & tf, T & rotation, const TTangentFunctionRotationMethod & rotation_method = RotationInvariant,
			const TTangentFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename Point>
                static void computeTangentFunctionPolyLine ( const Container <Point> *points, typename TTangentFunction <typename Point::Type>::Type & tf, const TTangentFunctionRotationMethod & rotation_method = RotationInvariant,
			const TTangentFunctionScaleMethod & scale_method = ScaleInvariant );


        private:

                template <typename T>
                static T compareTangentFunctions ( typename TTangentFunction <T>::Type & tf1,  typename TTangentFunction <T>::Type & tf2 ,  const TTangentFunctionRotationMethod & rotation_method = RotationInvariant );

                template <typename T>
                static void rotateKeys ( typename TTangentFunction <T>::Type & tf );
};

#include "TangentFunction.hpp"

#endif
