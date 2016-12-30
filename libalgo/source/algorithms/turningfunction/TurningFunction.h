// Description: Compute turning function

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

#ifndef TurningFunction_H
#define TurningFunction_H

#include <map>
#include <ostream>
#include <iostream>


#include "libalgo/source/structures/list/Container.h"


//Forward declarations
template <typename T>
class Face;

template <typename T>
class HalfEdge;

class TIndexList;

//New user type definition
template <typename T>
struct TTurningFunction
{
        typedef std::map <T, T> Type;
};

//Rotation dependent/invariant
typedef enum
{
        RotationInvariant = 0,
        RotationDependent,
} TTurningFunctionRotationMethod;


//Scale dependent/invariant
typedef enum
{
        ScaleInvariant = 0, 
	ScaleDependent
} TTurningFunctionScaleMethod;



//Compute turning function for the Face
class TurningFunction
{
        public:
                template <typename T>
                static T compare2FacesUsingTurningFunction ( const Face <T> *p1, const Face <T> *p2, const TTurningFunctionRotationMethod & rotation_method = RotationInvariant,
                                const TTurningFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename Point>
                static typename Point::Type compare2PolyLinesUsingTurningFunction ( const Container <Point> &pl1, const Container <Point> &pl2, const TTurningFunctionRotationMethod & rotation_method = RotationInvariant,
                                const TTurningFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename T>
                static void computeTurningFunctionFace ( const Face <T> *face, typename TTurningFunction <T>::Type & tf, T & rotation, const TTurningFunctionRotationMethod & rotation_method = RotationInvariant,
                                const TTurningFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename Point>
                static void computeTurningFunctionPolyLine ( const Container <Point> &points, typename TTurningFunction <typename Point::Type>::Type & tf, const TTurningFunctionRotationMethod & rotation_method = RotationInvariant,
                                const TTurningFunctionScaleMethod & scale_method = ScaleInvariant );

                template <typename T>
                static void print ( std::ostream * output, const typename TTurningFunction <T>::Type & tf );

                template <typename T>
                friend void operator << ( std::ostream & output, const typename TTurningFunction <T>::Type & tf ) { tf.print ( &output, tf );}

        public:

                template <typename T>
                static T compareTurningFunctions ( const typename TTurningFunction <T>::Type & tf1,  const typename TTurningFunction <T>::Type & tf2,  const TTurningFunctionRotationMethod & rotation_method = RotationInvariant );
};

#include "TurningFunction.hpp"

#endif
