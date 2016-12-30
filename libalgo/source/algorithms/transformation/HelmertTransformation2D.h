// Description: 2D Helmert weighted / non weighted transformation

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


#ifndef HelmertTransformation2D_H
#define HelmertTransformation2D_H

#include <vector>

#include "libalgo/source/structures/list/Container.h"

//Transformation key
template <typename T>
struct TTransformationKeyHelmert2D
{
        T x_mass_local, y_mass_local;		//Centre of mass: local system
        T x_mass_global, y_mass_global;		//Centre of mass: global system
        T c1, c2;				//Both transformation coefficients
        T J;					// (X^T * X)^-1
        T k;					// (X^T * X)^-1

        TTransformationKeyHelmert2D () : x_mass_local ( 0.0 ), y_mass_local ( 0.0 ), x_mass_global ( 0.0 ),
                y_mass_global ( 0.0 ), c1 ( 1.0 ), c2 ( 1.0 ), J ( 0.0 ), k ( 0.0 ) {}
};

//Weights of points
template <typename T>
struct TWeights
{
        typedef std::vector <T> Type;
};


//2D Helmert (weighted / non weighted ) transformation
class HelmertTransformation2D
{
        public:

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              TTransformationKeyHelmert2D <typename Point1::Type> &key_helmert, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              typename TWeights <typename Point1::Type> ::Type & weights, TTransformationKeyHelmert2D <typename Point1::Type> &key_helmert, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, TTransformationKeyHelmert2D <typename Point1::Type> &key_helmert );

                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, typename TWeights <typename Point1::Type> ::Type & weights,
                                              TTransformationKeyHelmert2D <typename Point1::Type> &key_helmert );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transform ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                        const TTransformationKeyHelmert2D <typename Point1::Type> & key_helmert );

		template <typename T>
		static Matrix <T> getTransformKey(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W);

		template <typename T>
		static void getTransformKey2(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W, Matrix <T> &A, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &C);
		


};

#include "HelmertTransformation2D.hpp"

#endif

