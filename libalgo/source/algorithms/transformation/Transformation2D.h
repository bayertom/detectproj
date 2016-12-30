// Description: Several equations related to 2D transformation

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


#ifndef Transformation2D_H
#define Transformation2D_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/transformation/HelmertTransformation2D.h"
#include "libalgo/source/algorithms/transformation/HomotheticTransformation2D.h"


//Forward declaration
class TIndexList;


//Residuals
template <typename T>
struct TResidual
{
        T res_x, res_y, res_xy;
};


//List of res for each global point
template <typename T>
struct TResiduals
{
        typedef std::vector < TResidual <T> > Type;
};


//Accuracy characteristics
template <typename T>
struct TAccuracyCharacteristics
{
        T std_dev;					//Standard deviation
        typename TResiduals <T>::Type res;		//Residuals for each point
        typename TItemsList <T>::Type q_xx;		//Diagonal matrix A*inv(A'*W*A)*A' for each point (diagonal)

        TAccuracyCharacteristics () : std_dev ( 0 ), res ( 0 ), q_xx ( 0 ) {}
};


//Pair <deviation, point_index > of an identical point
template <typename T>
struct TDevIndexPair
{
        typedef std::pair <T, unsigned int> Type;
};


//List of pairs of deviatioans and identical points
template <typename T>
struct TDevIndexPairs
{
        typedef std::vector <typename TDevIndexPair <T> ::Type> Type;
};


//Several equations related to 2D transformation
class Transformation2D
{
        public:

		template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              TTransformationKeyHelmert2D <typename Point1::Type> & key, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              TTransformationKeyHomothetic2D <typename Point1::Type> & key, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              typename TWeights <typename Point1::Type> ::Type & weights, TTransformationKeyHelmert2D <typename Point1::Type> & key, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                                              typename TWeights <typename Point1::Type> ::Type & weights, TTransformationKeyHomothetic2D <typename Point1::Type> & key, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transform ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points, const TTransformationKeyHelmert2D <typename Point1::Type> & key );

                template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static void transform ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points, const TTransformationKeyHomothetic2D <typename Point1::Type> & key );

               
                template <typename Point1, typename Point2, typename Point3, typename TKey, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static TAccuracyCharacteristics <typename Point1::Type> getAccuracyCharacteristics ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points,
                                const Container <Point3 *, destructable3> &transformed_points, TKey & key );

                template <typename Point1, typename Point2, typename Point3, typename TKey, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
                static TAccuracyCharacteristics <typename Point1::Type> getAccuracyCharacteristics ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points,
                                const Container <Point3 *, destructable3> &transformed_points, TKey & key, typename TWeights <typename Point1::Type> ::Type & weights );

		
                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, TTransformationKeyHelmert2D <typename Point1::Type> &key );

                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, TTransformationKeyHomothetic2D <typename Point1::Type> &key );

                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, typename TWeights <typename Point1::Type> ::Type & weights,
                                              TTransformationKeyHelmert2D <typename Point1::Type> &key );

                template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
                static void getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, typename TWeights <typename Point1::Type> ::Type & weights,
                                              TTransformationKeyHomothetic2D <typename Point1::Type> &key );


                template <typename Point1, typename Point2>
                static void rearrangePoints ( const Container <Point1 *> &global_source, const Container <Point2 *> &local_source, Container <Point1 *> &global_destination,
                                              Container <Point2 *> &local_destination, const typename TDevIndexPairs <typename Point1::Type> ::Type & pairs );

                
};

#include "Transformation2D.hpp"

#endif
