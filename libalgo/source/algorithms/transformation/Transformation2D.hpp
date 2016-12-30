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


#ifndef Transformation2D_HPP
#define Transformation2D_HPP

#include <cmath>
#include <ctime>
#include <algorithm>

#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"
#include "libalgo/source/algorithms/pointellipseposition/PointEllipsePosition.h"
#include "libalgo/source/algorithms/cartdistortion/CartDistortion.h"

#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByY.h"
#include "libalgo/source/comparators/sortPointPairsByResiduals.h"
#include "libalgo/source/comparators/sortPointPairsByIndices.h"
#include "libalgo/source/comparators/getResidualXY.h"



template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                TTransformationKeyHelmert2D <typename Point1::Type> & key, const bool print_exception, std::ostream * output )
{
        //Compute Helmert transformation, overloaded function for Helmert key
        HelmertTransformation2D::transformPoints ( global_points, local_points, transformed_points, key );
}


template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                TTransformationKeyHomothetic2D <typename Point1::Type> & key, const bool print_exception, std::ostream * output )
{
        //Compute Homothetic transformation, overloaded function for Homothetic key
        HomotheticTransformation2D::transformPoints ( global_points, local_points, transformed_points, key );
}


template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                typename TWeights <typename Point1::Type> ::Type & weights, TTransformationKeyHelmert2D <typename Point1::Type> & key, const bool print_exception, std::ostream * output )
{
        //Compute weighted Helmert transformation, overloaded function for Helmert key
        HelmertTransformation2D::transformPoints ( global_points, local_points, transformed_points, weights, key );
}


template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transformPoints ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points,
                typename TWeights <typename Point1::Type> ::Type & weights, TTransformationKeyHomothetic2D <typename Point1::Type> & key, const bool print_exception, std::ostream * output )
{
        //Compute weighted Homothetic transformation, overloaded function for Homothetic key
        HomotheticTransformation2D::transformPoints ( global_points, local_points, transformed_points, weights, key );
}


template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
void Transformation2D::getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, TTransformationKeyHelmert2D <typename Point1::Type> &key )
{
        //Compute Helmert transformation key, overloaded function for Helmert key
        HelmertTransformation2D::getTransformKey ( global_points, local_points, key );
}


template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
void Transformation2D::getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, TTransformationKeyHomothetic2D <typename Point1::Type> &key )
{
        //Compute Homothetic transformation key, overloaded overloaded function for Homothetic key
        HomotheticTransformation2D::getTransformKey ( global_points, local_points, key );
}


template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
void Transformation2D::getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, typename TWeights <typename Point1::Type> ::Type & weights,
                TTransformationKeyHelmert2D <typename Point1::Type> &key )
{
        //Compute weighted Helmert transformation key, overloaded function for Helmert key
        HelmertTransformation2D::getTransformKey ( global_points, local_points, weights, key );
}


template <typename Point1, typename Point2, TDestructable destructable, TDestructable destructable2>
void Transformation2D::getTransformKey ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, typename TWeights <typename Point1::Type> ::Type & weights,
                TTransformationKeyHomothetic2D <typename Point1::Type> &key )
{
        //Compute weoghted Homothetic transformation key, overloaded overloaded function for Homothetic key
        HomotheticTransformation2D::getTransformKey ( global_points, local_points, weights, key );
}


template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transform ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points, const TTransformationKeyHelmert2D <typename Point1::Type> & key )
{
        //Transform points using Helmert transformation, overloaded overloaded function for Helmert key
        HelmertTransformation2D::transform ( global_points, local_points, transformed_points, key );
}


template <typename Point1, typename Point2, typename Point3, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
void Transformation2D::transform ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points, Container <Point3 *, destructable3> &transformed_points, const TTransformationKeyHomothetic2D <typename Point1::Type> & key )
{
        //Transform points using Homothetic transformation, overloaded overloaded function for Homothetic key
        HomotheticTransformation2D::transform ( global_points, local_points, transformed_points, key );
}

template <typename Point1, typename Point2>
void Transformation2D::rearrangePoints ( const Container <Point1 *> &global_source, const Container <Point2 *> &local_source, Container <Point1 *> &global_destination,
                Container <Point2 *> &local_destination, const typename TDevIndexPairs <typename Point1::Type> ::Type & pairs )
{
        //Rearrange points from destination container according to index of best pairs
        const unsigned int k_pairs_best = pairs.size();

        //Clear containers
        global_destination.clear();
        local_destination.clear();

        //Process all k-best pairs: add global points->
        for ( unsigned int j = 0; j < k_pairs_best; j++ )
        {
                //Add k-best point to the list
                global_destination.push_back ( global_source [ pairs[j].second ] ->clone() );
        }

        //Process all k-best pairs: add local points
        for ( unsigned int j = 0; j < k_pairs_best; j++ )
        {
                //Add k-best point to the list
                local_destination.push_back ( local_source [ pairs[j].second ] ->clone() );
        }
}


template <typename Point1, typename Point2, typename Point3, typename TKey, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
TAccuracyCharacteristics <typename Point1::Type> Transformation2D::getAccuracyCharacteristics ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points,
                const Container <Point3 *, destructable3> &transformed_points, TKey & key )
{
        //Compute accuracy characteristic for non-weighted transformation
        typename TWeights <typename Point1::Type> ::Type weights ( global_points.size(), 1 );
        return getAccuracyCharacteristics ( global_points, local_points, transformed_points, key, weights );
}


template <typename Point1, typename Point2, typename Point3, typename TKey, TDestructable destructable, TDestructable destructable2, TDestructable destructable3>
TAccuracyCharacteristics <typename Point1::Type> Transformation2D::getAccuracyCharacteristics ( const Container <Point1 *, destructable> &global_points, const Container <Point2 *, destructable2> &local_points,
                const Container <Point3 *, destructable3> &transformed_points, TKey & key, typename TWeights <typename Point1::Type> ::Type & weights )
{
        //Compute accuracy characteristic for weighted transformation
        typename Point1::Type residuals_root = 0, weights_sum = 0;
        TAccuracyCharacteristics <typename Point1::Type> accuracy_char;

        //Total points
        const unsigned int n_global = global_points.size(),
                           n_transformed = transformed_points.size();

        //Different size of both lists
        if ( n_global != n_transformed )
        {
                throw BadDataException ( "BadDataException: can not compute standard deviation ( Helmert transformation ), ", "different size of local and global lists of points." );
        }

        //Compute residuals for each point
        for ( unsigned int i = 0; i < n_global; i++ )
        {
                TResidual <typename Point1::Type> residuals;

                //Compute residuals res_x, res_y
                residuals.res_x = global_points [i] -> getX() - transformed_points [i] -> getX();
                residuals.res_y = global_points [i] -> getY() - transformed_points [i] -> getY();
                residuals.res_xy = sqrt ( residuals.res_x * residuals.res_x + residuals.res_y * residuals.res_y );

                //Add to the list
                accuracy_char.res.push_back ( residuals );

                //Compute root of residuals
                residuals_root += weights[i] * residuals.res_xy * residuals.res_xy;

                //Sum of weights
                weights_sum += weights[i];
        }

        //Compute standard deviation
        accuracy_char.std_dev =  sqrt ( residuals_root / ( 2 * n_global - 4 ) );

        //Compute accuracy characteristic for each point
        for ( unsigned int i = 0; i < n_global; i++ )
        {
                //Compute accuracy and q_ee
                typename Point1::Type x_red_local = local_points [i]->getX() - key.x_mass_local;
                typename Point1::Type y_red_local = local_points [i]->getY() - key.y_mass_local;
                typename Point1::Type q_xx = ( x_red_local * x_red_local + y_red_local * y_red_local ) / key.J + 1.0 / key.k ;

                //Add to the list
                accuracy_char.q_xx.push_back ( q_xx );
        }

        return accuracy_char;
}


#endif
