// Description: Outliers detection using the least squares and M-estimators

// Copyright (c) 2010 - 2014
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


#ifndef Outliers_H
#define Outliers_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/transformation/Transformation2D.h"

//Scheme for M-estimators
typedef enum
{
	ShiftsScheme = 0, 
	ScaleScheme,
	ScaleShiftsScheme, 
	SimilarityScheme,
} TMEstimatorsScheme;


//Weight function for M-estimators
typedef enum
{
	HuberFunction = 0,
	AndrewFunction,
	TukeyFunction,
	YangFunction,
	DanishFunction,
	DanishFunction2,
} TMEstimatorsWeightFunction;


//Several equations related to 2D transformation
class Outliers
{
	public:
		template <typename Point1, typename Point2, typename TKey>
                static void findOutliersLTS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest,
                                Container <Point2 *> &local_points_dest, TKey & min_key,  typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const typename Point1::Type perc_ratio = 0.8 );

		template <typename T>
		static void findOutliersME(const Matrix <T> &P, const Matrix <T> &Q, const T k, const T tol, const TMEstimatorsScheme me_scheme, TMEstimatorsWeightFunction me_weight_function, const unsigned int max_iter, Matrix <T> &W, Matrix <unsigned int> &I, Matrix <T> &Eps, T &f_init, T &f, unsigned int &iter);

                template <typename Point1, typename Point2, typename TKey>
                static void findOutliersIRLS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
                                Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest,  TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs );

		template <typename Point1, typename Point2>
		static void findOutliersME(const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, const typename Point1::Type k, const typename Point1::Type tol,
			const TMEstimatorsScheme me_scheme, TMEstimatorsWeightFunction me_weight_function, const unsigned int max_iter, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, typename Point1::Type &f_init, typename Point1::Type &f, unsigned int &iter);


	private:
		template <typename T>
                static void createKBestPairsOfPoints ( const TAccuracyCharacteristics <T> &deviations, typename TDevIndexPairs<T>::Type & point_pairs, const float perc_ratio );



};

#include "Outliers.hpp"

#endif