// Description: 2D Helmert weighted / non weighted transformation

// Copyright (c) 2015 - 2016
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

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TTransformationKeyHelmert2D.h"


//2D Helmert (weighted / non weighted ) transformation
class HelmertTransformation2D
{
        public:

                template <typename T, typename Point1, typename Point2, typename Point3>
                static void transformPoints ( const TVector <Point1> &global_points, const TVector <Point2> &local_points, TVector <Point3> &transformed_points, TTransformationKeyHelmert2D <T> &key_helmert);

                template <typename T, typename Point1, typename Point2>
                static void getTransformKey (const TVector <Point1> &global_points, const TVector <Point2> &local_points, const TVector <T> weights, TTransformationKeyHelmert2D <T> &key_helmert );
  
                template <typename T, typename Point1, typename Point2, typename Point3>
                static void transform (const TVector <Point1> &global_points, const TVector <Point2> &local_points, TVector <Point3> &transformed_points, const TTransformationKeyHelmert2D <T> & key_helmert);

		//template <typename T>
		//static Matrix <T> getTransformKey(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W);

		//template <typename T>
		//static void getTransformKey2(const Matrix <T> &P, const Matrix <T> &Q, const Matrix <T> &W, Matrix <T> &A, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &C);
		
};

#include "HelmertTransformation2D.hpp"

#endif

