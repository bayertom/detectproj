// Description: Basic graph algorithms

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


#ifndef GraphAlgorithms_H
#define GraphAlgorithms_H

#include "libalgo/source/structures/list/Container.h"
#include "libalgo/source/structures/graph/GraphL.h"
#include "libalgo/source/structures/graph/GraphM.h"
#include "libalgo/source/structures/matrix/Matrix.h"

//Forward declaration
template <typename T>
class Node3DCartesian;

template <typename T>
class Face;


//Graph algorithms
class GraphAlgorithms
{
        public:

                template <typename T>
                static void floydWarshall ( const GraphM <T> &g, Matrix <T> &D, Matrix <unsigned int> &P );

                template <typename T>
                static GraphL <T> mst ( GraphL <T> &g, T & weight );

                template <typename T>
                static void bestBipartiteMatching ( Matrix <T> C, Matrix <unsigned short> &M, T & cost );

                template <typename Point>
                static void createKNNGraph ( const Container <Point *> *points, const unsigned int k, GraphM <typename Point::Type> &g );

                template <typename Point>
                static void createNNNGraph ( const Container <Point *> *points, GraphM <typename Point::Type> &g );

                template <typename Point>
                static void createGabrielGraph ( const Container <Point *> *points, GraphM <typename Point::Type> &g );

                template <typename Point>
                static void createSphereOfInfulenceGraph ( const Container <Point *> *points, GraphM <typename Point::Type> &g );

                template <typename T>
                static void FaceToWDistanceGraph ( const Face <T> *f, GraphM <T> &g, const bool normalized = false );

                template <typename T>
                static void FaceToInnerDistanceGraph ( const Face <T> *f, GraphM <T> &g, const bool normalized = false );


        private:

                static unsigned int findSet ( const unsigned int x, std::vector<unsigned int> &parents );

                template <typename T>
                static void bestBipartiteMatchingPhase1 ( Matrix <T> &CR, unsigned short & phase );

                template <typename T>
                static void bestBipartiteMatchingPhase2 ( const Matrix <T> &C, Matrix <unsigned short> &M, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, unsigned short & phase );

                static void bestBipartiteMatchingPhase3 ( const Matrix <unsigned short> &M, Matrix <unsigned short> &C_C, unsigned int min_rows, unsigned short & phase );

                template <typename T>
                static void bestBipartiteMatchingPhase4 ( const Matrix <T> &C, Matrix <unsigned short> &M, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, int & Z0_r, int & Z0_c, unsigned short & phase );

                static void bestBipartiteMatchingPhase5 ( Matrix <unsigned short> &M, Matrix <unsigned short> &P, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, int & Z0_r, int & Z0_c, unsigned short & phase );

                template <typename T>
                static void bestBipartiteMatchingPhase6 ( Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C, unsigned short & phase );

                template <typename T>
                static void findFirstUncoveredZero ( const Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C, int & row, int & col );

                template <typename T>
                static T findMinUncoveredValue ( const Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C );

                static int findStarredZeroInRow ( const Matrix <unsigned short> &M, const int row );
                static int findStarredZeroInColumn ( const Matrix <unsigned short> &M, const int col );
                static int findPrimedZeroInRow ( const Matrix <unsigned short> &M, const int row );
                static void clearAllCovers ( Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C );
                static void convertPath ( Matrix <unsigned short> &M, const Matrix <unsigned short> &P, const int count );
                static void erasePrimedZeroes ( Matrix <unsigned short> &M );
};

#include "GraphAlgorithms.hpp"

#endif
