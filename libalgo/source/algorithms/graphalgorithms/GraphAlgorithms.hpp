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


#ifndef GraphAlgorithms_HPP
#define GraphAlgorithms_HPP

#include <algorithm>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/tree/KDTree.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"
#include "libalgo/source/algorithms/pointfaceposition/PointFacePosition.h"
#include "libalgo/source/algorithms/linelineposition/LineLinePosition.h"
#include "libalgo/source/algorithms/faceperimeter/FacePerimeter.h"
#include "libalgo/source/algorithms/matrixoperations/MatrixOperations.h"

#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByY.h"

#include "libalgo/source/exceptions/MathMatrixNotSquareException.h"


//Set namespace
using namespace MatrixOperations;


template <typename T>
void GraphAlgorithms::floydWarshall ( const GraphM <T> &g, Matrix <T> &D, Matrix <unsigned int> &P )
{
        //Perform Floyd-Warshall algorithm
        const unsigned int m1 = D.rows();
        const unsigned int n1 = D.cols();
        const unsigned int m2 = P.rows();
        const unsigned int n2 = P.cols();

        //Matrices are not squares
        if ( m1 != n1 || m2 != n2 )
        {
                throw MathMatrixNotSquareException <Matrix <T> > ( "MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not perform Floyd-Warshall algorithm; (rows_row_index, columns_row_index):  ", D );
        }

        //Different matrices size
        if ( m1 != m2 )
        {
                throw MathMatrixDifferentSizeException <Matrix <T> > ( "MathMatrixDifferentSizeException: ", " invalid dimension of the matrices, can not perform Floyd-Warshall algorithm; (rows_row_index columns_row_index):  ", D, P );
        }

        //Initialization D = W
        D = g.getW();

        //Create matrix of predecessors P and set distances to infinity
        for ( unsigned int i = 0; i < m1; i++ )
        {
                for ( unsigned int j = 0; j < m1; j++ )
                {
                        //Non diagonal vertices connected by edges
                        if ( ( i != j ) && ( D ( i, j ) < MAX_FLOAT ) )
                        {
                                P ( i, j ) = i + 1;
                        }

                        //Diagonal items
                        else
                        {
                                P ( i, j ) = 0;
                        }
                }
        }

        //Floyd-Warshall algorithm
        for ( unsigned int k = 0; k < m1; k++ )
        {
                for ( unsigned int i = 0; i < m1; i++ )
                {
                        for ( unsigned int j = 0; j < m1; j++ )
                        {
                                //Test Euclidian distance
                                if ( D ( i, j ) > D ( i, k ) + D ( k, j ) )
                                {
                                        //Calculate new distances
                                        D ( i, j ) = D ( i, k ) + D ( k, j );

                                        //Set predecessors
                                        P ( i, j ) = P ( k, j );
                                }
                        }
                }
        }
}


template <typename T>
GraphL <T> GraphAlgorithms::mst ( GraphL <T> & g, T & weight )
{
        //Minimum spanning tree by Kruskal algorithm
        const unsigned int n_vertices = g.getVertices().size(),
                           n_edges = g.getEdges().size();
        GraphL <T> mst_tree;

        //Set weight to zero
        weight = 0;

        //Sort vertices according to vertex number
        sort ( g.getVertices().begin(), g.getVertices().end() );

        //Create vector of parents
        std::vector <unsigned int> parents ( g.getVertices() [n_vertices - 1] + 1 );

        for ( unsigned int  i = 0; i < n_vertices; i++ )
        {
                parents[g.getVertices() [i]] = g.getVertices() [i];
        }

        //Sort edges according to weights
        std::sort ( g.getEdges().begin(), g.getEdges().end() );

        //Process all edges according to weight
        for ( unsigned int i = 0; i < n_edges; i++ )
        {
                //Find set u
                unsigned int fs_u = findSet ( g.getEdges() [i].getStartPoint(), parents );

                //Find set v
                unsigned int fs_v = findSet ( g.getEdges() [i].getEndPoint(), parents );

                //fs_u != fs_v, add edge to the tree
                if ( fs_u != fs_v )
                {
                        //Compute new cost
                        weight += g.getEdges() [i].getWeight();

                        //Add to the mst tree
                        mst_tree.getEdges().push_back ( g.getEdges() [i] );

                        //Merge subtrees: set the same root vertex for both subtrees
                        parents[fs_u] = parents[fs_v];
                }
        }

        //Return mst tree
        return mst_tree;
}


inline unsigned int GraphAlgorithms::findSet ( const unsigned int x, std::vector <unsigned int> &parents )
{
        //Find root vertex using the recursive approach
        if ( x != parents[x] )
        {
                parents[x] = findSet ( parents[x], parents );
        }

        return parents[x];
}


template <typename T>
void GraphAlgorithms::bestBipartiteMatching ( Matrix <T> C, Matrix <unsigned short> &M, T & cost )
{
        //Best bipartite matching using Hungarian algorihm (Kuhn-Munkres)
        const unsigned int m = C.rows(),
                           n = C.cols();

        //Find max value in cost matrix
        const T max_val = MatrixOperations::max ( C );

        //Get max dimension
        const unsigned int dim = ( std::max ) ( m, n );

        //Copy of C
        Matrix <T> CR ( dim, dim, max_val );

        //Create mask matrix
        Matrix <unsigned short> MR ( dim, dim );

        //Create other matrices
        Matrix <unsigned short> P ( 2 * dim, 2 );
        Matrix <unsigned short> R_C ( dim, 1 );	//Covered rows
        Matrix <unsigned short> C_C ( 1, dim ); //Covered columns

        //Copy C to CR
        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        CR ( i, j ) = C ( i, j );
                }
        }

        //Set minimum rows as terminate condition
        const unsigned int min_rows = dim;

        //Perform repeatdly the following six phases
        unsigned short phase = 1;
        int Z0_r = 0, Z0_c = 0;

        for ( ; ; )
        {
                if ( phase == 1 )
                {
                        bestBipartiteMatchingPhase1 ( CR, phase );
                }

                else if ( phase == 2 )
                {
                        bestBipartiteMatchingPhase2 ( CR, M, R_C, C_C, phase );
                }

                else if ( phase == 3 )
                {
                        bestBipartiteMatchingPhase3 ( MR, C_C, min_rows, phase );
                }

                else if ( phase == 4 )
                {
                        bestBipartiteMatchingPhase4 ( CR, MR, R_C, C_C, Z0_r, Z0_c, phase );
                }

                else if ( phase == 5 )
                {
                        bestBipartiteMatchingPhase5 ( MR, P, R_C, C_C, Z0_r, Z0_c, phase );
                }

                else if ( phase == 6 )
                {
                        bestBipartiteMatchingPhase6 ( CR, R_C, C_C,  phase );
                }

                else
                {
                        break;
                }
        }

        //Trim matching matrix to rectangular (if necessary)
        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int  j = 0; j < n; j++ )
                {
                        M ( i, j ) = MR ( i, j );
                }
        }

        //Calculate cost of the matching
        cost = 0;

        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int  j = 0; j < n; j++ )
                {
                        if ( M ( i, j ) == 1 )
                        {
                                cost += C ( i, j );
                        }
                }
        }
}


template <typename Point>
void GraphAlgorithms::createKNNGraph ( const Container <Point *> *points, const unsigned int k, GraphM <typename Point::Type> &g )
{
        //Create K-NN point graph of the dataset in matrix representation
        const unsigned int n = points->size();

        const typename Point::Type dx = ( * std::max_element ( points->begin(), points->end(), sortPointsByX () ) )->getX()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByX () ) )->getX();
        const typename Point::Type dy = ( * std::max_element ( points->begin(), points->end(), sortPointsByY () ) )->getY()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByY () ) )->getY();
        const unsigned int min_point_id = ( * std::min_element ( points->begin(), points->end(), sortPointsByID <Point>() ) )->getPointID();

        const typename Point::Type norm =  std::max ( dx, dy );

        //Get matrices
        Matrix <int> V ( n, n );
        Matrix <typename Point::Type> W ( n, n );

        //Create KD tree
        KDTree2D <Point> tree;
        tree.createKDTree2D ( points );

        //Process all points
        for ( unsigned int i = 0; i < points ->size(); i++ )
        {
                //All k-nearest neighbours container, non-destructable
                Container <Point *, NonDestructable > knn;

                //Find all k-narest neighbours
                tree.findAllKNN ( ( *points ) [i], &knn, k );

                //Create edges of the graph
                for ( unsigned int j = 1; j < k; j++ )
                {
                        const unsigned int j_nn = knn[j]->getPointID() - min_point_id;
                        typename Point::Type dist = EuclDistance::getEuclDistance2D ( ( *points ) [i], ( *points ) [j_nn] );

                        //Set matrix V
                        V ( i, j_nn ) = 1;

                        //Set matrix W
                        W ( i, j_nn ) = dist / norm;
                }
        }

        //Set results
        g.setV ( V );
        g.setW ( W );
}


template <typename Point>
void GraphAlgorithms::createNNNGraph ( const Container <Point *> *points, GraphM <typename Point::Type> &g )
{
        //Create Natural NN graph of the dataset in matrix representation, every w-distannce is weighted by the
        //length of the twin edge shared by 2 adjacent Voronoi cells
        const unsigned int n = points->size();

        const typename Point::Type dx = ( * std::max_element ( points->begin(), points->end(), sortPointsByX () ) )->getX()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByX () ) )->getX();
        const typename Point::Type dy = ( * std::max_element ( points->begin(), points->end(), sortPointsByY () ) )->getY()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByY () ) )->getY();
        const unsigned int min_point_id = ( * std::min_element ( points->begin(), points->end(), sortPointsByID <Point *> () ) )->getPointID();

        const typename Point::Type norm =  std::max ( dx, dy );

        //Get matrices
        Matrix <int> V ( n, n );
        Matrix <typename Point::Type> W ( n, n );

        //Process all Voronoi cells
        typename TItemsList <Point *>::Type ::const_iterator i_points = points->begin();

        for ( unsigned int i = 0; i_points != points->end() ; i_points ++, i++ )
        {
                //Get actual Voronoi cell
                VoronoiCell <typename Point::Type> *vor_cell = dynamic_cast < VoronoiCell <typename Point::Type> * > ( ( *i_points ) -> getFace() );

                //Voronoi cell exists and it is bounded
                if ( ( vor_cell != NULL ) && ( vor_cell->getBounded() ) )
                {
                        //Get half edge of the Voronoi cell
                        const HalfEdge <typename Point::Type> * e_start = vor_cell->getHalfEdge();
                        HalfEdge <typename Point::Type> * e = const_cast <HalfEdge <typename Point::Type> * > ( e_start );

                        //Find all neighbouring Voronoi cells and compute graph
                        do
                        {
                                //There is neighbouring edge of enough length (and maybe Voronoi cell)
                                if ( ( e->getTwinEdge() != NULL ) && ( EuclDistance::getEuclDistance2D ( e->getPoint(), e->getNextEdge()-> getPoint() ) > MIN_POSITION_DIFF ) )
                                {
                                        //Get neighbouring Voronoi cell
                                        VoronoiCell <typename Point::Type> *vor_cell_adjacent = dynamic_cast < VoronoiCell <typename Point::Type > * > ( ( e -> getTwinEdge()->getFace() ) );

                                        //There is neighbouring and bounded Voronoi cell
                                        if ( ( vor_cell_adjacent != NULL ) && ( vor_cell_adjacent->getBounded() ) )
                                        {
                                                const unsigned int j_nnn = vor_cell_adjacent->getGenerator()->getPointID() - min_point_id;

                                                //Compute distance beween 2 generators: actual Voronoi cell and adjacent Voronoi cell
                                                const typename Point::Type dist_gen_gen = EuclDistance::getEuclDistance2D ( ( *points ) [i], ( *points ) [j_nnn] );

                                                //Compute length of the twin edge
                                                //const typename Point::Type e_twin_length = EuclDistance::getEuclDistance2D(e->getTwinEdge()->getPoint(), e->getTwinEdge()->getNextEdge()->getPoint());

                                                //Set matrix V
                                                V ( i, j_nnn ) = 1;

                                                //Set matrix W
                                                W ( i, j_nnn ) = dist_gen_gen /** e_twin_length */ / ( norm /** norm*/ );
                                        }
                                }

                                e = e->getNextEdge();
                        }
                        while ( e != e_start );
                }
        }

        //Set results
        g.setV ( V );
        g.setW ( W );
}


template <typename Point>
void GraphAlgorithms::createGabrielGraph ( const Container <Point *> *points, GraphM <typename Point::Type> &g )
{
        //Create Gabriel graph of the dataset in matrix representation
        const unsigned int n = points->size();

        const typename Point::Type dx = ( * std::max_element ( points->begin(), points->end(), sortPointsByX () ) )->getX()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByX () ) )->getX();
        const typename Point::Type dy = ( * std::max_element ( points->begin(), points->end(), sortPointsByY () ) )->getY()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByY () ) )->getY();
        const typename Point::Type norm =  std::max ( dx, dy );


        //Get matrices
        Matrix <int> V ( n, n );
        Matrix <typename Point::Type> W ( n, n );

        for ( unsigned int i = 0; i < n; i++ )
        {
                for ( unsigned int j = i + 1; j < n; j++ )
                {
                        //Test edge (points[i], points[j])
                        typename Point::Type xc = 0.5 * ( ( *points ) [i]->getX() + ( *points ) [j]->getX() );
                        typename Point::Type yc = 0.5 * ( ( *points ) [i]->getY() + ( *points ) [j]->getY() );
                        typename Point::Type r = 0.5 * EuclDistance::getEuclDistance2D ( ( *points ) [i], ( *points ) [j] );

                        //Empty circle test
                        bool empty_circle = true;

                        for ( int k = 0; k < n; k++ )
                        {
                                if ( ( k != i ) && ( k != j ) )
                                {
                                        //There is a point of the dataset inside the circle
                                        if ( EuclDistance::getEuclDistance ( xc, yc, ( typename Point::Type ) 0, ( *points ) [k]->getX(), ( *points ) [k]->getY(), ( typename Point::Type ) 0 ) < r )
                                        {
                                                empty_circle = false;
                                                break;
                                        }
                                }
                        }

                        //Add edge to the graph only if a circle is empty
                        if ( empty_circle )
                        {
                                //Set matrix V
                                V ( i, j ) = 1;
                                V ( j, i ) = 1;

                                //Set matrix W
                                W ( i, j ) = 2 * r / norm;
                                W ( j, i ) = W ( i, j ) / norm;
                        }
                }
        }

        //Set results
        g.setV ( V );
        g.setW ( W );
}


template <typename Point>
void GraphAlgorithms::createSphereOfInfulenceGraph ( const Container <Point* > *points, GraphM <typename Point::Type> &g )
{
        //Create Sphere of influence graph of the dataset in matrix representation
        const unsigned int n = points->size();

        const typename Point::Type dx = ( * std::max_element ( points->begin(), points->end(), sortPointsByX () ) )->getX()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByX () ) )->getX();
        const typename Point::Type dy = ( * std::max_element ( points->begin(), points->end(), sortPointsByY () ) )->getY()  -
                                        ( * std::min_element ( points->begin(), points->end(), sortPointsByY () ) )->getY();
        const typename Point::Type norm =  std::max ( dx, dy );

        //Get matrices
        Matrix <int> V ( n, n );
        Matrix <typename Point::Type> W ( n, n );

        //Create KD tree
        KDTree2D <Point> tree;
        tree.createKDTree2D ( const_cast < Container <Point* > *> ( points ) );

        //Find nearest neighbour to each point from the dataset
        std::vector <Point *> knn;

        for ( unsigned int i = 0; i < n; i++ )
        {
                //Find nearest neighbour point to the point
                Container <Point *, NonDestructable > knn_temp;
                tree.findAllKNN ( ( *points ) [i], &knn_temp, 2 );
                knn.push_back ( knn_temp[1] );
        }

        //Process all points
        for ( unsigned int i = 0; i < n; i++ )
        {
                //Find NN to points[i]
                Point *  nn1 = knn[i];

                //Sphere s1(points[i], d(points[i], nn1))
                typename Point::Type r1 = EuclDistance::getEuclDistance2D ( ( *points ) [i], nn1 );

                for ( unsigned int j = i + 1; j < n; j++ )
                {
                        //Find NN to points[j]
                        Point * nn2 = knn[j];

                        //Shere s2(points[j], d(points[j], nn2))
                        typename Point::Type r2 = EuclDistance::getEuclDistance2D ( ( *points ) [j], nn2 );

                        //Spheres s1 and s2 overlap
                        typename Point::Type dist = EuclDistance::getEuclDistance2D ( ( *points ) [i], ( *points ) [j] );

                        if ( dist + ARGUMENT_ROUND_ERROR < r1 + r2 )
                        {
                                //Set matrix V
                                V ( i, j ) = 1;
                                V ( j, i ) = 1;

                                //Set matrix W
                                W ( i, j ) = dist / norm;
                                W ( j, i ) = dist / norm;
                        }
                }
        }

        //Set results
        g.setV ( V );
        g.setW ( W );
}



template <typename T>
void GraphAlgorithms::FaceToWDistanceGraph ( const Face <T> *f, GraphM <T> &g, const bool normalized )
{
        //Convert face to complete graph
        //	normalized = false: w-distances are unnormalized
        //	normalized = true: w-distances are normalized
        //Test, if face is simple
        if ( f != NULL )
        {
                T area = 0;
                const HalfEdge <T> *e_start = f->getHalfEdge();

                //There is a valid edge
                if ( e_start != NULL )
                {
                        //Create graph matrices
                        unsigned int n = f->getVerticesCount();
                        Matrix <unsigned int> V ( n, n, 0 );
                        Matrix <T> W ( n, n, 0 );

                        //Get perimeter of the Face (if normalized)
                        T perimeter = 1;

                        if ( normalized )
                        {
                                perimeter = FacePerimeter::getFacePerimeter ( f );
                        }

                        //Proces all edges of the Face
                        HalfEdge <T> *e1 = const_cast <HalfEdge <T> *> ( e_start );
                        unsigned int i = 0;

                        do
                        {
                                //Get start point of the first edge ( tested, if inner / outer )
                                const Node3DCartesian <T> *p1 = e1->getPoint();

                                //Proces all edges of the Face
                                HalfEdge <T> *e2 = e1->getNextEdge();

                                unsigned int j = ( e2 == e_start ? 0 : i + 1 );

                                do
                                {
                                        //Get start point of the first edge ( tested, if inner / outer )
                                        const Node3DCartesian <T> *p2 = e2->getPoint();

                                        //The first segment does not have to be too short
                                        T dist_p1_p2 =  EuclDistance::getEuclDistance2D ( p1, p2 );

                                        if ( dist_p1_p2 > MIN_POSITION_DIFF )
                                        {
                                                V ( i, j ) = 1;
                                                V ( j, i ) = 1;
                                                W ( i, j ) = dist_p1_p2 / perimeter;
                                                W ( j, i ) = W ( i, j );
                                        }

                                        //Increment j
                                        j++;

                                        //Increment edge
                                        e2 = e2->getNextEdge();

                                }
                                while ( e2 != e_start );

                                //Increment i
                                i++;

                                //Increment edge
                                e1 = e1->getNextEdge();

                        }
                        while ( e1 != e_start );

                        //Add matrices to the graph
                        g.setV ( V );
                        g.setW ( W );
                }
        }
}


template <typename T>
void GraphAlgorithms::FaceToInnerDistanceGraph ( const Face <T> *f, GraphM <T> &g, const bool normalized )
{
        //Convert face to inner distance graph
        //	normalized = false: w-distances are unnormalized
        //	normalized = true: w-distances are normalized

        //Test, if face is simple
        if ( f != NULL )
        {
                const HalfEdge <T> *e_start = f->getHalfEdge();

                //There is a valid edge
                if ( e_start != NULL )
                {
                        //Create graph matrices
                        unsigned int n = f->getVerticesCount();
                        Matrix <unsigned int> V ( n, n );
                        Matrix <T> W ( n, n, MAX_FLOAT, 0 );

                        //Get perimeter of the Face (if normalized)
                        T perimeter = 1;

                        if ( normalized )
                        {
                                perimeter = FacePerimeter::getFacePerimeter ( f );
                        }

                        //Proces all edges of the Face
                        HalfEdge <T> *e1 = const_cast <HalfEdge <T> *> ( e_start );
                        unsigned int i = 0;

                        do
                        {
                                //Get start point of the first inner/outer edge
                                const Node3DCartesian <T> *p1 = e1->getPoint();

                                //Do not compute edge starts FROM duplicate point (i.e from second point with same coordinates)
                                const T dist_p1_p1_prev =  EuclDistance::getEuclDistance2D ( p1, e1->getPreviousEdge()->getPoint() );

                                if ( dist_p1_p1_prev > MIN_POSITION_DIFF )
                                {
                                        //Proces all edges of the Face
                                        HalfEdge <T> *e2 = e1->getNextEdge();
                                        unsigned int j = ( e2 == e_start ? 0 : i + 1 );

                                        while ( e2 != e_start )
                                        {
                                                //Get end point of the first  inner/outer edge
                                                const Node3DCartesian <T> *p2 = e2->getPoint();

                                                //Compute inner / outer edge length and test, and do not compute edge end AT duplicate point
                                                const T dist_p1_p2 =  EuclDistance::getEuclDistance2D ( p1, p2 );
                                                const T dist_p2_p2_prev =  EuclDistance::getEuclDistance2D ( p2, e2->getPreviousEdge()->getPoint() );

                                                //Inner/outer does not have to be too short and end point does not to be duplica compute weight matrix FROM duplicate point (i.e from second point with same coordinates)
                                                if ( ( dist_p1_p2 > MIN_POSITION_DIFF ) && ( dist_p2_p2_prev > MIN_POSITION_DIFF ) )
                                                {
                                                        //Test intersections of the inner/outer edge with the segments of the face
                                                        HalfEdge <T> *e3 = const_cast <HalfEdge <T> *> ( e_start );

                                                        //Get start point of the face edge
                                                        const Node3DCartesian <T> *p3 = e3->getPoint();

                                                        //Tested inner/outer edge does not have incident vertices with start point of the the face edge
                                                        bool intersection_exists = false;

                                                        do
                                                        {
                                                                //Get end point of the face edge
                                                                const Node3DCartesian <T> *p4 = e3->getNextEdge()-> getPoint();

                                                                //Face edge does not have to be too short
                                                                if ( EuclDistance::getEuclDistance2D ( p3, p4 ) > MIN_POSITION_DIFF )
                                                                {
                                                                        //Tested inner/outer edge does not have incident vertices with end points of the the face edge
                                                                        if ( ( *p1 != *p3 ) && ( *p2 != *p3 ) && ( *p1 != *p4 ) && ( *p2 != *p4 ) )
                                                                        {
                                                                                //Test intersection of the inner / outer edge with the actual face edge
                                                                                double x_int, y_int;
                                                                                unsigned short t = LineLinePosition::get2LineSegmentsPosition ( p1, p2, p3, p4, x_int, y_int );

                                                                                //There is intersecton: stop checking
                                                                                if ( ( t > 0 ) && ( t < 4 ) )
                                                                                {
                                                                                        intersection_exists = true;
                                                                                        break;
                                                                                }
                                                                        }
                                                                }

                                                                //Assign point
                                                                p3 = p4;

                                                                //Increment edge
                                                                e3 = e3->getNextEdge();

                                                        }
                                                        while ( e3 != e_start ) ;


                                                        //Inner or outer edge
                                                        if ( !intersection_exists )
                                                        {
                                                                //Get mid-point of the tested edge
                                                                Node3DCartesian <T>  mid ( ( p1->getX() + p2->getX() ) / 2.0, ( p1 ->getY() + p2->getY() ) / 2.0 );

                                                                //Test position of the mid-point and face
                                                                unsigned short t = PointFacePosition::getPointFacePosition ( &mid, f );

                                                                //Outer edge: mid point outside polygon
                                                                if ( t != 1 )
                                                                {
                                                                        V ( i, j ) = 1;
                                                                        V ( j, i ) = 1;
                                                                        W ( i, j ) = dist_p1_p2 / perimeter;
                                                                        W ( j, i ) = W ( i, j );
                                                                }
                                                        }
                                                }

                                                //Increment index j
                                                j++;

                                                //Increment edge
                                                e2 = e2->getNextEdge();
                                        }
                                }

                                //Increment index i
                                i++;

                                //Increment edge
                                e1 = e1->getNextEdge();

                        }
                        while ( e1 != e_start );

                        //Add matrices to the graph
                        g.setV ( V );
                        g.setW ( W );
                }
        }
}


template <typename T>
void GraphAlgorithms::bestBipartiteMatchingPhase1 ( Matrix <T> &CR, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm: phase 1, subtract min from each  item of the row
        const unsigned int m = CR.rows();

        //Find min in each row
        for ( unsigned int i = 0; i < m ; i++ )
        {
                //Initialize min
                T min = CR ( i, 0 );

                //Process row, find min
                for ( unsigned int j = 0; j < m; j++ )
                {
                        if ( CR ( i, j ) < min )
                        {
                                min = CR ( i, j );
                        }
                }

                //Subtract min from each item of the row
                for ( unsigned int j = 0 ; j < m; j++ )
                {
                        CR ( i, j ) = CR ( i, j ) - min;
                }
        }

        //Set phase
        phase = 2;
}


template <typename T>
void GraphAlgorithms::bestBipartiteMatchingPhase2 ( const Matrix < T> &C, Matrix <unsigned short> &M, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm: phase 2, fonf no starred zeros and set row and col as starred
        const unsigned int m = M.rows();;

        //Find zero items
        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int j = 0; j < m; j++ )
                {
                        //Found no starred zero, set row and col containing zero as stared
                        if ( ( C ( i, j ) == 0 ) && ( R_C ( i, 0 ) == 0 ) && ( C_C ( 0, j ) == 0 ) )
                        {
                                M ( i, j ) = 1;
                                R_C ( i, 0 ) = 1;
                                C_C ( 0, j ) = 1;
                        }
                }
        }

        //Set zero for both matrices
        for ( unsigned int i = 0; i < m; i++ )
        {
                R_C ( i, 0 ) = 0;
                C_C ( 0, i ) = 0;
        }

        //Set phase
        phase = 3;
}


inline void GraphAlgorithms::bestBipartiteMatchingPhase3 ( const Matrix <unsigned short> &M, Matrix <unsigned short> &C_C, unsigned int min_rows, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm: terminate or continue
        const unsigned int m = M.rows();

        //Get covering
        for ( unsigned int i = 0; i < m; i++ )
        {
                for ( unsigned int j = 0; j < m; j++ )
                {
                        if ( M ( i, j ) == 1 )
                        {
                                C_C ( 0, j ) = 1;
                        }
                }
        }

        //Check sum of C_C covers
        unsigned int sum = 0;

        for ( unsigned int j = 0; j < m; j++ )
        {
                sum += C_C ( 0, j );
        }

        //C_C  covers >= min_rows: terminate computation
        if ( sum >= min_rows )
        {
                phase = 7;
        }

        //Continue in computation
        else
        {
                phase = 4;
        }
}


template <typename T>
void GraphAlgorithms::bestBipartiteMatchingPhase4 ( const Matrix <T> &C, Matrix <unsigned short> &M, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, int & Z0_r, int & Z0_c, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm: find a noncovered zeros and prime them.
        //const unsigned int n = M.cols();

        //Find non covered zeros
        for ( ; ; )
        {
                int row = -1, col = -1;

                //Find row of the first non covered zero
                findFirstUncoveredZero ( C, R_C, C_C, row, col );

                //No row has been found
                if ( row == -1 )
                {
                        //Continue using 6
                        phase = 6;

                        //Set Z0_r, Z0_c as zero
                        Z0_r = 0;
                        Z0_c = 0;

                        //Stop
                        break;
                }

                //Row has been found
                else
                {
                        //Prime zero
                        M ( row, col ) = 2;

                        //Find stared zero in actual row
                        const int zero_column = findStarredZeroInRow ( M, row );

                        if ( zero_column >= 0 )
                        {
                                //Cover this row
                                R_C ( row, 0 ) = 1;

                                //Uncover this col
                                C_C ( 0, zero_column ) = 0;
                        }

                        //Stared zero has not been found
                        else
                        {
                                Z0_r = row;
                                Z0_c = col;

                                //Continue using 5
                                phase = 5;

                                //Stop
                                break;
                        }
                }
        }
}


inline void GraphAlgorithms::bestBipartiteMatchingPhase5 ( Matrix <unsigned short> &M, Matrix <unsigned short> &P, Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C, int & Z0_r, int & Z0_c, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm:  find col of primed zero in starred zero row,
        unsigned int row_index = 0;

        //Initialize path matrix
        P ( row_index, 0 ) = Z0_r;
        P ( row_index, 1 ) = Z0_c;

        //Find col of primed zero in starred zero row
        for ( ; ; )
        {
                //Find starred zero in row
                const int row = findStarredZeroInColumn ( M, P ( row_index, 1 ) );

                //No row has been found
                if ( row < 0 )
                {
                        break;
                }

                //Row has been found
                else
                {
                        //Increment row index
                        row_index++;

                        //Remember row of the starred zero
                        P ( row_index, 0 ) = row;

                        //Assign index
                        P ( row_index, 1 ) = P ( row_index - 1, 1 );

                        //Find col of primed zero in starred zero row
                        const int col = findPrimedZeroInRow ( M, P ( row_index, 0 ) );

                        //No col has been found
                        //if ( col < 0 )
                        {
                                //        break;
                        }

                        //Increment row index
                        row_index ++;

                        //Remember col
                        P ( row_index, 1 ) = col;

                        //Assign index
                        P ( row_index, 0 ) = P ( row_index - 1, 0 );
                }
        }

        //Set starred zeros in the path as unstarred, set all primed zeros as starred
        convertPath ( M, P, row_index );

        //Clear all covers
        clearAllCovers ( R_C, C_C );

        //Remove all primes
        erasePrimedZeroes ( M );

        //Set phase to 3
        phase = 3;
}


template <typename T>
void GraphAlgorithms::bestBipartiteMatchingPhase6 ( Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C, unsigned short & phase )
{
        //Best bipartite matching using Hungarian algorihm: Add the minimum uncovered value to every element of covered row, subtract  the minimum uncovered value from every element of each uncovered col.
        const unsigned int n = C.cols();

        //Find minimum value
        T min  = findMinUncoveredValue ( C, R_C, C_C );

        //Process all items
        for ( unsigned int i = 0; i < n; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        //Add the minimum uncovered value to every element of covered row
                        if ( R_C ( i, 0 ) == 1 )
                        {
                                C ( i, j ) += min;
                        }

                        //Subtract  the minimum uncovered value from every element of each uncovered col.
                        if ( C_C ( 0, j ) == 0 )
                        {
                                C ( i, j ) -= min;
                        }
                }
        }

        //Set phase to 4
        phase = 4;
}


template <typename T>
void GraphAlgorithms::findFirstUncoveredZero ( const Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C, int & row, int & col )
{

        //Find row of the first non covered zero
        const unsigned int n = C.cols();

        for ( unsigned int i = 0; i < n; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        if ( ( C ( i, j ) == 0 ) && ( R_C ( i, 0 ) == 0 ) && ( C_C ( 0, j ) == 0 ) )
                        {
                                row = i;
                                col = j;

                                //Was found, break
                                return;
                        }
                }
        }
}


template <typename T>
T GraphAlgorithms::findMinUncoveredValue ( const Matrix <T> &C, const Matrix <unsigned short> &R_C, const Matrix <unsigned short> &C_C )
{
        //Find minimum uncovered value in C
        const unsigned int n = C.cols();

        //Initialize min
        T min = MAX_FLOAT;

        //Process all items
        for ( unsigned int i = 0; i < n; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        //Unvovered value
                        if ( ( R_C ( i, 0 ) == 0 ) && ( C_C ( 0, j ) == 0 ) )
                        {
                                if ( min > C ( i, j ) )
                                {
                                        min = C ( i, j );
                                }
                        }
                }
        }

        return min;
}


inline int GraphAlgorithms::findStarredZeroInRow ( const Matrix <unsigned short> &M, const int row )
{
        //Find col of the stared zero in the row
        const unsigned int n = M.cols();

        //Process all items in the row
        for ( unsigned int j = 0; j < n; j++ )
        {
                if ( M ( row, j ) == 1 )
                {
                        return j;
                }
        }

        return -1;
}


inline int GraphAlgorithms::findStarredZeroInColumn ( const Matrix <unsigned short> &M, const int col )
{
        //Find row of the stared zero in the col
        unsigned int n = M.cols();

        //Process all items in the col
        for ( unsigned int i = 0; i < n; i++ )
        {
                if ( M ( i, col ) == 1 )
                {
                        return i;
                }
        }

        return -1;
}


inline int GraphAlgorithms::findPrimedZeroInRow ( const Matrix <unsigned short> &M, const int row )
{
        //Find col of the primed zero in the row
        const unsigned int n = M.cols();

        //Process all items in the row
        for ( unsigned int j = 0; j < n; j++ )
        {
                if ( M ( row, j ) == 2 )
                {
                        return j;
                }
        }

        return -1;
}


inline void GraphAlgorithms::clearAllCovers ( Matrix <unsigned short> &R_C, Matrix <unsigned short> &C_C )
{
        //Clear all covers
        const unsigned int n = C_C.cols();

        for ( unsigned int i = 0; i < n; i++ )
        {
                R_C ( i, 0 ) = 0;
                C_C ( 0, i ) = 0;
        }
}


inline void GraphAlgorithms::convertPath ( Matrix <unsigned short> &M, const Matrix <unsigned short> &P, const int row_index )
{
        //Set all the starred zeros in the path  to unstarred and set all primed zeros to starred
        for ( unsigned int i = 0; i <= ( unsigned int ) row_index; i++ )
        {
                //Set starred zeros to  unstarred
                if ( M ( P ( i, 0 ), P ( i, 1 ) ) == 1 )
                {
                        M ( P ( i, 0 ), P ( i, 1 ) ) = 0;
                }

                //Set primed zeroes to starred
                else
                {
                        M ( P ( i, 0 ), P ( i, 1 ) ) = 1;
                }
        }
}


inline void GraphAlgorithms::erasePrimedZeroes ( Matrix <unsigned short> &M )
{
        //Erase all primes
        const unsigned int n = M.cols();

        //Process all items
        for ( unsigned int i = 0; i < n; i++ )
        {
                for ( unsigned int j = 0; j < n; j++ )
                {
                        if ( M ( i, j ) == 2 )
                        {
                                M ( i, j ) = 0;
                        }
                }
        }
}


#endif
