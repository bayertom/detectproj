// Description: 2D Delaunay triangle by incremental insertion + Lawson Oriented Walk, topological model, de Berg et al, 2000

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


#ifndef DT2D_HPP
#define DT2D_HPP

#include <vector>
#include <list>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "libalgo/source/structures/line/HalfEdge.h"

#include "libalgo/source/algorithms/lawsonorientedwalk/LawsonOrientedWalk.h"
#include "libalgo/source/algorithms/convexquadrilateral/ConvexQuadrilateral.h"
#include "libalgo/source/algorithms/swappingcriteria/SwappingCriteria.h"

#include "libalgo/source/comparators/removeSimplexHalfEdge.h"
#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByY.h"
#include "libalgo/source/comparators/sortPointsToBins.h"
#include "libalgo/source/comparators/isEqualPointByPlanarCoordinates.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"
#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
void DT2D::DT ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges_dt, const bool print_message, const bool print_exception, std::ostream * output )
{
        // Create 2D Delaunay triangulation using incremental insertion method
        const unsigned int nodes_count_before = nl.size();

        // Remove duplicate points
        nl.removeDuplicateElements ( nl.begin(), nl.end(), sortPointsByX (), isEqualPointByPlanarCoordinates <Node3DCartesian <T> *> () );

        // Get nodes count after deletion of duplicated points
        unsigned int nodes_count_after = nl.size();

        //Print info
        if ( print_message )
        {
                *output << "> Starting DT, please wait... ";
                *output << nodes_count_after << " points, " << ( nodes_count_before - nodes_count_after ) << " removed.";
        }

        try
        {
                //There are at least 3 points
                if ( nodes_count_after > 2 )
                {
                        //Create simplex triangle
                        createSimplexTriangle ( nl, half_edges_dt );

                        //Increment nodes count
                        nodes_count_after += 3;

                        //Starting half edge using for searching
                        HalfEdge <T> *e_heuristic = half_edges_dt [0];

                        // Insert all points into triangulation using incremental method
                        for ( unsigned int i = 3; i < nodes_count_after; i++ )	// Jump over simplex
                        {
                                DTInsertPoint ( nl [i], nl [0], nl [1], nl [2], &e_heuristic, half_edges_dt );
                        }

                        //Remove triangles having simplex points
                        removeSimplexTriangles ( nl, half_edges_dt );
                }

                //Sort points by ID again
                std::sort ( nl.begin(), nl.end(), sortPointsByID <Node3DCartesian <T> *> () );

                //Print results
                if ( print_message )
                {
                        *output << " Completed." << std::endl;
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                //Clear half edges
                half_edges_dt.clear();

                //Print exception
                if ( print_exception )
                {
                        error.printException ( output );
                }

                throw;
        }
}


template <typename T>
void DT2D::createSimplexTriangle ( Container <Node3DCartesian <T> *> &pl, Container <HalfEdge <T> *> &half_edges_dt )
{
        //Create simplex triangle
        Node3DCartesian <T> *s1 = NULL;
        Node3DCartesian <T> *s2 = NULL;
        Node3DCartesian <T> *s3 = NULL;
        HalfEdge <T> *he1 = NULL;
        HalfEdge <T> *he2 = NULL;
        HalfEdge <T> *he3 = NULL;

        //Find extremal nodes
        const Node3DCartesian <T> *p_x_max = * std::max_element ( pl.begin(), pl.end(), sortPointsByX () );
        const Node3DCartesian <T> *p_x_min = * std::min_element ( pl.begin(), pl.end(), sortPointsByX () );
        const Node3DCartesian <T> *p_y_max = * std::max_element ( pl.begin(), pl.end(), sortPointsByY () );
        const Node3DCartesian <T> *p_y_min = * std::min_element ( pl.begin(), pl.end(), sortPointsByY () );

        // Find width and height of the min-max box
        const T	dx = p_x_max->getX() - p_x_min->getX();
        const T	dy = p_y_max->getY() - p_y_min->getY();

        //All points are colinear
        if ( dx == 0 || dy == 0 )
                throw BadDataException ( "BadDataException: All input points are colinear. ", "Can not contruct Delaunay triangulation. " );

        // Coordinates of center of the min-max box
        const T	xc = 0.5 * ( p_x_max->getX() + p_x_min->getX() );
        const T	yc = 0.5 * ( p_y_max->getY() + p_y_min->getY() );;

        // Create simplex, new  3 simplex point
        s1 = new Node3DCartesian <T> ( ( T ) ( xc + 3 * ( std::max ) ( dx, dy ) ), ( T ) yc );
        s2 = new Node3DCartesian <T> ( ( T ) xc, ( T ) ( yc + 3 * ( std::max ) ( dx, dy ) ) );
        s3 = new Node3DCartesian <T> ( ( T ) ( xc - 3 * ( std::max ) ( dx, dy ) ), ( T ) ( yc - 3 * ( std::max ) ( dx, dy ) ) );

        // Create first simplex triangle (3 half edges)
        he3 = new HalfEdge <T> ( s3, NULL, NULL );
        he2 = new HalfEdge <T> ( s2, he3, NULL );
        he1 = new HalfEdge <T> ( s1, he2, NULL );

        //Link first and last edge of the triangle
        he3->setNextEdge ( he1 );

        // Add 3 simplex points to list
        pl.insert ( pl.begin(), s1 );
        pl.insert ( pl.begin(), s2 );
        pl.insert ( pl.begin(), s3 );

        // Add 3 half edges of the triangle to list
        half_edges_dt.push_back ( he1 );
        half_edges_dt.push_back ( he2 );
        half_edges_dt.push_back ( he3 );

        // Sort points to bin (do not sort 3 simplex points)
        pl.sort ( pl.begin() += 3, pl.end(), sortPointsToBins <Node3DCartesian<T> *> ( p_x_min->getX(), p_y_min->getY(), p_x_max->getX(), p_y_max->getY(), ( unsigned int ) ( 0.1 * ( sqrt ( ( double ) pl.size() ) ) ) ) );
}


template <typename T>
void DT2D::DTInsertPoint ( Node3DCartesian <T> *p, const Node3DCartesian <T> *s1, const Node3DCartesian <T> *s2, const Node3DCartesian <T> *s3, HalfEdge <T> **e1, Container <HalfEdge <T> *>  &half_edges_dt )
{
        // One step of the Delaunay triangulation, incremental insertion by de Berg (2001)
        short status = -1;

        //HalfEdges
        HalfEdge <T> *e2 = NULL;
        HalfEdge <T> *e3 = NULL;
        HalfEdge <T> *e4 = NULL;
        HalfEdge <T> *e5 = NULL;
        HalfEdge <T> *e6 = NULL;
        HalfEdge <T> *e31 = NULL;
        HalfEdge <T> *e21 = NULL;
        HalfEdge <T> *e12 = NULL;
        HalfEdge <T> *e32 = NULL;
        HalfEdge <T> *e23 = NULL;
        HalfEdge <T> *e13 = NULL;
        HalfEdge <T> *e53 = NULL;
        HalfEdge <T> *e44 = NULL;
        HalfEdge <T> *e63 = NULL;

        try
        {
                // Test, if point lies inside triangle
		const T max_steps = 5 * half_edges_dt.size();
                *e1 = LawsonOrientedWalk::findFaceWalk2 ( p, &status, *e1, max_steps);

		//Boundary triangle found, add point
                if ( *e1 != NULL )
                {
                        // Remaining edges of the triangle inside the added point lies
                        e2 = ( *e1 )->getNextEdge();
                        e3 = e2->getNextEdge();

                        // Point lies inside the triangle: split into 3 triangles
                        if ( status == 1 )
                        {
                                // Create first new triangle T1, twin edges set after creation
                                e31 = new HalfEdge <T> ( p, *e1, NULL );
                                e21 = new HalfEdge <T> ( e2->getPoint(), e31, NULL );
                                ( *e1 )->setNextEdge ( e21 );

                                // Create second new triangle T2, twin edges set after creation
                                e12 = new HalfEdge <T> ( p, e2, NULL );
                                e32 = new HalfEdge <T> ( e3->getPoint(), e12, NULL );
                                e2->setNextEdge ( e32 );

                                // Create third new triangle T3, twin edges set after creation
                                e23 = new HalfEdge <T> ( p, e3, NULL );
                                e13 = new HalfEdge <T> ( ( *e1 )->getPoint(), e23, NULL );
                                e3->setNextEdge ( e13 );

                                // Set twin edges in T1, T2, T3
                                e12->setTwinEdge ( e21 );
                                e21->setTwinEdge ( e12 );
                                e13->setTwinEdge ( e31 );
                                e31->setTwinEdge ( e13 );
                                e23->setTwinEdge ( e32 );
                                e32->setTwinEdge ( e23 );

                                // Add new edges into list
                                half_edges_dt.push_back ( e21 );
                                half_edges_dt.push_back ( e12 );
                                half_edges_dt.push_back ( e31 );
                                half_edges_dt.push_back ( e13 );
                                half_edges_dt.push_back ( e32 );
                                half_edges_dt.push_back ( e23 );

                                // Legalize triangle T1
                                if ( ( *e1 )->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, *e1, s1, s2, s3 );
                                }

                                // Legalize triangle T2
                                if ( e2->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e2, s1, s2, s3 );
                                }

                                // Legalize triangle T3
                                if ( e3->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e3, s1, s2, s3 );
                                }
                        }

                        // Point lies on the edge of the triangle: split into 4 triangles
                        else if ( status == 2 )
                        {
                                // Find adjacent triangle
                                e4 = ( *e1 )->getTwinEdge();
                                e5 = e4->getNextEdge();
                                e6 = e5->getNextEdge();

                                // Create first new triangle T1, twin edges set after creation
                                e21 = new HalfEdge < T > ( p, e3, NULL );
                                ( *e1 )->setNextEdge ( e21 );

                                // Create second new triangle T2, OK
                                e12 = new HalfEdge <T> ( p, e2, e4 );
                                e32 = new HalfEdge <T> ( e3->getPoint(), e12, e21 );
                                e2->setNextEdge ( e32 );

                                // Create third new triangle T3, twin edges set after creation
                                e53 = new HalfEdge <T> ( p, e6, NULL );
                                e4->setNextEdge ( e53 );

                                // Create fourth new triangle T4, OK
                                e44 = new HalfEdge <T> ( p, e5, *e1 );
                                e63 = new HalfEdge <T> ( e6->getPoint(), e44, e53 );
                                e5->setNextEdge ( e63 );

                                // Set twin edges in T1, T3
                                e21->setTwinEdge ( e32 );
                                ( *e1 )->setTwinEdge ( e44 );
                                e53->setTwinEdge ( e63 );
                                e4->setTwinEdge ( e12 );

                                // Add new edges into list
                                half_edges_dt.push_back ( e21 );
                                half_edges_dt.push_back ( e12 );
                                half_edges_dt.push_back ( e32 );
                                half_edges_dt.push_back ( e53 );
                                half_edges_dt.push_back ( e63 );
                                half_edges_dt.push_back ( e44 );

                                // Legalize triangle T1
                                if ( e3->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e3, s1, s2, s3 );
                                }

                                // Legalize triangle T4
                                if ( e5->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e5, s1, s2, s3 );
                                }

                                // Legalize triangle T3
                                if ( e6->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e6, s1, s2, s3 );
                                }

                                // Legalize triangle T2
                                if ( e2->getTwinEdge() != NULL )
                                {
                                        legalizeTriangle ( p, e2, s1, s2, s3 );
                                }
                        }
                }

		//Boundary triangle did not find
		else
		{
			//Re-initialize start edge
			*e1 = half_edges_dt[0];
		}
        }

        //Throw exception
        catch ( MathZeroDevisionException <T> & error )
        {
                //Set NULL pointers
                if ( e1 != NULL ) ( *e1 ) -> setNextEdge ( NULL );

                if ( e2 != NULL ) e2 ->setNextEdge ( NULL );

                if ( e3 != NULL ) e3 ->setNextEdge ( NULL );

                if ( e4 != NULL ) e4->setNextEdge ( NULL );

                if ( e5 != NULL ) e5->setNextEdge ( NULL );

                //Free memory
                if ( e31 != NULL ) delete e31;

                if ( e21 != NULL ) delete e21;

                if ( e12 != NULL ) delete e12;

                if ( e32 != NULL ) delete e32;

                if ( e23 != NULL ) delete e23;

                if ( e13 != NULL ) delete e13;

                if ( e53 != NULL ) delete e53;

                if ( e44 != NULL ) delete e44;

                if ( e63 != NULL ) delete e63;

                //Throw exception
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: ", "Delaunay triangulation: Can not create new triangles for inserted point p.", error.getArg() );
        }
}


template <typename T>
void DT2D::legalizeTriangle ( Node3DCartesian <T> *p, HalfEdge <T> *e11, const Node3DCartesian <T> *s1, const Node3DCartesian <T> *s2, const Node3DCartesian <T> *s3 )
{
        //Recursive legalization procedure: defined in Computational Geometry, de Berg et al, pp 196
        HalfEdge <T> *e21 = e11->getTwinEdge();		// First edge in triangle adjacent to legalized
        HalfEdge <T> *e22 = e21->getNextEdge();		// Second edge in triangle adjacent to legalized
        HalfEdge <T> *e23 = e22->getNextEdge();		// Third edge in triangle adjacent to legalized
        HalfEdge <T> *e12 = e11->getNextEdge();         // Second edge in  legalized triangle
        HalfEdge <T> *e13 = e12->getNextEdge();		// Third edge in  legalized triangle

        //Get vertices of the quadrilateral
        const Node3DCartesian <T> *p1 = e21->getPoint();
        const Node3DCartesian <T> *p2 = e22->getPoint();
        const Node3DCartesian <T> *p3 = e23->getPoint();

        //Test, if points p1, p2, p3, p4 are identical to simplex points s1, s2, s3
        //Returns number of point to which is point p identical (1 - 3),
        //We test each vertex  of the quadrilateral to each simplex point
        const unsigned short b11 = ( p1 == s1 );
        const unsigned short b12 = 2 * ( p1 == s2 );
        const unsigned short b13 = 3 * ( p1 == s3 );
        const unsigned short b21 = ( p2 == s1 );
        const unsigned short b22 = 2 * ( p2 == s2 );
        const unsigned short b23 = 3 * ( p2 == s3 );
        const unsigned short b31 = ( p3 == s1 );
        const unsigned short b32 = 2 * ( p3 == s2 );
        const unsigned short b33 = 3 * ( p3 == s3 );
        const unsigned short b41 = ( p == s1 );
        const unsigned short b42 = 2 * ( p == s2 );
        const unsigned short b43 = 3 * ( p == s3 );

        //Point p1 is identical with any of point s1, s2, s3
        const unsigned short b1 = b11 + b12 + b13;

        //Point p2 is identical with any of point s1, s2, s3
        const unsigned short b2 = b21 + b22 + b23;

        //Point p3 is identical with any of point s1, s2, s3
        const unsigned short b3 = b31 + b32 + b33;

        //Point p4 is identical with any of point s1, s2, s3
        const unsigned short b4 = b41 + b42 + b43;

        //Define rules
        bool illegal_triangle = false;

        //Rule 1: (p1 and p2) are simplex points
        if ( ( b1 > 0 ) && ( b2 > 0 ) )
        {
                //Edge is legal, we want to keep the large triangle
                illegal_triangle = false;
        }

        //Rule2: no point is simplex: test Delaunay condition
        else if ( ( b1 == 0 ) && ( b2 == 0 ) && ( b3 == 0 ) && ( b4 == 0 ) )
        {
                //Swap diagonal in convex quadrilateral using Cline-Renka
                illegal_triangle = SwappingCriteria::getClineRenka ( p, p1, p2, p3 ) ;
        }

        //Rule 3: One point is a simplex
        else if ( ( b1 > 0 ) && ( b2 + b3 + b4 == 0 ) || ( b2 > 0 ) && ( b1 + b3 + b4 == 0 ) || ( b3 > 0 ) && ( b1 + b2 + b4 == 0 ) || ( b4 > 0 ) && ( b1 + b2 + b3 == 0 ) )
        {
                //Point p1 or p2 represent a simplex
                if ( b1 > 0 || b2 > 0 )
                {
                        illegal_triangle = ( ConvexQuadrilateral::isStrictlyConvex ( p1, p, p2, p3 ) == 1 );
                }
        }

        //Rule 4 : (p1 or p2) and (p3 or p4) are simplex points
        else if ( ( b1 * b3 > 0 ) || ( b2 * b4 > 0 ) || ( b1 * b4 > 0 ) || ( b2 * b3 > 0 ) )
        {
                //Points p1, p3 are simplex
                if ( b1 * b3 > 0 )
                {
                        //Swap diagonal
                        if ( b1 > b3 )
                        {
                                illegal_triangle = ( ConvexQuadrilateral::isStrictlyConvex ( p1, p, p2, p3 ) == 1 );
                        }
                }

                //Points p2, p3  are simplex
                if ( b2 * b3 > 0 )
                {
                        //Swap diagonal
                        if ( b2 > b3 )
                        {
                                illegal_triangle = ( ConvexQuadrilateral::isStrictlyConvex ( p1, p, p2, p3 ) == 1 );
                        }
                }

                //Points p1, p4  are simplex
                if ( b1 * b4 > 0 )
                {
                        //Swap diagonal
                        if ( b1 > b4 )
                        {
                                illegal_triangle = ( ConvexQuadrilateral::isStrictlyConvex ( p1, p, p2, p3 ) == 1 );
                        }
                }

                //Points p2, p4 are simplex
                if ( b2 * b4 > 0 )
                {
                        //Swap diagonal
                        if ( b2 > b4 )
                        {
                                illegal_triangle = ( ConvexQuadrilateral::isStrictlyConvex ( p1, p, p2, p3 ) == 1 );
                        }
                }
        }

        // Illegal triangle, swap diagonal
        if ( illegal_triangle )
        {
                // Swap diagonal in (e11,e12,e13) and (e21,e22,e23)
                swapDiagonal ( e11, e12, e13, e21, e22, e23 );

                // Legalize adjacent triangles
                if ( e22->getTwinEdge() != NULL )
                {
                        legalizeTriangle ( p, e22, s1, s2, s3 );			// Triangle adjacent to first newly created
                }

                if ( e23->getTwinEdge() != NULL )
                {
                        legalizeTriangle ( p, e23, s1, s2, s3 );			// Triangle adjacent to second newly created
                }

        }
}


template <typename T>
void DT2D::swapDiagonal ( HalfEdge <T> *e_twin, HalfEdge <T> *e12, HalfEdge <T> *e13, HalfEdge <T> *e21, HalfEdge <T> *e22, HalfEdge <T> *e23 )
{
        // Swap diagonal in triangles (e11,e12,e13) and (e21,e22,e23)
        e_twin->setPoint ( e23->getPoint() );		// Change start point
        e_twin->setNextEdge ( e13 );			// Change next edge
        e_twin->setTwinEdge ( e21 );			// Change twin edge

        e13->setNextEdge ( e22 );			// Change next edge
        e22->setNextEdge ( e_twin );			// Change next edge

        // Swap edges in triangle adjacent to legalized
        e21->setPoint ( e13->getPoint() );		// Change start point
        e21->setNextEdge ( e23 );			// Change next edge
        e21->setTwinEdge ( e_twin );			// Change twin edge

        e23->setNextEdge ( e12 );			// Change next edge
        e12->setNextEdge ( e21 );			// Change next edge
}



template <typename T>
void DT2D::removeSimplexTriangles ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges_dt )
{
        // Remove simplex triangles
        const Node3DCartesian <T> *s1 = nl [0];  //Simplex s1
        const Node3DCartesian <T> *s2 = nl [1];  //Simplex s2
        const Node3DCartesian <T> *s3 = nl [2];  //Simplex s3

        // Check all edges
        unsigned int n = half_edges_dt.size();

        for ( unsigned int i = 0; i < n; i++ )
        {
                const Node3DCartesian <T> *p = half_edges_dt [i]->getPoint();

                // Remember simplex triangles
                if ( p == s1 || p == s2 || p == s3 )
                {
                        // Get edges belonging to simplex triangle
                        HalfEdge <T> *e2 = half_edges_dt [i]->getNextEdge();
                        HalfEdge <T> *e3 = e2->getNextEdge();

                        // Mark edges belonging to simplex triangle
                        half_edges_dt [i]->setEdgeAsSimplex ( true );
                        e2->setEdgeAsSimplex ( true );
                        e3->setEdgeAsSimplex ( true );

                        // Get twin edges of the simplex triangle
                        HalfEdge <T> *et1 = half_edges_dt [i]->getTwinEdge();
                        HalfEdge <T> *et2 = e2->getTwinEdge();
                        HalfEdge <T> *et3 = e3->getTwinEdge();

                        // Remove pointers from adjacent edges to simplex edges
                        if ( et1 != NULL )
                        {
                                et1->setTwinEdge ( NULL );
                        }

                        if ( et2 != NULL )
                        {
                                et2->setTwinEdge ( NULL );
                        }

                        if ( et3 != NULL )
                        {
                                et3->setTwinEdge ( NULL );
                        }
                }
        }

        // Remove simplex
        typename  TItemsList <HalfEdge<T> *>::Type ::iterator i_half_edges = std::remove_if ( half_edges_dt.begin(), half_edges_dt.end(), removeSimplexHalfEdge() );
        half_edges_dt.erase ( i_half_edges, half_edges_dt.end() );

        //Remove three simplex points
        delete s1;
        delete s2;
        delete s3;

        //Remove first three items of vector
        typename  TItemsList <Node3DCartesian <T> *>::Type ::iterator i_points = nl.begin();
        i_points = nl.erase ( i_points );
        i_points = nl.erase ( i_points );
        i_points = nl.erase ( i_points );

        // Set all edges as non simplex
        n = half_edges_dt.size();

        for ( unsigned int i = 0; i < n; i++ )
        {
                half_edges_dt [i] ->setEdgeAsSimplex ( false );
        }
}

#endif
