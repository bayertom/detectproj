// Description: Constrained 2D Delaunay triangulation (Sloan, 1992)

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


#ifndef CDT2D_HPP
#define CDT2D_HPP

#include <iomanip>
#include <fstream>
#include <cmath>

#include "libalgo/source/structures/line/Edge.h"

#include "libalgo/source/structures/point/Node3DCartesian.h"

#include "libalgo/source/algorithms/lawsonorientedwalk/LawsonOrientedWalk.h"
#include "libalgo/source/algorithms/linelineposition/LineLinePosition.h"
#include "libalgo/source/algorithms/swappingcriteria/SwappingCriteria.h"
#include "libalgo/source/algorithms/convexquadrilateral/ConvexQuadrilateral.h"

#include "DT2D.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
void CDT2D::CDT ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> & half_edges, Container <Edge < Node3DCartesian <T> > *> &constrained_edges,  const bool print_message, const bool print_exception, std::ostream * output )
{
        // Constrained Delaunay triangulation

        //Create DT
        DT2D::DT ( nl, half_edges, print_message, print_exception );

        //Initialize edge
        HalfEdge <T> *e_heuristic = ( *half_edges ) [0];

        //Print info
        if ( print_message )
        {
                *output << "> Starting CDT, please wait..." << std::endl;
                *output << "Total: " << half_edges->size() << " constrained edges." << std::endl;
        }

        //Throw exception
        try
        {
                //Process all constrained edges
                while ( !constrained_edges->empty() )
                {
                        // Get constrained edge
                        Edge <Node3DCartesian <T> > *e_constrained = constrained_edges->front();

                        // Run constrained Delaunay triangulation for this edge
                        insertConstrainedEdge ( e_constrained, e_heuristic, constrained_edges );

                        // Pop this edge
                        constrained_edges->pop_front();
                }

                //Print results
                if ( print_message )
                {
                        *output << "Completed..." << std::endl;
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                //Clear constrained edges
                constrained_edges->clear();

                //Print exception
                if ( print_exception )
                {
                        error.printException ( output );
                }

                //Throw exception
                throw;
        }
}


template <typename T>
void CDT2D::insertConstrainedEdge ( Edge <Node3DCartesian <T> > *e_constrained, HalfEdge <T> *e_heuristic, Container <Edge <Node3DCartesian <T> > *> &constrained_edges )
{
        //Insert one constrained edge
        typename THalfEdgesCDTList <T>::Type half_edges_created;
        typename THalfEdgesCDTList <T>::Type half_edges_removed;

        // Constrained Delaunay triangulation by Sloan(1992)
        findIntersectedEdges ( e_constrained, e_heuristic, half_edges_removed, constrained_edges );

        // Remove edges intersecting constrained edge
        removeIntersectedEdges ( e_constrained, half_edges_removed, half_edges_created );

        // Restore Delaunay triangulation over edge
        restoreDT ( e_constrained, half_edges_removed, half_edges_created );
}


template <typename T>
void CDT2D::findIntersectedEdges ( Edge <Node3DCartesian <T> > *e_constrained, HalfEdge <T> *he, typename THalfEdgesCDTList <T>::Type & half_edges_removed, Container <Edge <Node3DCartesian <T> > *> &constrained_edges )
{

        /* Find edges intersection constrained edge ;
        Used heuristic: finding starts from the last triangle defined by edge e */

        bool left = true;	// Searching direction
        short status = -1;

        // Is a constrained edge already in triangulation?
        he = LawsonOrientedWalk::findTriangleWalk ( &e_constrained->getP1(), &status, he, 0 );

        // If found edge has an opposite orientation (n1==n4)
        if ( ( * he->getPoint() ) != ( e_constrained->getP1() ) )
        {
                // Switch orientation
                he = he->getTwinEdge();
        }

        // Get start and end point of the edge
        const Node3DCartesian <T> *n1 = &e_constrained->getP1();
        Node3DCartesian <T> *n2 = &e_constrained->getP2();

        try
        {
                // Found point has same start point as the constrainend edge
                if ( he != NULL )
                {
                        // Loop over all edges intersecting constrained edge
                        for ( ;; )
                        {
                                // Next edge
                                he = he->getNextEdge();

                                // Intersect edge e constrained edge?
                                const HalfEdge <T> *e_next = he->getNextEdge();

                                // Get nodes of the next analyzed edge
                                const Node3DCartesian <T> *n3 = he->getPoint();
                                const Node3DCartesian <T> *n4 = e_next->getPoint();

                                const short inters2 = LineLinePosition::get2LineSegmentsPosition ( n1, n2, n3, n4, 0 );

                                // Constrained edge intersects some edge
                                if ( inters2 == 1 )
                                {
                                        // Remove intersected edge from triangulation
                                        half_edges_removed.push_back ( he );

                                        // Set incident edge of the adjacent triangle
                                        he = he->getTwinEdge();
                                }

                                // Constrained edge and Delaunay edge intersect in 1 end point: brake constrained edge in 2 parts
                                else if ( inters2 == 2 )
                                {
                                        // Get start point of intersected edge
                                        const HalfEdge <T> *e_next = he->getNextEdge();
                                        Node3DCartesian <T> *n = e_next->getPoint();

                                        // Correct end point of the first constrained edge
                                        e_constrained->setP2 ( n );

                                        // Create new constrained edge
                                        Edge <Node3DCartesian <T> > *ce2 = new Edge <Node3DCartesian <T> > ( n, n2 );

                                        // Add next part of the broken constrained edge into list
                                        constrained_edges->push_back ( ce2 );

                                        // Stop searching
                                        break;
                                }


                                // Constrained edge and Delaunay edge intersect in 2 endpoints: brake constrained edge in 2 edges
                                else if ( inters2 == 3 )
                                {
                                        // Get start point end end point
                                        const Node3DCartesian <T> *ns = he->getPoint();
                                        const Node3DCartesian <T> *ne = he->getNextEdge()->getPoint();	// We found end point of the constrained edge, this edge was not collinear

                                        // Is it endpoint of the edge
                                        if ( ns == n2 || ne == n2 )
                                        {
                                                break;	// We found endpoint
                                        }

                                        // We are again in start point of the constrained edge n1
                                        else
                                        {
                                                // Switch to the left triangle
                                                if ( left )
                                                {
                                                        // Good way, switch the triangle
                                                        if ( he->getTwinEdge() != NULL )
                                                        {
                                                                // Switch triangle to left incident triangle
                                                                he = he->getTwinEdge();
                                                        }

                                                        // Wrong way, change direction to the right
                                                        else
                                                        {
                                                                left = false;
                                                        }
                                                }

                                                // Switch to the right triangle
                                                if ( !left )
                                                {
                                                        // Go to next edge
                                                        he = he->getNextEdge();

                                                        // Switch to the right incident triangle
                                                        he = he->getTwinEdge();

                                                        // Go to the next edge with start point different from n1
                                                        he = he->getNextEdge();
                                                }
                                        }
                                }
                        }
                }
        }

        //Throw an exception
        catch ( MathZeroDevisionException <T> & e )
        {
                //Clear half edges
                half_edges_removed.clear();

                //Throw exception
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: ", "Constrained Delaunay triangulation: Can not find intersected edges." );
        }
}


template <typename T>
void CDT2D::removeIntersectedEdges ( const Edge <Node3DCartesian <T> > *e_constrained, typename THalfEdgesCDTList <T>::Type & half_edges_removed, typename THalfEdgesCDTList <T>::Type & half_edges_created )
{
        //Remove edges intersecting constrained edge

        //Get start and end point of the edge
        const Node3DCartesian <T> *n1 = &e_constrained->getP1();
        const Node3DCartesian <T> *n2 = &e_constrained->getP2();

        try
        {
                while ( !half_edges_removed.empty() )
                {
                        // Take the first edge
                        HalfEdge <T> *e_swap = half_edges_removed.front();

                        // Pop this edge
                        half_edges_removed.pop_front();

                        // Test of convexity for quadrilateral
                        HalfEdge <T> *e12 = e_swap->getNextEdge();
                        HalfEdge <T> *e13 = e12->getNextEdge();
                        HalfEdge <T> *e21 = e_swap->getTwinEdge();
                        HalfEdge <T> *e22 = e21->getNextEdge();
                        HalfEdge <T> *e23 = e22->getNextEdge();

                        // Get nodes, counterclockwise set of nodes
                        Node3DCartesian <T> *n11 = e_swap->getPoint();
                        Node3DCartesian <T> *n23 = e23->getPoint();
                        Node3DCartesian <T> *n12 = e12->getPoint();
                        Node3DCartesian <T> *n13 = e13->getPoint();

                        //Is quadrilateral strictly convex?
                        if ( ConvexQuadrilateral::isStrictlyConvex ( static_cast< Point3DCartesian <T> * > ( n11 ), static_cast < Point3DCartesian <T> * > ( n23 ), static_cast < Point3DCartesian <T> * > ( n12 ), static_cast < Point3DCartesian <T> * > ( n13 ) ) == 1, true )
                        {
                                // Swap diagonal
                                DT2D::swapDiagonal ( e_swap, e12, e13, e21, e22, e23 );

                                // Intersect edge e constrained edge?
                                const HalfEdge <T> *e_next = e_swap->getNextEdge();

                                // Get nodes of the next analyzed edge
                                const Node3DCartesian <T> *n3 = e_swap->getPoint();
                                const Node3DCartesian <T> *n4 = e_next->getPoint();

                                // Test swapped diagonal
                                const short inters = LineLinePosition::get2LineSegmentsPosition ( n1, n2, n3, n4, 0 );

                                if ( inters == 1 )
                                {
                                        // Add swapped edge again to list of removed edges
                                        half_edges_removed.push_back ( e_swap );
                                }

                                // Add swapped edge to list of newly created edges
                                else
                                {
                                        half_edges_created.push_back ( e_swap );
                                }
                        }

                        // Non convex
                        else
                        {
                                // Add this edge again to list of removed edges
                                half_edges_removed.push_back ( e_swap );
                        }

                }
        }

        //Throw exception
        catch ( MathZeroDevisionException <T> & e )
        {
                //Clear half edges
                half_edges_removed.clear();
                half_edges_created.clear();

                //Throw exception
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: ", "Constrained Delaunay triangulation: Can not remove intersected edges." );
        }
}


template <typename T>
void CDT2D::restoreDT ( const Edge <Node3DCartesian <T> > *e_constrained, typename THalfEdgesCDTList <T>::Type & half_edges_removed, typename THalfEdgesCDTList <T>::Type & half_edges_created )
{
        // Restore Delaunay triangulation
        bool swap = true;	// Swap identificator
        typename THalfEdgesCDTList <T>::Type::iterator i_half_edges_created;

        // Get start and end point of the edge
        const Node3DCartesian <T> *n1 = &e_constrained->getP1();
        const Node3DCartesian <T> *n2 = &e_constrained->getP2();

        // Loop until no further swap is available
        i_half_edges_created = half_edges_created.begin();

        while ( swap )
        {
                swap = false;	// Set no swap will be done

                // Loop all newly created edges
                while ( i_half_edges_created != half_edges_created.end() )
                {
                        // Take edge e
                        HalfEdge <T> *e_swap = *i_half_edges_created;
                        HalfEdge <T> *e_swap_next = e_swap->getNextEdge();

                        // Is edge e = e_constrained?
                        const Node3DCartesian <T> *n3 = e_swap->getPoint();
                        const Node3DCartesian <T> *n4 = e_swap_next->getPoint();

                        // Is edge constrained edge?
                        if ( ! ( ( ( n1 == n3 ) && ( n2 == n4 ) ) || ( ( n1 == n4 ) && ( n2 == n3 ) ) ) )
                        {
                                // Legalize triangle
                                if ( legalizeConstrainedTriangle ( e_swap ) )
                                {
                                        swap = true;
                                }
                        }

                        //Increment iterator
                        i_half_edges_created++;
                }

        }

        // Empty list of newly created half-edges
        half_edges_created.clear();
        half_edges_removed.clear();
}


template <typename T>
bool CDT2D::legalizeConstrainedTriangle ( HalfEdge <T> *e_twin )
{

        // Legalize triangle in CDT (without recursive legalization) */
        HalfEdge <T> *e21 = e_twin->getTwinEdge();		/* First edge in triangle adjacent to legalized */
        HalfEdge <T> *e22 = e21->getNextEdge();		/* Second edge in triangle adjacent to legalized */
        HalfEdge <T> *e23 = e22->getNextEdge();		/* Third edge in triangle adjacent to legalized */

        // Adjacent triangle
        HalfEdge <T> *e12 = e_twin->getNextEdge();
        HalfEdge <T> *e13 = e12->getNextEdge();

        // Illegal triangle?
        if ( SwappingCriteria::getClineRenka ( e13->getPoint(), e21->getPoint(), e22->getPoint(), e23->getPoint() ) )
        {
                //Swap diagonal in (e11,e12,e13) and (e21,e22,e23)
                DT2D::swapDiagonal ( e_twin, e12, e13, e21, e22, e23 );

                // HalfEdge <T> has been swapped
                return true;
        }

        // Legal triangle
        else
        {
                //HalfEdge <T> has not been not swapped
                return false;
        }
}

#endif
