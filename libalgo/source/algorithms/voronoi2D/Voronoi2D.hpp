// Description: Compute 2D Voronoi digaram using incremental algorithm with full topology

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


#ifndef Voronoi2D_HPP
#define Voronoi2D_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <queue>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/line/HalfEdge.h"

#include "libalgo/source/structures/face/VoronoiCell.h"

#include "libalgo/source/algorithms/circle3points/Circle3Points.h"
#include "libalgo/source/algorithms/faceoverlay/FaceOverlay.h"
#include "libalgo/source/algorithms/triangulations2D/DT2D.h"


template <typename T>
void Voronoi2D::VD ( Container <Node3DCartesian <T> *> &points, Container <Node3DCartesian <T> *> &vor_points, Container <HalfEdge <T> *> &hl_dt, Container <HalfEdge <T> *> &hl_vor, Container <VoronoiCell <T> *> &vl,
                     const TVoronoiCellsType cells_type, const TVoronoiDiagramMethod vor_diagram_method, const bool print_message, const bool print_exception, std::ostream * output )
{
        //Create 2D Voronoi diagrams using duality between DT <-> VT
        try
        {
                //Print info
                if ( print_message )
                {
                        *output << "> Starting VD, please wait... ";
                        *output << "Cells type: ";
                        cells_type == AllCells ? *output << "all cells. \n" : *output << "unbounded cells. \n";
                        *output << "Voronoi diagram type: ";
                        vor_diagram_method == TopologicApproach ? *output << "topologic approach. \n" : *output << "error free approach. \n";
                }

                //Perform Delaunay triangulation
                DT2D::DT ( points, hl_dt, print_message, print_exception, output );
                //DXFExport::exportDTToDXF ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\out\\dt_error_merc.dxf", hl_dt );

                //At least 3 Delaunay edges (one triangle) were created
                if ( hl_dt.size() > 3 )
                {
                        //Construct Voronoi diagram cell by cell from each halfedge: edge -> twin -> previous
                        createVoronoiCells ( hl_dt, hl_vor, vor_points, vl, cells_type, vor_diagram_method );

                        //Correct Voronoi cells for topologic approach
                        correctTopologyInVoronoiCells ( hl_vor, vor_points );

                        //Remove unbounded Voronoi cells
                        if ( cells_type == BoundedCells || cells_type == AppropriateBoundedCells )
                        {
                                removeUnboundedVoronoiCells ( vl );
                        }

                        //Print info
                        if ( print_message )
                        {
                                *output << " Completed." << std::endl;
                        }
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                //Clear lists
                vor_points.clear();
                hl_dt.clear();
                hl_vor.clear();
                vl.clear();

                //Print exception
                if ( print_exception )
                {
                        error.printException();
                }

                throw error;
        }
}


template <typename T>
void Voronoi2D:: createVoronoiCells ( Container <HalfEdge <T> *> &hl_dt, Container <HalfEdge <T> *> &hl_vor, Container <Node3DCartesian <T> *> &vor_points, Container <VoronoiCell <T> *> &vl,
                                      const TVoronoiCellsType cells_type, const TVoronoiDiagramMethod vor_diagram_method )
{
        //Create Voronoi Diagram with full topolgy from Delaunay edges
        std::queue <HalfEdge <T> *> Q;

        //Add first edge to the list
        Q.push ( hl_dt [0] );

        //Process all edges int the queue
        while ( ! Q.empty() )
        {
                bool switch_direction = false, infinity_vertex = false;
                unsigned int vertices_count = 0;
                HalfEdge <T> *vor_edge = NULL, *vor_edge_previous = NULL, * vor_edge_first = NULL, *vor_edge_boundary = NULL;

                //Remember start edge
                HalfEdge <T> *e = Q.front();
                HalfEdge <T> *e_start = e;

                //Remove element on the top
                Q.pop();

                //Edge has not been processed
                if ( !e_start->isSimplexEdge() )
                {
                        //Get vertices of the Deluanay triangle
                        Node3DCartesian <T> *q1 = e->getPoint();
                        Node3DCartesian <T> *q2 = e->getNextEdge() -> getPoint();
                        Node3DCartesian <T> *q3 = e->getNextEdge()->getNextEdge()->getPoint();

                        //Remember start nodes
                        Node3DCartesian <T> *q1_start = q1;
                        Node3DCartesian <T> *q3_start = q3;

                        //Create new Voronoi Cell
                        VoronoiCell <T> *voronoi_cell_new = new VoronoiCell <T> ( q2, NULL,  1 );

                        //Walk from the edge e_start = e in direction: e -> twin -> previous, where
                        do
                        {
                                T x_center_circle = 0, y_center_circle = 0, r = 0;

                                vertices_count ++;

                                try
                                {
                                        //Compute center of the inscribed circle
                                        Circle3Points::getCentreAndDiameterCircle ( q1, q2, q3, x_center_circle, y_center_circle, r );

                                        //Test, if created Voronoi edge is not too long ( such a vertex forms Voronoi cell of inappropriate shape )
                                        if ( ( cells_type == AppropriateBoundedCells ) && ( r > 2 * max ( EuclDistance::getEuclDistance2D ( q1, q2 ),
                                                        max ( EuclDistance::getEuclDistance2D ( q2, q3 ), EuclDistance::getEuclDistance2D ( q3, q1 ) ) ) ) )
                                        {
                                                //Voronoi edge is too long: set vertex as infinity, cell will be considered as unbounded
                                                throw BadDataException ( "BadDataException: Inappropriate Voronoi vertex, too accute angle", " throw this vertex." ) ;
                                        }
                                }

                                //Points identical or nearby to colinear, center of the circle goes to infinity
                                catch ( Exception & error )
                                {
                                        //Set cell as unbounded
                                        infinity_vertex = true;
                                }

                                //Set this edge as processed (all edges on the same side of the triangle as edge e_start)
                                e->setEdgeAsSimplex ( true );

                                //Get twin edge
                                HalfEdge <T> *e_twin = e->getTwinEdge();
                                Node3DCartesian <T> *center = NULL;

                                //Add visited unprocessed twin edges to the queue
                                if ( ( e_twin != NULL ) && ( !e_twin->isSimplexEdge() ) )
                                {
                                        Q.push ( e_twin );
                                }

                                //Test if all Voronoi edges sharing this node = center have already been created: 3 Voronoi edges
                                const HalfEdge <T> *e21_dual = e->getNextEdge()->getDualEdge();
                                const HalfEdge <T> *e31_dual = e->getNextEdge()->getNextEdge()->getDualEdge();

                                //Keep topologic consistency of Voronoi cells: if computed Voronoi point too close to another
                                //Voronoi point, use this point
                                if ( ( vor_diagram_method == TopologicApproach ) && ( !switch_direction ) && ( vor_edge_previous != NULL ) &&
                                                ( EuclDistance::getEuclDistance2D ( vor_edge_previous->getPoint()->getX(), vor_edge_previous->getPoint()->getY(), x_center_circle, y_center_circle ) < MIN_POSITION_DIFF ) )
                                {
                                        center = vor_edge_previous->getPoint();
                                }
                                else if ( ( vor_diagram_method == TopologicApproach ) && ( vor_edge_first != NULL ) &&
                                                ( EuclDistance::getEuclDistance2D ( vor_edge_first->getPoint()->getX(), vor_edge_first->getPoint()->getY(), x_center_circle, y_center_circle ) < MIN_POSITION_DIFF ) )
                                {
                                        center = vor_edge_first->getPoint();
                                }

                                //Get start point of all Voronoi edges sharing vertex ( = center) if exist
                                else if ( e21_dual != NULL ) center = e21_dual->getPoint();
                                else if ( e31_dual != NULL ) center = e31_dual->getPoint();

                                //Otherwise create new Voronoi point
                                else
                                {
                                        center = new Node3DCartesian <T> ( x_center_circle,  y_center_circle );
                                        vor_points.push_back ( center );
                                }

                                //Create new Voronoi edge starting at this point
                                vor_edge = new HalfEdge <T> ( center, NULL, NULL );

                                //Set actual cell for HalfEdge
                                vor_edge->setFace ( voronoi_cell_new );

                                //Add new new Voronoi edge to the list
                                hl_vor.push_back ( vor_edge );

                                //Remember first created Voronoi edge (we link this edge with the last created Voronoi edge)
                                if ( vor_edge_first == NULL )
                                {
                                        vor_edge_first = vor_edge;
                                }

                                //Link actual voronoi edge and previous/first edge according to the direction of processing
                                if ( !switch_direction )
                                {
                                        //Link actual Voronoi edge and the previous Voronoi edge: add actual Voronoi edge after the previous edge
                                        if ( vor_edge_previous != NULL )
                                        {
                                                vor_edge_previous->setNextEdge ( vor_edge );
                                                vor_edge->setPreviousEdge ( vor_edge_previous );
                                        }
                                }

                                //Link actual Voronoi edge and the first Voronoi edge: we have already found first Voronoi edge
                                else
                                {
                                        vor_edge->setNextEdge ( vor_edge_first );
                                        vor_edge_first->setPreviousEdge ( vor_edge );

                                        //Set this edge as first Voronoi edge
                                        vor_edge_first = vor_edge;
                                }

                                //Set dual Delaunay edge e to this Voronoi edge and this Voronoi edge as dual to Delaunay edge e
                                vor_edge->setDualEdge ( e );
                                e->setDualEdge ( vor_edge );

                                //Set dual edge for twin edge to actual Voronoi edge e (in both directions of processing)
                                if ( e_twin != NULL )
                                {
                                        //Is there a twin Voronoi edge dual to e_twin ?
                                        if ( ( e_twin->getDualEdge() != NULL ) )
                                        {
                                                //Set twin Voronoi edge to actual Voronoi edge
                                                vor_edge->setTwinEdge ( e_twin->getDualEdge() );
                                                e_twin->getDualEdge()->setTwinEdge ( vor_edge );
                                        }
                                }

                                //There is no twin edge, cell is unbounded (test only in common direction of processing) : switch direction
                                else if ( !switch_direction )
                                {
                                        //Switch direction of the processing: add segments to the beginning
                                        switch_direction = true;

                                        //Assign start points: we start from the start triangle in the different direction
                                        q1 = q1_start;
                                        q3 = q3_start;

                                        //Assign start edge: we continue searching from the start edge in diffrerent direction
                                        //(add new Voronoi edge before the first created Voronoi edge)
                                        e = e_start;

                                        //Remember boundary edge: Voronoi edge assigned to Delaunay edge with e_twin = NULL;
                                        //We close Voronoi cell according to this edge
                                        vor_edge_boundary = vor_edge;
                                }

                                //Set new start edge and start point: in any case
                                //Common direction of the processing: we will add new Voronoi edge after the last created Voronoi edge
                                if ( !switch_direction )
                                {
                                        //Store actual Voronoi edge as previous
                                        vor_edge_previous = vor_edge;

                                        //Increment actual edge
                                        e = e_twin->getNextEdge()->getNextEdge();

                                        //Assign circle points
                                        q3 = q1;
                                        q1 = e->getPoint();
                                }

                                //Changed direction of the processing: we add new Voronoi edge before the first created Voronoi edge)
                                else
                                {
                                        //Increment actual edge
                                        e = e->getNextEdge()->getTwinEdge();

                                        //We found a boundary triangle, we processed both parts of an unbounded Voronoi cell
                                        if ( e == NULL )
                                        {
                                                break;
                                        }

                                        //Assign circle points
                                        q1 = q3;
                                        q3 = e ->getNextEdge()->getNextEdge()->getPoint();
                                }

                        }
                        while ( e != e_start );

                        //Close Voronoi edges in bounded cell
                        if ( !switch_direction )
                        {
                                vor_edge_first->setPreviousEdge ( vor_edge );
                                vor_edge->setNextEdge ( vor_edge_first );
                        }

                        //Pseudoclose Voronoi edges in unbounded cell
                        else
                        {
                                vor_edge_first->setPreviousEdge ( vor_edge_boundary );
                                vor_edge_boundary->setNextEdge ( vor_edge_first );
                        }

                        //Set properties of Voronoi Cell
                        voronoi_cell_new->setHalfEdge ( vor_edge );
                        voronoi_cell_new->setVerticesCount ( vertices_count );
                        voronoi_cell_new->setBounded ( !switch_direction && !infinity_vertex );

                        //Link Voronoi Cell and generator
                        q2->setFace ( voronoi_cell_new );

                        //Add Voronoi cell to the list
                        vl.push_back ( voronoi_cell_new );
                }
        }
}


template <typename T>
void Voronoi2D::correctTopologyInVoronoiCells ( Container <HalfEdge <T> *> &hl_vor, Container <Node3DCartesian <T> *> &vor_points )
{
        //Corect topology in Voronoi diagram
        const unsigned int n = hl_vor.size();

        for ( unsigned int i = 0; i < n; i++ )
        {
                //Get start point n1 of the edge (first Voronoi cell vc1 )
                Node3DCartesian <T> *n1 = hl_vor [i]->getPoint() ;

                //Test twin edge
                if ( hl_vor [i]->getTwinEdge() != NULL )
                {
                        //Get start point n2 of the edge (second Voronoi cell vc2 sharing vc1 : n1 = n2)
                        Node3DCartesian <T> *n2 = hl_vor [i]->getTwinEdge()->getNextEdge()->getPoint();

                        //Get start point n3 of the edge (third Voronoi cell vc3 sharing vc2 and vc3: n1 = n2 = n3)
                        Node3DCartesian <T> *n3 = NULL;

                        if ( hl_vor [i]->getPreviousEdge()->getTwinEdge() != NULL )
                        {
                                n3 = hl_vor [i]->getPreviousEdge()->getTwinEdge()->getPoint();
                        }

                        //Shared Voronoi edges in 3 adjacent cells do not have the same start point
                        if ( ( n1 != n2 )  || ( n3 != NULL ) && ( n1 != n3 || n2 != n3 ) )
                        {
                                //Correct start points of all edges
                                for ( unsigned int j = 0; j < n; j++ )
                                {
                                        if ( ( hl_vor [j]->getPoint() == n2 ) || ( ( n3 != NULL ) && ( hl_vor [j]->getPoint() == n3 ) ) )
                                        {
                                                hl_vor [j]->setPoint ( n1 );
                                        }
                                }
                        }
                }
        }
}


template <typename T>
void Voronoi2D::removeUnboundedVoronoiCells ( Container <VoronoiCell <T> *> &vl )
{
        //Remove unbounded cells
        typename TItemsList <VoronoiCell <T> *>::Type ::iterator i_cells = vl.begin();

        //Process all cells
        while ( i_cells != vl.end() )
        {
                //Unbounded cell, delete
                if ( ! ( *i_cells )->getBounded() )
                {
                        //Remove adjacency to adjacent cells
                        ( *i_cells )->removeAdjacency();

                        //Delete Face
                        delete ( *i_cells );

                        //Remove Face from list
                        i_cells = vl.erase ( i_cells );

                        //Jump to the next iteration
                        continue;
                }

                i_cells ++;
        }
}


template <typename T>
void Voronoi2D::mergeVoronoiCellAndAdjacentCells ( const VoronoiCell <T> *voronoi_cell, Face <T> ** output_face, Container <Node3DCartesian <T> *>  &intersections, Container <HalfEdge <T> *> &hl )
{
        //Merge Voronoi cell with all adjacent Voronoi cells and creates new face using overlay operations
        bool f1_simple = false, f2_simple = false;
        std::list <VoronoiCell<T> *> adjacent_cells;

        //Get half edge of the Voronoi cell
        const HalfEdge <T> * e_start = voronoi_cell->getHalfEdge();
        HalfEdge <T> * e = const_cast <HalfEdge <T> * > ( e_start );

        //Find all neighbouring Voronoi cells and add to the list
        do
        {
                //There is a twin edgeand edge is long enough
                if ( ( e->getTwinEdge() != NULL ) && ( EuclDistance::getEuclDistance2D ( e->getPoint(), e->getNextEdge()-> getPoint() ) > MIN_POSITION_DIFF ) )
                {
                        VoronoiCell <T> *voronoi_cell_adjacent = static_cast <VoronoiCell <T> *> ( e->getTwinEdge()->getFace() );

                        //Merge only bounded adjacent Voronoi cells
                        if ( ( voronoi_cell_adjacent != NULL ) && ( voronoi_cell_adjacent->getBounded() ) )
                        {
                                adjacent_cells.push_back ( voronoi_cell_adjacent );
                        }
                }

                e = e->getNextEdge();
        }
        while ( e != e_start );

        //Create results list and first Voronoi cell
        Container < Face <T> *> results;
        results.push_back ( ( dynamic_cast <Face <T> *> ( const_cast < VoronoiCell <T> *> ( voronoi_cell ) ) ) );

        //Process all cells
        bool switch_direction = false;

        while ( !adjacent_cells.empty() )
        {
                VoronoiCell <T> *voronoi_cell_next = NULL;

                //Get actual cell
                switch_direction ? voronoi_cell_next = adjacent_cells.back() : voronoi_cell_next = adjacent_cells.front();

                //Merge both faces
                try
                {
                        FaceOverlay::createOverlay ( results[results.size() - 1], voronoi_cell_next, Union, results, intersections, hl, f1_simple, f2_simple );
                }

                //One input face is not simple: last overlay operation processed cells in wrong direction
                catch ( BadDataException & error )
                {
                        //Bad result of the last merge operation: rare case, when only some adjacent Voronoi cells completely surround Voronoi cell
                        //Otherwise we delete unbounded adjacent cell and continue
                        if ( !f1_simple )
                        {
                                //Switch direction of the processing of Voronoi cells
                                switch_direction = true;

                                //Erase last result of the overlay operation: non valid result
                                if ( results.size() > 1 )
                                {
                                        results.pop_back();
                                        continue;
                                }

                                //Merged Voronoi cell is not simple: stop processing
                                break;
                        }
                }

                //Remove firs/last Voronoi cell according to the direction
                switch_direction ? adjacent_cells.pop_back() : adjacent_cells.pop_front();
        }

        //Assign result: no merging has been performed
        if ( results.size() == 1 )
        {
                *output_face = new Face <T> ( ( dynamic_cast <Face <T> *> ( const_cast < VoronoiCell <T> *> ( voronoi_cell ) ) ) );
        }

        //At least one merging has been performed
        else
        {
                *output_face = dynamic_cast <Face <T> *> ( results[results.size() - 1] ) ;

                //Remove last element (result) from the list not to be deleted by the destructor
                results.pop_back();
        }

        //*output_face = ( results.size() > 1 ? ( dynamic_cast <Face <T> *> ( results[results.size() - 1] ) )  : new Face <T> ( ( dynamic_cast <Face <T> *> ( const_cast < VoronoiCell <T> *> ( voronoi_cell ) ) ) ) );

        //Remove first element (input Voronoi cell) from the list not to be deleted by the destructor
        results.pop_front();
}

#endif
