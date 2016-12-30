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


#ifndef Voronoi2D_H
#define Voronoi2D_H


#include "libalgo/source/structures/list/Container.h"


//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class HalfEdge;

template <typename T>
class VoronoiCell;

template <typename T>
class HalfEdgesList;


//Set type of generated VoronoiCells: all cells, only bounded cells, only bounded cells of the appropriate shape
typedef enum
{
        AllCells = 0,
        BoundedCells,
        AppropriateBoundedCells,
} TVoronoiCellsType;


//Set type of generated Voronoi diagram: error free approach or topologic approach by Okabe
typedef enum
{
        ErrorFreeApproach = 0, 
	TopologicApproach
} TVoronoiDiagramMethod;


//Construct 2D Voronoi diagrams using incremental construction
class Voronoi2D
{
        public:

                template <typename T>
                static void VD ( Container <Node3DCartesian <T> *> &points, Container <Node3DCartesian <T> *> &vor_points, Container <HalfEdge <T> *> &hl_dt, Container <HalfEdge <T> *> &hl_vor, Container <VoronoiCell <T> *> &vl, const TVoronoiCellsType cells_type = BoundedCells, const TVoronoiDiagramMethod vor_diagram_method = TopologicApproach, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static void mergeVoronoiCellAndAdjacentCells ( const VoronoiCell <T> *voronoi_cell, Face <T> ** output_face, Container <Node3DCartesian <T> *>  &intersections, Container <HalfEdge <T> *> &hl );

        private:

                template <typename T>
                static void createVoronoiCells ( Container <HalfEdge <T> *> &hl_dt, Container <HalfEdge <T> *> &hl_vor, Container <Node3DCartesian <T> *> &vor_points, Container <VoronoiCell <T> *> &vl,
                                                 const TVoronoiCellsType cells_type, const TVoronoiDiagramMethod vor_diagram_method );

                template <typename T>
                static void correctTopologyInVoronoiCells ( Container <HalfEdge <T> *> &hl_vor, Container <Node3DCartesian <T> *> &vor_points );

                template <typename T>
                static void removeUnboundedVoronoiCells ( Container <VoronoiCell <T> *> &vl );

};

#include "Voronoi2D.hpp"

#endif
