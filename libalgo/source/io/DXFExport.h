// Description: Export lines, points, polygons to 2D/3D DXF file

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


#ifndef DXFExport_H
#define DXFExport_H

#include <stdio.h>

#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class VoronoiCell;

template <typename T>
struct TMeridiansListF;

template <typename T>
struct TParallelsListF;


//Export to DXF file
class DXFExport
{
        public:
                template <typename T>
                static void exportDTToDXF ( const char * file_name, const Container <HalfEdge <T> *> &hl );

                template <typename T>
                static void exportVDToDXF ( const char * file_name, const Container <VoronoiCell <T> *> &vl, const bool process_all_cells );

                template <typename T, TDestructable destructable>
                static void exportFacesToDXF ( const char * file_name, const Container <Face <T> *, destructable> &faces );

                template <typename Point, TDestructable destructable>
		static void exportPointsToDXF(const char * file_name, const Container <Point *, destructable> &points, const typename Point::Type font_height );

                template <typename Point>
                static void exportGraticuleToDXF ( const char * file_name, const typename TMeridiansListF <typename Point::Type> ::Type & meridians, const typename TParallelsListF <typename Point::Type> ::Type & parallels,
                                                   const Container <Point *> &points, const typename Point::Type font_height, const typename Point::Type lat_step, const typename Point::Type lon_step );

        private:

                static void createHeaderSection ( FILE * file );
                static void createTableSection ( FILE * file );
                static void endTableSection ( FILE * file );
                static void createLayerSection ( FILE * file, const char * layer_name, const unsigned int color );
                static void createEntitySection ( FILE * file );

                template <typename T>
                static void createLine ( FILE * file, const char * layer_name, const T x1, const T y1, const T z1, const T x2, const T y2, const T z2 );

                template <typename T>
                static void createPoint ( FILE * file, const char * layer_name, const T x, const T y, const T z );

                template <typename T>
                static void createText ( FILE * file, const char * layer_name, const char * text, const T x, const T y, const T z, const T rotation, const T height );

                template <typename T>
                static void processHalfEdges ( FILE * file, const Container <HalfEdge <T> *> &hl, const char * layer_name );

                template <typename T>
                static void processFace ( FILE * file, const Face <T> * face, const char * layer_name_edges );

                template <typename T>
                static void processVoronoiCells ( FILE * file, const Container <VoronoiCell <T> *> &vl, const char * layer_name_edges, const char * layer_name_generators , const bool process_all_cells );

                template <typename GraticulePart, typename Point>
                static void processGraticuleElements ( FILE * file, GraticulePart & part, const Container <Point*> & points, const char * layer_name_edges, const char * layer_name_generators, const typename Point::Type font_height, const typename Point::Type step );

                static void endHeaderSection ( FILE * file );
};

#include "DXFExport.hpp"

#endif
