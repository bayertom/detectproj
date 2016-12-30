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



#ifndef DXFExport_HPP
#define DXFExport_HPP

#include <string.h>
#include <stdlib.h>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/VoronoiCell.h"

#include "libalgo/source/algorithms/bearing/Bearing.h"

#include "libalgo/source/exceptions/FileWriteException.h"


template <typename T>
void DXFExport::exportDTToDXF ( const char * file_name, const Container <HalfEdge <T> *> &hl )
{
        //Export Delaunay triangulation to DXF
        FILE * file = NULL;

        try
        {
                if ( ( file = fopen ( file_name, "w" ) ) != NULL )
                {
                        //Create header section
                        createHeaderSection ( file );

                        //Create table section
                        createTableSection ( file );

                        //Create layer
                        createLayerSection ( file, "DT", 3 );

                        //End table header
                        endTableSection ( file );

                        //Create entity section
                        createEntitySection ( file );

                        //Process half edges
                        processHalfEdges ( file, hl, "DT" );

                        //End header section
                        endHeaderSection ( file );

                        //Close file
                        fclose ( file );
                }

                //Throw exception
                else
                {
                        //Can not open file
                        throw std::ios_base::failure ( "Error: can not open the file. " );
                }
        }

        //Any error has appeared
        catch ( std::ios_base::failure & )
        {
                //Close file
                fclose ( file );

                //Throw exception
                throw FileWriteException ( "FileWriteException: can not write the file: ", file_name );
        }
}


template <typename T>
void DXFExport::exportVDToDXF ( const char * file_name, const Container <VoronoiCell <T> *> &vl, const bool process_all_cells )
{
        //Export Voronoi diagrams to DXF
        FILE * file = NULL;

        try
        {
                if ( ( file = fopen ( file_name, "w" ) ) != NULL )
                {
                        //Create header section
                        createHeaderSection ( file );

                        //Create table section
                        createTableSection ( file );

                        //Create layer for Voronoi points
                        createLayerSection ( file, "VD", 3 );

                        //Create layer for Voronoi points
                        createLayerSection ( file, "Generators", 3 );

                        //Create layer for Voronoi generators
                        //createLayerSection(file, "Voronoi_generators", 3);

                        //End table header
                        endTableSection ( file );

                        //Create entity section
                        createEntitySection ( file );

                        //Process Voronoi cells
                        processVoronoiCells ( file, vl, "VD", "Generators" , process_all_cells );

                        //End header section
                        endHeaderSection ( file );

                        //Close file
                        fclose ( file );
                }

                //Throw exception
                else
                {
                        //Can not open file
                        throw std::ios_base::failure ( "Error: can not open the file. " );
                }
        }

        //Any error has appeared
        catch ( std::ios_base::failure & )
        {
                //Close file
                fclose ( file );

                //Throw exception
                throw FileWriteException ( "FileWriteException: can not write the file:", file_name );
        }
}


template <typename T, TDestructable destructable>
void DXFExport::exportFacesToDXF ( const char * file_name, const Container <Face <T> *, destructable> &faces )
{
        //Export list of faces to DXF
        FILE * file = NULL;

        try
        {
                if ( ( file = fopen ( file_name, "w" ) ) != NULL )
                {
                        //Create header section
                        createHeaderSection ( file );

                        //Create table section
                        createTableSection ( file );

                        //Create layer for Voronoi points
                        createLayerSection ( file, "Faces", 3 );

                        //End table header
                        endTableSection ( file );

                        //Create entity section
                        createEntitySection ( file );

                        //Process all faces
                        for ( unsigned int i = 0; i < faces->size(); i++ )
                        {
                                processFace ( file, faces [i], "Faces" );
                        }


                        //End header section
                        endHeaderSection ( file );

                        //Close file
                        fclose ( file );
                }

                //Throw exception
                else
                {
                        //Can not open file
                        throw std::ios_base::failure ( "Error: can not open the file. " );
                }
        }

        //Any error has appeared
        catch ( std::ios_base::failure & )
        {
                //Close file
                fclose ( file );

                //Throw exception
                throw FileWriteException ( "FileWriteException: can not write the file: ", file_name );
        }
}


template <typename Point, TDestructable destructable>
void DXFExport::exportPointsToDXF(const char * file_name, const Container <Point *, destructable> &points, const typename Point::Type font_height)
{
        //Export points to DXF
        FILE * file = NULL;

        try
        {
                if ( ( file = fopen ( file_name, "w" ) ) != NULL )
                {
                        //Create header section
                        createHeaderSection ( file );

                        //Create table section
                        createTableSection ( file );

                        //Create layer for points
                        createLayerSection ( file, "Points", 3 );

                        //Create layer for point labels
                        createLayerSection ( file, "Point_labels", 3 );

                        //End table header
                        endTableSection ( file );

                        //Create entity section
                        createEntitySection ( file );

                        //Process all points
                        for ( unsigned int i = 0; i < points.size(); i++ )
                        {
                                //Create point
                                createPoint ( file, "Points", points [i]->getX(), points [i]->getY(), points [i]->getZ() );

                                //Create text label for each point
                                const char *point_label =  points [i]->getPointLabel().c_str();
                                if ( point_label != NULL )
                                {
					createText(file, "Point_labels", point_label, points[i]->getX() + 0.5 * font_height, points[i]->getY() - 0.5 * font_height, points[i]->getZ(), (typename Point::Type) 0, (typename Point::Type) font_height);
                                }

                                else
                                {
                                        char point_id_text [255];
                                        sprintf ( point_id_text, "%d", point_label );
					createText(file, "Point_labels", point_id_text, points[i]->getX() + 0.5 * font_height, points[i]->getY() - 0.5 * font_height, points[i]->getZ(), (typename Point::Type) 0, (typename Point::Type) font_height);
                                }
                        }

                        //End header section
                        endHeaderSection ( file );

                        //Close file
                        fclose ( file );
                }

                //Throw exception
                else
                {
                        //Can not open file
                        throw std::ios_base::failure ( "Error: can not open the file. " );
                }
        }

        //Any error has appeared
        catch ( std::ios_base::failure & )
        {
                //Close file
                fclose ( file );

                //Throw exception
                throw FileWriteException ( "FileWriteException: can not write the file", file_name );
        }
}


template <typename Point>
void DXFExport::exportGraticuleToDXF ( const char * file_name, const typename TMeridiansListF <typename Point::Type> ::Type & meridians,
                                       const typename TParallelsListF <typename Point::Type> ::Type & parallels, const Container <Point *> &points,
                                       const typename Point::Type font_height, const typename Point::Type lat_step, const typename Point::Type lon_step )
{
        //Export generated meridians and parallels to DXF
        FILE * file = NULL;

        try
        {
                if ( ( file = fopen ( file_name, "w" ) ) != NULL )
                {
                        //Create header section
                        createHeaderSection ( file );

                        //Create table section
                        createTableSection ( file );

                        //Create layer for meridians
                        createLayerSection ( file, "Meridians", 3 );

                        //Create layer for parallels
                        createLayerSection ( file, "Parallels", 3 );

                        //Create layer for meridian labels
                        createLayerSection ( file, "Meridian_labels", 3 );

                        //Create layer for parallel labels
                        createLayerSection ( file, "Parallel_labels", 3 );

                        //End table header
                        endTableSection ( file );

                        //Create entity section
                        createEntitySection ( file );

                        //Process all meridians
                        typename TMeridiansListF <typename Point::Type> ::Type::const_iterator i_meridians = meridians.begin();

                        for ( i_meridians = meridians.begin(); i_meridians != meridians.end(); ++i_meridians )
                        {
                                processGraticuleElements ( file, *i_meridians, points, "Meridians", "Meridian_labels", font_height, lon_step );
                        }

                        //Process all parallels
                        typename TParallelsListF <typename Point::Type> ::Type::const_iterator i_parallels = parallels.begin();

                        for ( i_parallels = parallels.begin(); i_parallels != parallels.end(); ++i_parallels )
                        {
                                processGraticuleElements ( file, *i_parallels, points, "Parallels", "Parallel_labels", font_height, lat_step );
                        }

                        //End header section
                        endHeaderSection ( file );

                        //Close file
                        fclose ( file );
                }

                //Throw exception
                else
                {
                        //Can not open file
                        throw std::ios_base::failure ( "Error: can not open the file. " );
                }
        }

        //Any error has appeared
        catch ( std::ios_base::failure & )
        {
                //Close file
                fclose ( file );

                //Throw exception
                throw FileWriteException ( "FileWriteException: can not write the file: ", file_name );
        }
}


inline void DXFExport::createHeaderSection ( FILE * file )
{
        //Create header section
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "SECTION\n", 8, 1, file );
        fwrite ( "2\n", 2, 1, file );
        fwrite ( "HEADER\n", 7, 1, file );
        fwrite ( "9\n", 2, 1, file );
        fwrite ( "$ACADVER\n", 9, 1, file );
        fwrite ( "1\n", 2, 1, file );
        fwrite ( "AC1006\n", 7, 1, file );
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "ENDSEC\n", 7, 1, file );
}


inline void DXFExport::endHeaderSection ( FILE * file )
{
        //Create end of the header section
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "ENDSEC\n", 7, 1, file );
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "EOF\n", 4, 1, file );
}


inline void DXFExport::createTableSection ( FILE * file )
{
        //Create table section
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "SECTION\n", 8, 1, file );
        fwrite ( "2\n", 2, 1, file );
        fwrite ( "TABLES\n", 7, 1, file );
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "TABLE\n", 6, 1, file );
        fwrite ( "2\n", 2, 1, file );
        fwrite ( "LAYER\n", 6, 1, file );
        fwrite ( "70\n", 3, 1, file );
        fwrite ( "0\n", 2, 1, file );
}


inline void DXFExport::endTableSection ( FILE * file )
{
        //Write end of the table section
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "ENDTAB\n", 7, 1, file );
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "ENDSEC\n", 7, 1, file );
}


inline void DXFExport::createLayerSection ( FILE * file, const char * layer_name, const unsigned int color )
{
        //Add section for one layer
        char layer_name_line[MAX_TEXT_LENGTH];
        strcpy ( layer_name_line, layer_name );
        strcat ( layer_name_line, "\n" );
        char color_text[32];
        sprintf ( color_text, "%i", color );
        //itoa ( color, color_text, 10 );
        strcat ( color_text, "\n" );

        fwrite ( "0\n", 2, 1, file );
        fwrite ( "LAYER\n", 6, 1, file );
        fwrite ( "2\n", 2, 1, file );
        fwrite ( layer_name_line, strlen ( layer_name_line ), 1, file );
        fwrite ( "70\n", 3, 1, file );
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "62\n", 3, 1, file );
        fwrite ( color_text, strlen ( color_text ), 1, file );
        fwrite ( "6\n", 2, 1, file );
        fwrite ( "CONTINUOUS\n", 11, 1, file );
}


inline void DXFExport::createEntitySection ( FILE * file )
{
        //Create section for entities
        fwrite ( "0\n", 2, 1, file );
        fwrite ( "SECTION\n", 8, 1, file );
        fwrite ( "2\n", 2, 1, file );
        fwrite ( "ENTITIES\n", 9, 1, file );
}


template <typename T>
void DXFExport::createLine ( FILE * file, const char * layer_name, const T x1, const T y1, const T z1, const T x2, const T y2, const T z2 )
{
        //Write line to DXF file
        char layer_name_line[MAX_TEXT_LENGTH];
        strcpy ( layer_name_line, layer_name );
        strcat ( layer_name_line, "\n" );

        const char * entity_id = "0\n";
        const char * entity_name = "LINE\n";
        const char * level_id = "8\n";
        const char * xi_id = "10\n";
        const char * yi_id = "20\n";
        const char * zi_id = "30\n";
        const char * xii_id = "11\n";
        const char * yii_id = "21\n";
        const char * zii_id = "31\n";

        //Convert number to char
        char x1_text[4096];
        char y1_text[MAX_TEXT_LENGTH];
        char z1_text[MAX_TEXT_LENGTH];
        char x2_text[MAX_TEXT_LENGTH];
        char y2_text[MAX_TEXT_LENGTH];
        char z2_text[MAX_TEXT_LENGTH];

        //Float to char
        sprintf ( x1_text, "%f", x1 );
        sprintf ( y1_text, "%f", y1 );
        sprintf ( z1_text, "%f", z1 );
        sprintf ( x2_text, "%f", x2 );
        sprintf ( y2_text, "%f", y2 );
        sprintf ( z2_text, "%f", z2 );

        //Add end line char
        strcat ( x1_text, "\n" );
        strcat ( y1_text, "\n" );
        strcat ( z1_text, "\n" );
        strcat ( x2_text, "\n" );
        strcat ( y2_text, "\n" );
        strcat ( z2_text, "\n" );

        /* Add to file */
        fwrite ( entity_id, 2, 1, file );
        fwrite ( entity_name, 5, 1,  file );
        fwrite ( level_id, 2, 1, file );
        fwrite ( layer_name_line, strlen ( layer_name_line ), 1, file );
        fwrite ( xi_id, 3, 1, file );
        fwrite ( x1_text, strlen ( x1_text ), 1, file );
        fwrite ( yi_id, 3, 1, file );
        fwrite ( y1_text, strlen ( y1_text ), 1, file );
        fwrite ( zi_id, 3, 1, file );
        fwrite ( z1_text, strlen ( z1_text ), 1, file );
        fwrite ( xii_id, 3, 1, file );
        fwrite ( x2_text, strlen ( x2_text ), 1, file );
        fwrite ( yii_id, 3, 1, file );
        fwrite ( y2_text, strlen ( y2_text ), 1, file );
        fwrite ( zii_id, 3, 1, file );
        fwrite ( z2_text, strlen ( z2_text ), 1, file );
}


template <typename T>
void DXFExport::createPoint ( FILE * file, const char * layer_name, const T x, const T y, const T z )
{
        //Write point to DXF file
        char layer_name_line[MAX_TEXT_LENGTH];
        strcpy ( layer_name_line, layer_name );
        strcat ( layer_name_line, "\n" );

        const char * entity_id = "0\n";
        const char * entity_name = "POINT\n";
        const char * level_id = "8\n";
        const char * color_id = "62\n";
        const char * entity_color = "5\n";
        const char * xi_id = "10\n";
        const char * yi_id = "20\n";
        const char * zi_id = "30\n";

        //Convert number to char
        char x_text[MAX_TEXT_LENGTH];
        char y_text[MAX_TEXT_LENGTH];
        char z_text[MAX_TEXT_LENGTH];

        //Float to char
        sprintf ( x_text, "%f", x );
        sprintf ( y_text, "%f", y );
        sprintf ( z_text, "%f", z );

        //Add end line char
        strcat ( x_text, "\n" );
        strcat ( y_text, "\n" );
        strcat ( z_text, "\n" );

        /* Add to file */
        fwrite ( entity_id, 2, 1, file );
        fwrite ( entity_name, 6, 1,  file );
        fwrite ( level_id, 2, 1, file );
        fwrite ( layer_name_line, strlen ( layer_name_line ), 1, file );
        fwrite ( color_id, 3, 1, file );
        fwrite ( entity_color, 2, 1, file );
        fwrite ( xi_id, 3, 1, file );
        fwrite ( x_text, strlen ( x_text ), 1, file );
        fwrite ( yi_id, 3, 1, file );
        fwrite ( y_text, strlen ( y_text ), 1, file );
        fwrite ( zi_id, 3, 1, file );
        fwrite ( z_text, strlen ( z_text ), 1, file );
}


template <typename T>
void DXFExport::createText ( FILE * file, const char * layer_name, const char * text, const T x, const T y, const T z, const T rotation, const T height )
{
        //Create text
        char layer_name_line[MAX_TEXT_LENGTH];
        strcpy ( layer_name_line, layer_name );
        strcat ( layer_name_line, "\n" );

        const char * entity_id = "0\n";
        const char * entity_name = "TEXT\n";
        const char * style_id = "7\n";
        const char * text_style = "PNTNUM\n";
        const char * rotation_id = "50\n";
        const char * level_id = "8\n";
        const char * color_id = "62\n";
        const char * entity_color = "10\n";
        const char * xi_id = "10\n";
        const char * yi_id = "20\n";
        const char * zi_id = "30\n";
        const char * height_id = "40\n";
        const char * text_id = "1\n";

        //Convert number to char
        char x_text[MAX_TEXT_LENGTH];
        char y_text[MAX_TEXT_LENGTH];
        char z_text[MAX_TEXT_LENGTH];
        char text_rotation[MAX_TEXT_LENGTH];
        char font_height[MAX_TEXT_LENGTH];
        char text_text[MAX_TEXT_LENGTH];

        //Float to char
        sprintf ( x_text, "%f", x );
        sprintf ( y_text, "%f", y );
        sprintf ( z_text, "%d", ( T ) z );
        sprintf ( text_rotation, "%f", rotation );
        sprintf ( font_height, "%f", height );

        //Add end line char
        strcat ( x_text, "\n" );
        strcat ( y_text, "\n" );
        strcat ( z_text, "\n" );
        strcat ( text_rotation, "\n" );
        strcat ( font_height, "\n" );
        strcpy ( text_text, text );
        strcat ( text_text, "\n" );

        /* Add to file */
        fwrite ( entity_id, 2, 1, file );
        fwrite ( entity_name, 5, 1,  file );
        fwrite ( style_id, 2, 1,  file );
        fwrite ( text_style, 7, 1,  file );
        fwrite ( rotation_id, 3, 1,  file );
        fwrite ( text_rotation, strlen ( text_rotation ), 1,  file );
        fwrite ( level_id, 2, 1, file );
        fwrite ( layer_name_line, strlen ( layer_name_line ), 1, file );
        fwrite ( color_id, 3, 1, file );
        fwrite ( entity_color, 3, 1, file );
        fwrite ( xi_id, 3, 1, file );
        fwrite ( x_text, strlen ( x_text ), 1, file );
        fwrite ( yi_id, 3, 1, file );
        fwrite ( y_text, strlen ( y_text ), 1, file );
        fwrite ( zi_id, 3, 1, file );
        fwrite ( z_text, strlen ( z_text ), 1, file );
        fwrite ( height_id, 3, 1, file );
        fwrite ( font_height, strlen ( font_height ), 1, file );
        fwrite ( text_id, 2, 1, file );
        fwrite ( text_text, strlen ( text_text ), 1, file );
}



template <typename T>
void DXFExport::processHalfEdges ( FILE * file, const Container <HalfEdge <T> *> &hl, const char * layer_name )
{
        //Process all half edges
        const unsigned int n = hl.size();

        for ( unsigned int i = 0; i < n; i++ )
        {
                // Get edge
                HalfEdge <T> *e = hl [i];

                //Use simplex indentificator to eliminate T processing of the edge
                if ( !e->isSimplexEdge() )
                {
                        // Set twin edge as simplex, i. e. processed
                        if ( e->getTwinEdge() ) e->getTwinEdge()->setEdgeAsSimplex ( true );

                        // Get start point
                        T x1 = e->getPoint()->getX();
                        T y1 = e->getPoint()->getY();
                        T z1 = ( ( Point3DCartesian <T> * ) e->getPoint() )->getZ();

                        // Get end point
                        T x2 = e->getNextEdge()->getPoint()->getX();
                        T y2 = e->getNextEdge()->getPoint()->getY();
                        T z2 = ( ( Point3DCartesian <T> * ) e->getNextEdge()->getPoint() )->getZ();

                        //Create line
                        createLine ( file, layer_name, x1, y1, z1, x2, y2, z2 );
                }
        }

        // Set all edges as non simplex
        for ( unsigned int i = 0; i < n; i++ )
        {
                hl [i] ->setEdgeAsSimplex ( false );
        }
}


template <typename T>
void DXFExport::processFace ( FILE * file, const Face <T> * face, const char * layer_name_edges )
{
        //Process one face
        HalfEdge <T> *e = face->getHalfEdge();

        //Previous half edge
        HalfEdge <T> *e_prev = NULL;

        //Remember start edge
        HalfEdge <T> *e_start = e;

        //Proces all edges of the Voronoi cell
        do
        {
                //Get ordering of edges in Voronoi cell
                e_prev = e;

                //Increment edge
                e = e->getNextEdge();

                // Get start point
                T x1 = e_prev->getPoint()->getX();
                T y1 = e_prev->getPoint()->getY();
                T z1 = 0;

                // Get end point
                T x2 = e->getPoint()->getX();
                T y2 = e->getPoint()->getY();
                T z2 = 0;

                //Create Voronoi edge
                createLine ( file, layer_name_edges, x1, y1, z1, x2, y2, z2 );

        }
        while ( e != e_start );
}


template <typename T>
void DXFExport::processVoronoiCells ( FILE * file, const Container <VoronoiCell <T> *> &vl, const char * layer_name_edges, const char * layer_name_generators, const bool process_all_cells )
{
        //Process all Voronoi cells
        const unsigned int n = vl.size();

        for ( unsigned int i = 0; i < n; i++ )
        {
                // Get edge
                VoronoiCell <T> *v = vl [i];

                //Process only bounded
                if ( ( process_all_cells ) || ( ( v->getBounded() ) && ( !process_all_cells ) ) )
                {
                        const Face <T> *face = dynamic_cast <Face <T> * const > ( v );

                        if ( face != NULL )
                        {
                                processFace ( file, face, layer_name_edges );
                        }
                }

                //Create generator as point
                T x = v->getGenerator()->getX();
                T y = v->getGenerator()->getY();
                T z = 0;

                //Create point
                createPoint ( file, layer_name_generators, x, y, z );
        }
}


template <typename GraticulePart, typename Point>
void DXFExport::processGraticuleElements ( FILE * file, GraticulePart & part, const Container <Point*> & points, const char * layer_graticule_name, const char * layer_labels_name,
                const typename Point::Type font_height, const typename Point::Type step )
{
        //Export meridian or parallel (defined as a template parameter GraticulePart)
        TIndexList ind_points = part.getPointsIndices();

        //Process all points
        const unsigned int n = ind_points.size();

        for ( unsigned int i = 1; i < n; i++ )
        {
                //Get actual point
                Point * p = points [ind_points[i]];

                //Get previous point
                Point * p_previous = points [ind_points[i - 1]];

                //Create a line> part of the meridian or parallel
                createLine ( file, layer_graticule_name, p_previous->getX(), p_previous->getY(), p_previous->getZ(), p->getX(), p->getY(), p->getZ() );
        }

        //Create label
        char point_id_text [255];

        if ( n > 1 )
        {
                //Set accuracy depending on a step
                if ( step > 1.0 ) sprintf ( point_id_text, "%3.1f", part.getCoord () );
                else if ( step > 0.1 ) sprintf ( point_id_text, "%3.2f", part.getCoord () );
                else sprintf ( point_id_text, "%3.3f", part.getCoord () );

                //Compute bearing
                const typename Point::Type bearing = Bearing::getBearing ( points [ind_points[0.5 * n - 1 ]], points [ind_points[0.5 * n]] );

                //Create label for meridian/parallel
                createText ( file, layer_labels_name, point_id_text, points[ind_points[0.5 * n - 1 ]]->getX() + 0.5 * font_height * cos ( bearing * M_PI / 180 ), points [ind_points[0.5 * n - 1]]->getY() + 0.5 * font_height * sin ( bearing * M_PI / 180 ), points [0]->getZ(), bearing, font_height );
        }
}


#endif
