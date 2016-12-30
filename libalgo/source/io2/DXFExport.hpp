// Description: Export lines, points, polygons to 2D/3D DXF file

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



#ifndef DXFExport_HPP
#define DXFExport_HPP

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>


#include "libalgo/source/const2/Const.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/exceptions/FileWriteException.h"


template <typename T>
void DXFExport::exportGraticuleToDXF(const std::string &file_name, const TVector <Meridian <T> > &meridians, const TVector2D <Point3DCartesian<T> > & meridians_proj, const TVector <Parallel <T> > &parallels, const TVector2D <Point3DCartesian<T> > & parallels_proj, const TVector <Point3DCartesian<T> > & test_points, const TVector <Point3DCartesian<T> > & reference_points_proj, const T font_height, const T step)
{
	//Export the graticule formed by meridians/parallels into DXF file
	//Export test points, reference points, and projected points
	const unsigned int color_graticule = 8, color_test_points = 5, color_reference_points = 1;

	const std::string  level_meridians = "Meridians", level_parallels = "Parallels", level_meridian_labels = "Meridian_labels",
		level_parallel_labels = "Parallel_labels", level_test_points = "Test_points", level_test_point_labels = "Test_point_labels",
		level_proj_reference_points = "Reference_points_proj", level_proj_reference_point_labels = "Reference_points_proj_labels";
	
	std::ofstream file;

	try
	{
		file.open(file_name/*, ios::out*/);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for meridians
			createLayerSection(file, level_meridians, color_graticule);

			//Create layer for parallels
			createLayerSection(file, level_parallels, color_graticule);

			//Create layer for meridian labels
			createLayerSection(file, level_meridian_labels, color_graticule);

			//Create layer for parallel labels
			createLayerSection(file, level_parallel_labels, color_graticule);

			//Create layer for test points
			createLayerSection(file, level_test_points, color_test_points);

			//Create layer for test points
			createLayerSection(file, level_test_point_labels, color_test_points);

			//Create layer for projected reference points
			createLayerSection(file, level_proj_reference_points, color_reference_points);

			//Create layer for projected reference points
			createLayerSection(file, level_proj_reference_point_labels, color_reference_points);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all meridians
			const unsigned int nm = meridians_proj.size();
			for (unsigned int i = 0; i < nm; i++)
			{
				processGraticuleElements(file, meridians_proj[i], meridians[i].getLon(), level_meridians, level_meridian_labels, font_height, step, color_graticule);
			}

			//Process all parallels
			const unsigned int np = parallels_proj.size();
			for (unsigned int i = 0; i < np; i++)
			{
				processGraticuleElements(file, parallels_proj[i], parallels[i].getLat(), level_parallels, level_parallel_labels, font_height, step, color_graticule);
			}
			
			//Process all test points
			unsigned int index_test = 0;
			for (const Point3DCartesian <T> point : test_points)
			{
				//Create test point
				createPoint(file, level_test_points, point.getX(), point.getY(), point.getZ(), color_test_points);

				//Create test point label
				std::string point_label = std::to_string(++index_test);
				createText(file, level_test_point_labels, point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0.0, font_height, color_test_points);
			}

			//Process all projected reference points
			unsigned int index_reference = 0;
			for (const Point3DCartesian <T> point : reference_points_proj)
			{
				//Create reference point
				createPoint(file, level_proj_reference_points, point.getX(), point.getY(), point.getZ(), color_reference_points);

				//Create reference point label
				std::string point_label = std::to_string( ++index_reference);
				createText(file, level_proj_reference_point_labels, point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0.0, font_height, color_reference_points);
			}
			
			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


template <typename Point, typename T>
void DXFExport::exportPointsToDXF(const std::string &file_name, const TVector <Point> &points, const T font_height, const unsigned int color)
{
	//Export list of points to DXF file
	std::ofstream file;

	try
	{
		file.open(file_name);

		if (file.is_open())
		{
			//Create header section
			createHeaderSection(file);

			//Create table section
			createTableSection(file);

			//Create layer for points
			createLayerSection(file, "Points", color);

			//Create layer for point labels
			createLayerSection(file, "Point_labels", color);

			//End table header
			endTableSection(file);

			//Create entity section
			createEntitySection(file);

			//Process all points
			for (const Point point:points)
			{
				//Create point
				createPoint(file, "Points", point.getX(), point.getY(), point.getZ(), color);

				//Create point label
				char point_id_text[255];
				sprintf(point_id_text, "%d", point.getPointID());
				std::string point_label = point_id_text;
				createText(file, "Point_labels", point_label, point.getX() + 0.5 * font_height, point.getY() - 0.5 * font_height, point.getZ(), 0, font_height, color);
			}

			//End header section
			endHeaderSection(file);

			//Close file
			file.close();
		}

		//Throw exception
		else
		{
			//Can not open file
			throw std::ios_base::failure("Exception: can not open the file. ");
		}
	}

	//Any error has appeared
	catch (std::ios_base::failure &)
	{
		//Close file
		file.close();

		//Throw exception
		throw FileWriteException("FileWriteException: can not write the file: ", file_name);
	}
}


inline void DXFExport::createHeaderSection ( std::ofstream & file )
{
        //Create header section
	const std::string object_type_id = "0\n";
	const std::string object_type = "SECTION\n";
	const std::string header_id = "2\n";
	const std::string header_name = "HEADER\n";
	const std::string variable_name = "9\n";
	const std::string acad_version = "$ACADVER\n";
	const std::string acad_version_id =  "1\n";
	const std::string acad_version_val = "AC1006\n";
	const std::string section_terminator = "ENDSEC\n";
        
	/* Add to file */
	file << object_type_id;
        file << object_type;
        file << header_id;
        file << header_name;
        file << variable_name;
        file << acad_version;
        file << acad_version_id;
        file << acad_version_val;
        file << object_type_id;
        file << section_terminator;
}


inline void DXFExport::endHeaderSection (std::ofstream & file)
{
        //Create end of the header section
	const std::string entity_id = "0\n";
	const std::string section_terminator ="ENDSEC\n";
	const std::string file_terminator = "EOF\n";
	
	/* Add to file */
	file << entity_id;
        file << section_terminator;
        file << entity_id;
        file << file_terminator;
}


inline void DXFExport::createTableSection (std::ofstream & file )
{
        //Create table section
	const std::string object_type_id = "0\n";
	const std::string object_name = "SECTION\n";
	const std::string table_id = "2\n";
	const std::string table = "TABLES\n";
	const std::string table_name = "TABLE\n";
	const std::string layer_name = "LAYER\n";
	const std::string max_number_of_entries_id = "70\n";
	const std::string max_number_of_entries = "0\n";

        file << object_type_id;
        file << object_name;
        file << table_id;
        file << table;
        file << object_type_id;
        file << table_name;
	file << table_id;
        file << layer_name;
        file << max_number_of_entries_id;
        file << max_number_of_entries;
}


inline void DXFExport::endTableSection ( std::ofstream & file )
{
        //Write end of the table section
	const std::string object_type_id = "0\n";
	const std::string table_end_name = "ENDTAB\n";
	const std::string section_end_name = "ENDSEC\n";

	/* Add to file */
        file << object_type_id;
        file << table_end_name;
        file << object_type_id;
        file << section_end_name;
}


inline void DXFExport::createLayerSection (std::ofstream & file, const std::string &layer_name, const unsigned int color )
{
        //Add section for one layer
	const std::string object_type_id = "0\n";
	const std::string object_name = "LAYER\n";
	const std::string layer_name_id = "2\n";
	const std::string layer_flag_id = "70\n";
	const std::string layer_flag = "0\n";
	const std::string color_number_id = "62\n";
	const std::string line_type_id = "6\n";
	const std::string line_type_name = "CONTINUOUS\n";

	/* Add to file */
        file << object_type_id;
        file << object_name;
        file << layer_name_id;
	file << layer_name << '\n';
        file << layer_flag_id;
        file << layer_flag;
        file << color_number_id;
        file << color << '\n';
        file << line_type_id;
        file << line_type_name;
}


inline void DXFExport::createEntitySection ( std::ofstream & file )
{
        //Create section for entities
	const std::string object_type_id = "0\n";
	const std::string object_name = "SECTION\n";
	const std::string entity_name_id = "2\n";
	const std::string entity_name = "ENTITIES\n";

	/* Add to file */
        file << object_type_id;
        file << object_name;
        file << entity_name_id;
        file << entity_name;
}


template <typename T>
void DXFExport::createLine (std::ofstream & file, const std::string &layer_name, const T x1, const T y1, const T z1, const T x2, const T y2, const T z2 )
{
        //Write line to DXF file
	const std::string entity_id = "0\n";
	const std::string entity_name = "LINE\n";
	const std::string level_id = "8\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";
	const std::string xii_id = "11\n";
	const std::string yii_id = "21\n";
	const std::string zii_id = "31\n";

        /* Add to file */
        file << entity_id;
        file << entity_name;
        file << level_id;
	file << layer_name << '\n';
        file << xi_id;
        file << x1 << '\n';
	file << yi_id;
	file << y1 << '\n';
        file << zi_id;
        file << z1 << '\n';
        file << xii_id;
        file << x2 << '\n';
        file << yii_id;
        file << y2 << '\n';
        file << zii_id;
        file << z2 << '\n';
}


template <typename T>
void DXFExport::createPoint (std::ofstream & file, const std::string & layer_name, const T x, const T y, const T z, const int color )
{
        //Write point to DXF file
	const std::string entity_id = "0\n";
	const std::string entity_name = "POINT\n";
	const std::string level_id = "8\n";
	const std::string color_id = "62\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";

        /* Add to file */
        file << entity_id;
        file << entity_name;
        file << level_id;
        file << layer_name << '\n';
        file << color_id;
        file << color << '\n';
        file << xi_id;
        file << x << '\n';
        file << yi_id;
        file << y << '\n';
        file << zi_id;
        file << z << '\n';
}


template <typename T>
void DXFExport::createText (std::ofstream & file, const std::string & layer_name, const std::string & text, const T x, const T y, const T z, const T rotation, const T height, const unsigned int color )
{
        //Create text
	const std::string entity_id = "0\n";
	const std::string entity_name = "TEXT\n";
	const std::string style_id = "7\n";
	const std::string text_style = "PNTNUM\n";
	const std::string rotation_id = "50\n";
	const std::string level_id = "8\n";
	const std::string color_id = "62\n";
	const std::string xi_id = "10\n";
	const std::string yi_id = "20\n";
	const std::string zi_id = "30\n";
	const std::string height_id = "40\n";
	const std::string text_id = "1\n";

        /* Add to file */
        file << entity_id;
        file << entity_name;
        file << style_id;
        file << text_style;
        file << rotation_id;
        file << rotation << '\n';
        file << level_id;
        file << layer_name << '\n';
        file << color_id;
        file << color << '\n';
        file << xi_id;
        file << x << '\n';
        file << yi_id;
        file << y << '\n';
        file << zi_id;
        file << z << '\n';
        file << height_id;
        file << height << '\n';
        file << text_id;
        file << text << '\n';
}


template <typename GraticulePart, typename T>
void DXFExport::processGraticuleElements( std::ofstream & file, GraticulePart & part, const T val, const std::string &layer_graticule_name, const std::string &layer_labels_name, const T font_height, const T step, const unsigned int color)
{
	//Export meridian or parallel (defined as a template parameter GraticulePart)

	//Process all points
	const unsigned int n = part.size();

	for (unsigned int i = 1; i < n; i++)
	{
		//Get actual point
		Point3DCartesian<T> p = &part[i];

		//Get previous point
		Point3DCartesian<T>  p_previous = &part[i - 1];

		//Create a line> part of the meridian or parallel
		createLine(file, layer_graticule_name, p_previous.getX(), p_previous.getY(), p_previous.getZ(), p.getX(), p.getY(), p.getZ());
	}

	//Create label
	char point_id_char[255];

	if (n > 1)
	{
		//Set accuracy depending on a step
		if (step > 1.0) sprintf(point_id_char, "%3.1f", val);
		else if (step > 0.1) sprintf(point_id_char, "%3.2f", val);
		else sprintf(point_id_char, "%3.3f", val);

		std::string point_id_text = point_id_char;

		//Compute bearing
		Point3DCartesian<T>  p1 = &part[0.5 * n - 1], p2 =  &part[0.5 * n];
		const T bearing = atan2(p2.getY() - p1.getY(), p2.getX() - p1.getX()) * RO;

		//Create label for meridian/parallel
		createText(file, layer_labels_name, point_id_text, part[0.5 * n - 1].getX() + 0.5 * font_height * cos(bearing * M_PI / 180), 
			part[0.5 * n - 1].getY() + 0.5 * font_height * sin(bearing * M_PI / 180), part[0].getZ(), bearing, font_height, color);
	}
}


template <typename T>
std::string DXFExport::to_string(const T value, const unsigned short dec_places)
{
	//Convert to string,  set the decimal places
	std::ostringstream out;
	out << std::setprecision(dec_places) << value;
	return out.str();
}

#endif
