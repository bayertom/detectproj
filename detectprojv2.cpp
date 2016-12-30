// Description: Estimation of the unknown cartographic projection from the map
// Version 2.0, C++11 support

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


// *******************************************************************************************************************************************
//
// ************************* RUNNING THE SCRIPT **********************************************************************************************

// 1] Basic: set method, test and reference files, detection/optimization methods 
//
//    detectprojv2 +met=nlsm7 test_points.txt reference.txt
//
//    available methods: nlsm7, nlsm8, nmm7, nmm8, dem7, dem8
//
// 		nlsm7		Non-linear least squares optimization, 7 determined parameters (M7)
// 		nlsm8		Non-linear least squares optimization, 8 determined parameters (M8), rotation involved
// 		nmm7		Nelder-Mead optimization, 7 determined parameters (M7)
// 		nmm8		Nelder-Mead optimization, 8 determined parameters (M8), rotation involved
// 		dem7		Differential evolution optimization, 7 determined parameters (M7)
// 		dem8		Differential evolution optimization, 8 determined parameters (M8), rotation involved

// 2] Set method, test and reference files, detection/optimization methods, meridian/parallel increments
// 
//    Important for DXF with the generated meridians/parallels; setting dlat/dlon increments of meridians/parallels
//
//    detectprojv2 +met=nlsm7 +dlat=10 +dlon=10 test_points.txt reference.txt
//
// 		dlat		Latitude step between two parallels (dlat >=1)
// 		dlon		Longitude step between two meridians (dlon>=1)

// 3] Set method, test and reference files, detection/optimization methods, meridian/parallel increments, amount of exported graticules to DXF
//
//    detectprojv2 +met=nlsm7 +dlat=10 +dlon=10 +gr=20 test_points.txt reference.txt
//
// 		gr		Amount of best-fit graticules exported to the DXF file (gr<=90)
//
//
// Example:
//
//		detectprojv2.exe +met=nlsm7 +dlat=10 +dlon=10 +gr=30 e:\maps\WorldMaps\Seutter\test.txt e:\maps\WorldMaps\Seutter\reference.txt 
//
// *******************************************************************************************************************************************



#include <vector>
#include <memory>
#include <ostream>
#include <iostream>
#include <fstream>
#include <array>

#include "libalgo/source/types/TVector2D.h"
#include "libalgo/source/types/TListS.h"
#include "libalgo/source/types/TInterval.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/graticule/Meridian.h"
#include "libalgo/source/structures/graticule/Parallel.h"

#include "libalgo/source/structures/projection2/Projection.h"
#include "libalgo/source/structures/projection2/Projections.h"
#include "libalgo/source/algorithms/cartanalysis2/CartAnalysis.h"
#include "libalgo/source/algorithms/round/Round.h"
#include "libalgo/source/algorithms/graticule3/Graticule.h"
#include "libalgo/source/algorithms/projectiontoproj4/ProjectionToProj42.h"


#include "libalgo/source/algorithms/planeintersection/PlaneIntersection.h"
#include "libalgo/source/algorithms/sphereintersection/SphereIntersection.h"
#include "libalgo/source/algorithms/greatcircleintersection/GreatCircleIntersection.h"

#include "libalgo/source/io2/DXFExport.h"

#include "libalgo/source/io2/File.h"
#include "libalgo/source/io/Format.h"

#include "libalgo/source/exceptions/FileReadException.h"

#ifndef CPP11_SUPPORT				//C++11 Support enabled
#define CPP11_SUPPORT				1
#endif


void printHelp(std::ostream & output)
{
	//Print help
	output << "\nParameters: \n\n";
	output << "  test.txt          Set a test file (unknown projection, cartesian coordinates x, y).\n";
	output << "  reference.txt     Set a reference file (geographic coordinates lat, lon).\n\n";
	output << "  -?                Help. \n\n";
	output << "Commands: \n\n";
	output << "  +met=             Set a method of analysis: gs, de, mls (default met=mls).\n";
	output << "  +gr=              The number exported graticules \n";
	output << "  Examples: \n\n";
	output << "  detectproj.exe +met=m6 vogt_test vogt_ref.txt \n";
	output << "  detectproj.exe +met=m7 +gr=10 vogt_test vogt_ref.txt \n";
}


int main(int argc, char * argv[])
{
	unsigned int exported_graticules = 60;
	bool reference_set = false;
	double lat_step = 10, lon_step = 10;

	std::string test_file, reference_file;
	TAnalysisMethod method = NLSM7;

	//Process the command line
	while (--argc > 0)
	{
		//Get - (A parameter follows)
		if (*argv[argc] == '-')
		{
			//Process parameter after -
			for (char * parameter = argv[argc]; ; )
			{
				switch (*(++parameter))
				{
					//Print help and cancel the analysis
				case '?':
				{
					printHelp(std::cout);
					return 0;
				}

				//Terminate character \0 of the argument
				case '\0':
					break;

					//Throw exception
				default:
					throw Exception("Exception: Invalid parameter in command line!");
				}

				//Stop processing of the argument
				break;

			}
		}

		//Get + ( A value follows)
		else if (*argv[argc] == '+')
		{
			//Get new command line parameter
			char * attribute = const_cast <char *> (argv[argc] + 1), *value = NULL;
			char * command = attribute;

			//Find splitter: =
			for (; (*command != '=') && (*command != '\0'); command++);

			//We found splitter, trim command and copy to value
			if ((*command == '=') && (*command != '\0'))
			{
				*command = '\0';
				value = command + 1;
			}

			//Throw exception
			if (attribute == NULL || value == NULL)
				throw Exception("Exception: Invalid value in command line!");

			//Set detection method
			if (!strcmp("met", attribute))
			{
				if (!strcmp("nlsm7", value)) method = NLSM7;
				else if (!strcmp("nlsm8", value)) method = NLSM8;
				else if (!strcmp("dem7", value)) method = DEM7;
				else if (!strcmp("dem8", value)) method = DEM8;
				else if (!strcmp("nmm7", value)) method = NMM7;
				else if (!strcmp("nmm8", value)) method = NMM8;
				else throw Exception("Exception: Invalid analysis method type in command line!");
			}

			//Set the latitude increment of adjacent meridians in the graticule
			else if (!strcmp("dlat", attribute))
			{
				lat_step = std::max(std::min(atof(value), 90.0), 0.1);

			}

			//Set the longitude increment of adjacent meridians in the graticule
			else if (!strcmp("dlon", attribute))
			{
				lon_step = std::max(std::min(atof(value), 180.0), 0.1); 
			}


			//Export graticule of the first k samples into DXF file
			else if (!strcmp("gr", attribute))
				exported_graticules = atof(value);

			//Bad argument
			else
			{
				std::cout << attribute << '\n';
				throw Exception("Exception: Invalid attribute in command line!");
			}
		}

		//Process input / output file
		else
		{
			//Set reference file
			if (!reference_set)
			{
				reference_file = argv[argc];
				reference_set = true;
			}
			//Set test file
			else if (reference_set)
				test_file = argv[argc];

			//Throw exception
			else 
				throw Exception("Exception: Too many input files in command line!");
		}
	}


	//List of projections
	TListS <Projection<double>> projections;

	//Initialize projections
	Projections::init(projections);


	//Manually set the file
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\v2\\M8\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\v2\\M8\\reference.txt";


	try
	{
		//Load list of test points
		TVector <Point3DCartesian<double>> test_points;
		File::load2DPoints<Point3DCartesian<double>>(test_file.c_str(), test_points);

		//Load list of reference points
		TVector <Point3DGeographic<double>> reference_points;
		File::load2DPoints<Point3DGeographic<double>>(reference_file.c_str(), reference_points);

		//Different amount of points
		if (test_points.size() != reference_points.size())
			throw BadDataException("BadDataException: different amount of points in the test/reference files. ", "Analysis canceled.");

		//Not enough points
		if (test_points.size() < 4)
			throw BadDataException("BadDataException: not enough points (n < 4). ", "Analysis canceled.");


		//Set output file
		std::ostream * output = &std::cout;
		static std::ofstream output_file;

		//Create output protocol
		std::string output_file_text;
		output_file_text += test_file + ".log";
		output_file_text = "output.log";						//DEBUG
		output_file.open(output_file_text, std::fstream::app );

		output = &output_file;
		std::cout << "\nMap projection analysis: method = " << (int)method << "\n\n";
		*output << "\nMap projection analysis: method = " << (int)method << "\n\n";
		*output << test_file << '\n';							//DEBUG

		//Analyze projection
		TResults <double> results;
		CartAnalysis::analyzeProjection(test_points, reference_points, projections, results, method, *output);

		//Print results
		CartAnalysis::printResults(results, std::cout);
		CartAnalysis::printResults(results, *output);

		output_file.close();

		//Create graticules
		if (exported_graticules > 0)
		{
			std::cout << ">> Exporting points, graticules, creating Proj.4 strings. Please wait... \n";

			// Geographic extent of the analyzed territory
			double lat_min = (std::min_element(reference_points.begin(), reference_points.end(), sortPointsByLat()))->getLat();
			double lat_max = (std::max_element(reference_points.begin(), reference_points.end(), sortPointsByLat()))->getLat();
			double lon_min = (std::min_element(reference_points.begin(), reference_points.end(), sortPointsByLon()))->getLon();
			double lon_max = (std::max_element(reference_points.begin(), reference_points.end(), sortPointsByLon()))->getLon();

			//Get limits; stretch over the whole planishere, if necessarry
			if (((lon_min < MIN_LON + 20) || (lon_max > MAX_LON - 20)) && (lon_max - lon_min > 200))
			{
				lon_min = MIN_LON;
				lon_max = MAX_LON;
			}
			TInterval <double> lat_interval{ lat_min, lat_max };
			TInterval <double> lon_interval{ lon_min, lon_max };

			//Change step of  meridians/parallels
			const double dlat = lat_max - lat_min;
			const double dlon = lon_max - lon_min;
			const double lat_step = (dlat < 20.0 ? (dlat < 2.0 ? 0.1 : 1.0) : 10.0);
			const double lon_step = (dlon < 20.0 ? (dlon < 2.0 ? 0.1 : 1.0) : 10.0);

			//Create graticules
			unsigned int index = 0;
			for (auto res : results)
			{
				TVector <Meridian <double> > meridians;
				TVector <Parallel < double> > parallels;
				TVector2D <Point3DCartesian <double> > meridians_proj, parallels_proj;

				//Set font height
				const double font_height = 0.05 * res.second.proj->getR() * std::min(lat_step, lon_step) * M_PI / 180;

				//Get map rotation
				double alpha = res.second.map_rotation;

				//Create graticule
				Graticule::createGraticule(res.second.proj, lat_interval, lon_interval, lat_step, lon_step, 0.1 * lat_step, 0.1 * lon_step, alpha, meridians, meridians_proj, parallels, parallels_proj);

				//Stop export
				if (++index > exported_graticules)
					break;

				//Create file names
				std::string output_file_graticule = test_file + "_" + std::to_string(index) + "_" + res.second.proj->getName();
				std::string output_file_proj4 = output_file_graticule;
				output_file_graticule += "_grat.dxf";
				output_file_proj4 += "_proj4.bat";

				//Export graticule to DXF
				TVector <Point3DCartesian<double> > reference_points_proj;
				CartTransformation::latsLonstoXY(reference_points, res.second.proj, alpha, reference_points_proj);
				DXFExport::exportGraticuleToDXF(output_file_graticule.c_str(), meridians, meridians_proj, parallels, parallels_proj, test_points, reference_points_proj,font_height, std::min(lat_step, lon_step));
				std::cout << ".";

				//Create Proj.4 definition string
				std::string proj4_string = ProjectionToProj4::ProjectionToProj4String(res.second.proj);
				std::ofstream output_file(output_file_proj4);
				output_file << proj4_string;
				output_file.close();
			}
		}
	}

	//Some error occurs
	catch (Exception &e)
	{
		e.printException();
	}

	return 0;
}