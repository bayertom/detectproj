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
	/*
	std::array<std::array<double, 3>, 2> items = { {
	{ 1, 2, 3 },
	{ 4, 5, 6 },
	} };

	Matrix <double> A(items);
	A.print();

	std::array<std::array<double, 2>, 2> items2 = { {
		{ 10, 11 },
		{ 20, 21 },
		} };

	Matrix <double>B(items2);
	B.print();

	A(B, 1, 1);
	A.print();
	*/
	//Matrix <double> E = pinv1(D);
	//E.print();
	/*
	double latp = 50.0, lonp = 15.0;

	double lat1 = 51, lon1 = 16, lat2 = 51, lon2 = 14, lat3 = 49, lon3 = 14, lat4 = 49, lon4 = 16;

	TTransformedLongitudeDirection mode = NormalDirection2;

	double latt1 = CartTransformation::latToLatTrans(lat1, lon1, latp, lonp);
	double latt2 = CartTransformation::latToLatTrans(lat2, lon2, latp, lonp);
	double latt3 = CartTransformation::latToLatTrans(lat3, lon3, latp, lonp);
	double latt4 = CartTransformation::latToLatTrans(lat4, lon4, latp, lonp);

	double lont1 = CartTransformation::lonToLonTrans(lat1, lon1, latp, lonp, mode);
	double lont2 = CartTransformation::lonToLonTrans(lat2, lon2, latp, lonp, mode);
	double lont3 = CartTransformation::lonToLonTrans(lat3, lon3, latp, lonp, mode);
	double lont4 = CartTransformation::lonToLonTrans(lat4, lon4, latp, lonp, mode);

	double lat21 = CartTransformation::latTransToLat(latt1, lont1, latp, lonp, mode);
	double lat22 = CartTransformation::latTransToLat(latt2, lont2, latp, lonp, mode);
	double lat23 = CartTransformation::latTransToLat(latt3, lont3, latp, lonp, mode);
	double lat24 = CartTransformation::latTransToLat(latt4, lont4, latp, lonp, mode);

	double lon21 = CartTransformation::lonTransToLon(latt1, lont1, latp, lonp, mode);
	double lon22 = CartTransformation::lonTransToLon(latt2, lont2, latp, lonp, mode);
	double lon23 = CartTransformation::lonTransToLon(latt3, lont3, latp, lonp, mode);
	double lon24 = CartTransformation::lonTransToLon(latt4, lont4, latp, lonp, mode);

	double xxx = 62;
	*/
	/*
	double x1 = 0;
	double y1 = 0;
	double z1 = 0;

	double x2 = 4;
	double y2 = -1;
	double z2 = 1;
	
	double x3 = 4;
	double y3 = -1;
	double z3 = -1;

	double x4 = 0;
	double y4 = 0;
	double z4 = 0;

	double x5 = 4;
	double y5 = -1;
	double z5 = -1;

	double x6 = 4;
	double y6 = 1;
	double z6 = -1;

	double x0 = 8.0/3;
	double y0 = -1.0/3;
	double z0 = -1.0/3;*/

	/*
	double x1 = 0;
	double y1 = 0;
	double z1 = 1;

	double x2 = 1;
	double y2 = 0;
	double z2 = 0;

	double x3 = 0;
	double y3 = 1;
	double z3 = 0;

	double x4 = 0;
	double y4 = 0;
	double z4 = 0;

	double x5 = 1;
	double y5 = 2;
	double z5 = 0;

	double x6 = 0;
	double y6 = 0;
	double z6 = 1;

	double x0 = 0;
	double y0 = 0;
	double z0 = 0;


	double xi, yi, zi, tx, ty, tz;
	//bool res = PlaneIntersection::get2PlanesIntersection(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x0, y0, z0, xi, yi, zi, tx, ty, tz);

	double x0x = 37;

	double xc = 0, yc = 0, zc = 0, r = 2;
	double xa = 0, ya = 0, za = 0, ux = 1.0, uy = 1.0, uz = 1.0;
	double xi1, yi1, zi1, xi2, yi2, zi2;

	bool res2 = SphereIntersection::getSphereAndLineIntersection(xc, yc, zc, r, xa, ya, za, ux, uy, uz, xi1, yi1, zi1, xi2, yi2, zi2);


	Point3DGeographic<double> p1(90, 0), p2(-90, 0), p3(0, 90);
	Point3DGeographic<double> p4(10, 0), p5(10, 90), p6(10, 180);
	Point3DGeographic<double> i1, i2, pole(50,15);
	TTransformedLongitudeDirection lon_direction = NormalDirection2;
	GreatCircleIntersection::getGreatCirclePlainIntersection(p1, p2, p3, p4, p5, p6, i1, i2, pole, lon_direction);
	
	Matrix <double> A(3, 3);
	A(0, 0) = 1;
	A(0, 1) = 2;
	A(0, 2) = 3;
	A(2, 0) = 1;
	A(2, 1) = 0;
	A(2, 2) = 1;
	A(1, 0) = 4;
	A(1, 1) = 5;
	A(1, 2) = 6;
	double dd = MatrixOperations::det(A);


	x0x = 38;

	*/
	/*
	double latp = 89.771283020885448, lonp = 109.99618390644875;

	double lat = 0, lon1 = 105, lon2 = 115, lon3 = 125;

	double lt1 = CartTransformation::lonToLonTrans(lat, lon1, latp, lonp, NormalDirection2);
	double lt2 = CartTransformation::lonToLonTrans(lat, lon2, latp, lonp, NormalDirection2);
	double lt3 = CartTransformation::lonToLonTrans(lat, lon3, latp, lonp, NormalDirection2);
	
	double value = 55, multiple = 2.5;
	double xdown, xup;
	xdown = Round::roundToMultipleFloor(value, multiple) + multiple;
	xup = Round::roundToMultipleCeil(value, multiple) - multiple;

	//double xup = std::ceil(std::ceil(value / multiple) * multiple);
	/*
	TInterval<double> lat_interval1{ -86, -33 };
	TInterval<double> lat_interval2{ -86, -30 };
	TInterval<double> lat_interval3{ -80, -33 };

	TInterval<double> lon_interval1{ -176, -133 };
	TInterval<double> lon_interval2{ -176, -130 };
	TInterval<double> lon_interval3{ -170, -133 };

	Meridian <double> m1(15, lat_interval1, 10, 0, 0);
	Meridian <double> m2(15, lat_interval2, 10, 0, 0);
	Meridian <double> m3(15, lat_interval3, 10, 0, 0);

	Parallel <double> p1(15, lon_interval1, 10, 0, 0);
	Parallel <double> p2(15, lon_interval2, 10, 0, 0);
	Parallel <double> p3(15, lon_interval3, 10, 0, 0);
	*/
	/*
	
	TVector <Meridian <double> > meridians2;
	TVector <Parallel < double> > parallels2;
	TVector2D <Point3DCartesian <double> > meridians_proj2, parallels_proj2;

	TVector <Point3DCartesian <double> > test_points, reference_points_proj;

	//Set font height
	const double font_height2 = 100;
	const double lat_step2 = 15, lon_step2 = 30;
	const double alpha2 = 0.0;
	TInterval <double> lat_interval2{15,90 }, lon_interval2{ -180, 180 };

	//List of projections
	TListS <Projection<double>> projs;

	//Initialize projections
	Projections::init(projs);
	
	for (int i = 0; i < 1; i++)
	{
		std::shared_ptr <Projection<double>> proj2 = *(projs.begin());
		auto i_proj2 = projs.begin();
		Point3DGeographic <double> cart_pole(90, 0);
		//Point3DGeographic <double> cart_pole(50, 15);
		//Point3DGeographic <double> cart_pole(89.5, -174.48);
		(*i_proj2)->setLat1(0.0);
		//(*i_proj2)->setLon0(10);
		//(*i_proj2)->setLon0(-122.71);

		(*i_proj2)->setCartPole(cart_pole);


		double XX1 = (*i_proj2)->getX((*i_proj2)->getR(), (*i_proj2)->getLat1(), (*i_proj2)->getLat2(), 50, 15, 0, (*i_proj2)->getDx(), (*i_proj2)->getDy(), (*i_proj2)->getC());
		double YY1 = (*i_proj2)->getY((*i_proj2)->getR(), (*i_proj2)->getLat1(), (*i_proj2)->getLat2(), 50, 15, 0, (*i_proj2)->getDx(), (*i_proj2)->getDy(), (*i_proj2)->getC());

		/*
		double lat = 0;
		double lon = 150;// 89.999999999999801;
		double R = 6380;
		double lat1 = 10;


		double F = ((0.25*M_PI*M_PI - (lat / RO) *(lat / RO)) / (M_PI*abs(sin(lat / RO)) - 2 * abs(lat) / RO));
		double L = ((2 * lon / (M_PI*RO))*(2 * lon / (M_PI*RO)) - 1);

		double Y1 = R / L*lat / abs(lat)*(sqrt(F*F - L*(0.25*M_PI*M_PI - F*M_PI*abs(sin(lat / RO)) - (lon / RO)* (lon / RO))) - F);
		double X1 = R*lon / RO * sqrt(1.0 - (2.0 * Y1 / (M_PI * R)) * (2.0 * Y1 / (M_PI * R)));

		double X2 = R*sqrt(2.0 / ((1 + sin(lat1 / RO)) / 2)*(1 - sin(lat / RO)))*sin((1 + sin(lat1 / RO)) / 2 * lon / RO);
		double Y2 = 2 * R*sqrt ( (1 - sin(lat1 / RO)) / (1 + sin(lat1 / RO))) - R*sqrt(2 / ((1 + sin(lat1 / RO)) / 2)*(1 - sin(lat / RO)))*cos((1 + sin(lat1 / RO)) / 2 * lon / RO);
		
		//Reduce longitude
		double lonr = CartTransformation::redLon0(lon, (*i_proj2)->getLon0());

		double lat_trans = CartTransformation::latToLatTrans(lat, lonr, cart_pole.getLat(), cart_pole.getLon());
		double lon_trans = CartTransformation::lonToLonTrans(lat, lonr, cart_pole.getLat(), cart_pole.getLon(), NormalDirection);

		double XX1 = (*i_proj2)->getX((*i_proj2)->getR(), (*i_proj2)->getLat1(), (*i_proj2)->getLat2(), lat_trans, lon_trans, 0, (*i_proj2)->getDx(), (*i_proj2)->getDy(), (*i_proj2)->getC());
		double YY1 = (*i_proj2)->getY((*i_proj2)->getR(), (*i_proj2)->getLat1(), (*i_proj2)->getLat2(), lat_trans, lon_trans, 0, (*i_proj2)->getDx(), (*i_proj2)->getDy(), (*i_proj2)->getC());

		double XX2 = (*i_proj2)->getX(lat, lon+10);
		double YY2 = (*i_proj2)->getY(lat, lon+10);
		double XX3 = (*i_proj2)->getX(-lat, -lon);
		double YY3 = (*i_proj2)->getY(-lat, -lon);
		double XX4 = (*i_proj2)->getX(lat, -lon);
		double YY4 = (*i_proj2)->getY(lat, -lon);

		lat = lat * M_PI / 180;
		lon = lon * M_PI / 180;
		lat1 = 40 * M_PI / 180;
		double lat2 = 50 * M_PI / 180;
		double XXX = 2 * R / (sin(lat1) + sin(lat2))*sqrt(cos(lat1)*cos(lat1) + (sin(lat1) + sin(lat2))*(sin(lat1) - sin(lat)))*sin((sin(lat1) + sin(lat2)) / 2 * lon);
		double YYY = 2 * R / (sin(lat1) + sin(lat2))*sqrt(cos(lat1)*cos(lat1) + (sin(lat1) + sin(lat2))*(sin(lat1) - sin((lat1 + lat2) / 2))) - 2 * R / (sin(lat1) + sin(lat2))*sqrt(cos(lat1)*cos(lat1) + (sin(lat1) + sin(lat2))*(sin(lat1) - sin(lat)))*cos((sin(lat1) + sin(lat2)) / 2 * lon);
		*/
		
		//Create graticule
		Graticule::createGraticule(*i_proj2, lat_interval2, lon_interval2, lat_step2, lon_step2, 1.0, 1.0, alpha2, meridians2, meridians_proj2, parallels2, parallels_proj2);

		//Create file names
		int index2 = 0;
		std::string test_file2 = "testx.txt";
		std::string output_file_graticule2 = test_file2;
		std::string output_file_points_test2 = output_file_graticule2;  //Copy string without ID and projection name
		output_file_graticule2 += "_" + std::to_string(index2 + 1) + "_" + proj2->getName() + "_proj_" + std::to_string(i);
		std::string output_file_points_ref2 = output_file_graticule2;
		std::string output_file_proj42 = output_file_graticule2;
		output_file_graticule2 += "_grat.dxf";
		output_file_points_ref2 += "_points_ref.dxf";
		output_file_points_test2 += "_points_test.dxf";
		output_file_proj42 += "_proj4.bat";

		//Export graticule to DXF
		//const std::string &file_name, const TVector <Meridian <T> > &meridians, const TVector2D <Point3DCartesian<T> > & meridians_proj, const TVector <Parallel <T> > &parallels, const TVector2D <Point3DCartesian<T> > & parallels_proj, const TVector <Point3DCartesian<T> > & test_points, const TVector <Point3DCartesian<T> > & reference_points_proj, const T font_height, const T step
		//DXFExport::exportGraticuleToDXF(output_file_graticule.c_str(), meridians, meridians_proj, parallels, parallels_proj, test_points, reference_points_proj, font_height, std::min(lat_step, lon_step));

		DXFExport::exportGraticuleToDXF(output_file_graticule2.c_str(), meridians2, meridians_proj2, parallels2, parallels_proj2,test_points, reference_points_proj, lat_step2, font_height2);

		meridians2.clear();
		parallels2.clear();
		meridians_proj2.clear();
		parallels_proj2.clear();

		std::cout << i << '\n';
	}
	*/
	/*
	Matrix <double> A(3,3);
	A(0, 0) = 1;
	A(0, 1) = 2;
	A(0, 2) = 3;
	A(1, 0) = 4;
	A(1, 1) = 5;
	A(1, 2) = 6;
	A(2, 0) = 0;
	A(2, 1) = 0;
	A(2, 2) = 0;

	std::array<std::array<double, 6>, 6> items = { {

		{ 1.3757475400518885E-5, -2.847157448783027E-6, 0.0, 0.0, -2.8471620691817876E-6, 0.0 },
		{ -2.847157448783027E-6, 3.0588500226423854E-5, 0.0, 0.0, 3.058852823111269E-5, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ -2.8471620691817876E-6, 3.058852823111269E-5, 0.0, 0.0, 3.058855623584992E-5, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }

	} };

	Matrix <double> D(items);
	D.print();
	
	Matrix <double> E = pinv1(D);
	E.print();
	*/
	//**********************************************************************************************************

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

	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Mercator\\v2\\M7S\\test2.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Mercator\\v2\\M7S\\reference.txt";
	
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\AirRoutes\\v2\\M8\\NLS\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\AirRoutes\\v2\\M8\\NLS\\reference.txt";

	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\WorldAround\\v2\\M8\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\WorldAround\\v2\\M8\\reference.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\v2\\M8\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\v2\\M8\\reference.txt";
	
	test_file = "E:\\Tomas\\Java\\detectprojv2j\\test\\Others\\test_points_title_grat2.txt";
	reference_file = "E:\\Tomas\\Java\\detectprojv2j\\test\\Others\\reference_points_title_grat2.txt";

	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\Europe\\test_points_grat.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\Europe\\reference_points_grat.txt";


	//(M7, M8: NLS fail)
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Ortelius\\v2\\M8\\NLS\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Ortelius\\v2\\M8\\NLS\\reference.txt";

	//(M7, M8: NLS fail; 350 iterations required)
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Agriculture\\v2\\M8\\NLS\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Agriculture\\v2\\M8\\NLS\\reference.txt";
	
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Reichard\\v2\\M7S\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Reichard\\v2\\M7S\\reference.txt";
	
	//M7: NLS fail
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Ortelius\\v2\\M8\\NLS\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Ortelius\\v2\\M8\\NLS\\reference.txt";
	/*
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\DeLisle\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\DeLisle\\v2\\M7S\\reference.txt";
	
	//?
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Northern\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Northern\\v2\\M7S\\reference.txt";
	//?
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Hondius_world\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Hondius_world\\v2\\M7S\\reference.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Colton_world\\v2\\M7S\\testw.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Colton_world\\v2\\M7S\\referencew.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Colton_world\\v2\\M7S\\teste.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Colton_world\\v2\\M7S\\referencee.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Africa\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Africa\\v2\\M7S\\reference.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Asia\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Asia\\v2\\M7S\\reference.txt";
	
	//M8: NLS fail
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Europe\\v3\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Europe\\v3\\M7S\\reference.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Noah\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Noah\\v2\\M7S\\reference.txt";
	
	test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\SouthAfrica\\v2\\M7S\\test.txt";
	reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\SouthAfrica\\v2\\M7S\\reference.txt";
	*/
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\Safarik\\v2\\M7S\\test.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\Safarik\\v2\\M7S\\reference.txt";
	
	//M7, M8: NLS fail
	///test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Problems\\Mercator\\v2\\M7S\\testn.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Problems\\Mercator\\v2\\M7S\\referencen.txt";
	
	//test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Problems\\Mercator\\v2\\M7S\\tests.txt";
	//reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Problems\\Mercator\\v2\\M7S\\references.txt";
	
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


		/*
		//(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);
		std::shared_ptr <Projection<double>> proj = *(projections.begin());

		auto i_proj = projections.begin();

		//Jump to cea
		std::advance(i_proj, 9);

		p_coord_function <double> getX = (*i_proj)->getX();

		double lat = 50, lon = 15;
		Matrix <double> X(7, 1);
		X(0, 0) = 6380;
		X(1, 0) = 90;
		X(2, 0) = 0;
		X(3, 0) = 45;
		X(4, 0) = 45;
		X(5, 0) = 0;
		X(6, 0) = 0;

		Matrix <double> XT = trans(X);
		/*
		double fR = NumDerivative::getDerivative(FDiffM7<double>(lat, lon, getX, ReversedDirection), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, false);
		double flat1 = NumDerivative::getDerivative(FDiffM7<double>(lat, lon, getX, ReversedDirection), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, false);
		double flat2 = NumDerivative::getDerivative(FDiffM7<double>(lat, lon, getX, ReversedDirection), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, false);
		double flon0 = NumDerivative::getDerivative(FDiffM7<double>(lat, lon, getX, ReversedDirection), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, false);
		*/
		/*
		Matrix <double> X2(6, 1);
		const double R = 6380;
		X2(0, 0) = 90;
		X2(1, 0) = 0;
		X2(2, 0) = 45;
		X2(3, 0) = 45;
		X2(4, 0) = 0;
		X2(5, 0) = 0;

		Matrix <double> XT2 = trans(X2);
		/*
		double flat12 = NumDerivative::getDerivative(FDiffM8<double>(X(0,0), lat, lon, getX, ReversedDirection), XT2, FirstDerivative, VariableX3, NUM_DERIV_STEP, false);
		double flat22 = NumDerivative::getDerivative(FDiffM8<double>(X(0,0), lat, lon, getX, ReversedDirection), XT2, FirstDerivative, VariableX4, NUM_DERIV_STEP, false);
		double flon02 = NumDerivative::getDerivative(FDiffM8<double>(X(0,0), lat, lon, getX, ReversedDirection), XT2, FirstDerivative, VariableX5, NUM_DERIV_STEP, false);
		*/

		double xxx = 3.4466e06;


		//Eckert IV
		/*
		std::advance(i_proj, 16);
		//std::shared_ptr <Projection<double>> proj = *(projections.begin());

		//Eckert VI
		/*
		std::advance(i_proj, 18);

		//Hatano
		std::advance(i_proj,32);

		//MBTFPQ
		std::advance(i_proj, 44);

		//MBTFPS
		std::advance(i_proj, 45);



		//Mollweid
		//std::advance(i_proj, 47);

		//Wag IV
		//std::advance(i_proj, 72);

		//WINK II
		//std::advance(i_proj, 78);

		//Armadillo
		std::advance(i_proj, 6);

		std::shared_ptr <Projection<double>> proj2 = *i_proj;
		double x = proj2->getX(6380, 10, 10, -85, 15, 0, 0, 0, 0);
		double y = proj2->getY(6380, 10, 10, -85, 15, 0, 0, 0, 0);
		*/
	}

	//Some error occurs
	catch (Exception &e)
	{
		e.printException();
	}

	return 0;
}