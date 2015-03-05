// Description: Detection of the unknown cartographic projection using various methods

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

#include <iostream>
#include <ostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <time.h>

/*
//Memory leak analysis in VS 2010
#define  _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
//End memory leak analysis
*/

#include "libalgo/source/structures/point/Node3DCartesianProjected.h"
#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/structures/projection/Projection.h"

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/graticule/Graticule.h"
#include "libalgo/source/algorithms/ransac/Ransac.h"
#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/voronoi2D/Voronoi2D.h"
#include "libalgo/source/algorithms/faceoverlay/FaceOverlay.h"

#include "libalgo/source/algorithms/geneticalgorithms/DifferentialEvolution.h"

#include "libalgo/source/comparators/sortPointsByLat.h"

#include "libalgo/source/io/File.h"
#include "libalgo/source/io/DXFExport.h"

#include "libalgo/source/exceptions/Error.h"
#include "libalgo/source/exceptions/ErrorMathInvalidArgument.h"
#include "libalgo/source/exceptions/ErrorFileRead.h"

#include "libalgo/source/algorithms/transformation/HomotheticTransformation2D.h"
#include "libalgo/source/algorithms/turningfunction/TurningFunction.h"
#include "libalgo/source/algorithms/innerdistance/InnerDistance.h"
#include "libalgo/source/algorithms/projectiontoproj4/ProjectionToProj4.h"

//For tests
#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative2Var.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjJ.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjC.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV2S.h"
//#include "libalgo/source/algorithms/bisection/Bisection.h"

#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV2.h"

void printHeader ( std::ostream * output )
{
	//Print header
	*output << "  ************************************************************************" << std::endl;
	*output << "  *                                                                      *" << std::endl;
	*output << "  *               DETECTION OF THE CARTOGRAPHIC PROJECTION               *" << std::endl;
	*output << "  *                                                                      *" << std::endl;
	*output << "  *                (C)2010-2014 Tomas Bayer, v. 1.141012                 *" << std::endl;
	*output << "  *                       tomas.bayer@natur.cuni.cz                      *" << std::endl;
	*output << "  *                                                                      *" << std::endl;
	*output << "  ************************************************************************" << std::endl << std::endl << std::endl;
}


void printHelp ( std::ostream * output )
{
	//Print help
	*output << "\nParameters: \n\n";
	*output << "  test.txt          Set a test file 1 (unknown projection, cartesian coordinates).\n";
	*output << "  reference.txt     Set a reference file 2 (geographic coordinates).\n\n";
	*output << "  -h                Perform heuristic (throwing non-perspective samples before analysis_type).\n";
	*output << "  -n                Detect a projection in the normal aspect.\n";
	*output << "  -t                Detect a projection in the transverse aspect.\n";
	*output << "  -o                Detect a projection in the oblique aspect.\n";
	*output << "  -r                Remove outliers from both datasets, optimize transformation key. \n";
	*output << "  -c                Correct rotation of the sample (sample additionally rotated by an angle \n";
	*output << "                       of (87.5, 92.5)+0.5*k*M_PI deg).\n";
	*output << "  -?                Help. \n\n";
	*output << "Commands: \n\n";
	*output << "  +an=              Set a type of analysis: and, homt, helt, gntf, vdtf, all.\n";
	*output << "  +dlatp=           Set a value increment for the latitude of the cartographic \n";
	*output << "                       pole for both normal and transverse aspects (default dlatp = 10 deg).\n";
	*output << "  +dlonp=           Set a value increment for the longitude of the cartographic \n";
	*output << "                       pole for both normal and transverse aspects (default dlonp = 10 deg).\n";
	*output << "  +dlat=            Set the latitude increment of adjacent parallels in the graticule (deafault dlat = 10deg ).\n";
	*output << "  +dlon=            Set the longitude increment of adjacent meridians in the graticule (deafault dlat = 10deg ).\n";
	*output << "  +dlat0=           Set a value increment for the latitude of true \n";
	*output << "                       parallel (default dlat0 = 10 deg).\n";
	*output << "  +lon0=            Set a central meridian to lon0 (default lon0 = 0 deg) \n";
	*output << "  +proj=            Set a cartographic projection that will be analyzed (names are \n";
	*output << "                       in accordance with Proj.4, look to the projections.txt file).\n";
	*output << "  +latp=            Set a latitude of the the cartographic pole for analyzed projection.\n";
	*output << "  +lonp=            Set a longitude of the the cartographic pole for analyzed projection.\n";
	*output << "  +lat0=            Set an undistorted parallel latitude for analyzed projection.\n";
	*output << "  +met=             Set a method of analysis: gs, de, mls (default met=mls).\n";
	*output << "  +res=             Set number of printed results (default res=50).\n";
	*output << "  +gr=              The number exported graticules \n";
	*output << "  +sens=            Increase the heuristic sensitivity, more smaples will be added for \n";
	*output << "                       testing (default sens = 1).\n";
	*output << "  +dsens=           The heuristic sensitivity increment for each repetition, where sens \n";
	*output << "                       += dsens (default dsens = 2).\n";
	*output << "  +rep=             The number of repetition of the analysis_type with the increasing \n";
	*output << "                       sensitivity (default rep = 0).\n";
	*output << "  +match=           Set uncertainty regions for absolute matching: circle or Tissot Indicatrix\n";
	*output << "                       (match=circ, match=tiss). \n\n";
	*output << "  Examples: \n\n";
	*output << "  detectproj.exe +met=de +an=all vogt_test vogt_ref.txt \n";
	*output << "  detectproj.exe -n +met=de +an=all projections.txt vogt_test vogt_ref.txt \n";
	*output << "  detectproj.exe -n -t +an=cnd +an=homt projections.txt vogt_test vogt_ref.txt \n";
	*output << "  detectproj -o -h +met=mls +an=cnd +res=10 +sens=2 +dsens=2 +rep=3 +match=tiss points_test.txt \n";
	*output << "     points_ref.txt \n";
	*output << "  detectproj -o -h +met=de +an=all +dlatp=20 +dlonp=20 +dlat0=5 +res=10 +sens=2 +dsens=2 +rep=3 \n";
	*output << "     +match=circ points_test.txt points_ref.txt \n";
	*output << "  detectproj -o -h -r +met=de +an=all +res=10 +sens=2 +dsens=2 +rep=3 +proj=lcc +latp=50 +lonp=15 \n";
	*output << "     +lat0=50 +proj=eck5 +match=tiss points_test.txt points_ref.txt \n";
	*output << "  detectproj -o -h -r -c +met=de +an=all +res=10 +sens=2 +dsens=2 +rep=3 +proj=lcc +latp=50 \n";
	*output << "     +lonp=15 +lat0=50 +proj=eck5 +match=tiss points_test.txt points_ref.txt \n";
}

int main ( int argc, char * argv[] )
{
	//Initialize output
	std::ostream * output = &std::cout;
	static std::ofstream output_file;
	Container <Node3DCartesian <double> , NonDestructable > face_f;

	try
	{
		//Line arguments parameters
		TAnalysisParameters <double> analysis_parameters;
		TAnalyzedProjParametersList <double> ::Type analyzed_proj_parameters_list;
		TAnalyzedProjParameters <double> analyzed_proj_parameters;

		/*
		analysis_parameters.print_exceptions = false;
		analysis_parameters.analyze_normal_aspect = false;
		analysis_parameters.analyze_transverse_aspect = false;
		analysis_parameters.analyze_oblique_aspect = true;
		analysis_parameters.perform_heuristic = false;
		analysis_parameters.remove_outliers = false;
		analysis_parameters.correct_rotation = false;
		//analysis_parameters.analysis_method = SimplexRot2Method;
		analysis_parameters.analysis_method = SimplexRotMethod;
		//analysis_parameters.analysis_method = SimplexMethod;
		//analysis_parameters.analysis_method = SimplexShiftsMethod;
		//analysis_parameters.analysis_method = NonLinearLeastSquaresMethod;
		//analysis_parameters.analysis_method = NonLinearLeastSquaresRot2Method;
		//analysis_parameters.analysis_method = NonLinearLeastSquaresRotMethod;
		analysis_parameters.analysis_method = DifferentialEvolutionRotMethod;
		analysis_parameters.analysis_method = DifferentialEvolutionRot2Method;
		analysis_parameters.analysis_method = NonLinearLeastSquaresRot2Method;
		analysis_parameters.analysis_method = NonLinearLeastSquaresMethod;
		analysis_parameters.analysis_method = SimplexRot2Method;
		analysis_parameters.analysis_method = DifferentialEvolutionMethod;
		analysis_parameters.analysis_method = DifferentialEvolutionRot2Method;
		//analysis_parameters.analysis_method = NonLinearLeastSquaresMethod;
		analysis_parameters.analysis_method = NonLinearLeastSquaresRot2Method;
		//analysis_parameters.analysis_method = NonLinear;
		//analysis_parameters.analysis_method = SimplexRot2Method;
		analysis_parameters.projections_file = "projections_sinu2.txt";
		//analysis_parameters.analysis_method = DifferentialEvolutionRot2Method;

		analysis_parameters.analysis_type.a_cnd = false;
		analysis_parameters.analysis_type.a_gn_tf = false;
		analysis_parameters.analysis_type.a_helt = true;
		analysis_parameters.analysis_type.a_homt = true;
		analysis_parameters.analysis_type.a_vd_tf = false;
		analysis_parameters.analysis_type.a_vd_id = false;
		analysis_parameters.lat_step = 10;
		analysis_parameters.lon_step = 10;
		analysis_parameters.printed_results = 30;
		analysis_parameters.exported_graticule = 1;


		/*
		Container <Projection <double> *> projections1;
		projections1.load("projections_nicol.txt");
		Projection <double> *proj1 = projections1[0];
		Point3DGeographic <double> cart_pole1(90, 0);
		proj1->setCartPole(cart_pole1);
		proj1->setLat0(30.0);
		proj1->setR(1);

		//Point3DGeographic<double> p_temp(-89.997, -79.999);
		Point3DGeographic<double> p_temp(70, 10);

		const double x_temp = CartTransformation::latLonToX(&p_temp, proj1, analysis_parameters.print_exceptions);
		const double y_temp = CartTransformation::latLonToY(&p_temp, proj1, analysis_parameters.print_exceptions);

		
		/*
		unsigned short iterations = 0;
		const unsigned short max_iterations = 20;
		double  xmin, fmin, lon0;
		const double eps = 0.1, max_diff = 0.01;

		//Set limits for bisection
		Matrix <double> A(1, 1), B(1, 1), X(1,1);
		A(0, 0) = -180; B(0, 0) = 180;
	
		//Bisection::test(FAnalyzeProjV2<double>::functionResBis, A);
		*/
		/*
		double test = -180.3;
		if (test < MIN_LON)  test = MIN_LON - fmod(test, MIN_LON);

		double latp = 50, lonp = 0;
		double lat = 30, lon = -10;

		double lat_trans = CartTransformation::latToLatTrans(lat, lon, latp, lonp);
		double lont_trans = CartTransformation::lonToLonTrans(lat, lon, lat_trans, latp, lonp, ReversedDirection);
		/*
		// [latp =latp,  lonp = lonp + 180], Normaldirection2
		double latp2 = 50, lonp2 = lonp - 180;
		TTransformedLongtitudeDirection lon_direction = NormalDirection2;
		double lat_trans2 = CartTransformation::latToLatTrans(lat, lon, latp2, lonp2);
		double lont_trans2 = CartTransformation::lonToLonTrans(lat, lon, lat_trans2, latp2, lonp2, lon_direction);

		Container <Point3DGeographic<double> > points;
		points.load<Dim2D>("reference.txt");
		unsigned int n = points.size();

		for (unsigned int i = 0; i < n; i++)
		{
			double latp = 40, lonp = 0;
			double lat = 50;//points[i].getLat();
			double lon = 50;//points[i].getLon();

			double lat_trans = CartTransformation::latToLatTrans2(lat, lon, latp, lonp);
			double lont_trans = CartTransformation::lonToLonTrans2(lat, lon, lat_trans, latp, lonp);

			std::cout << lat_trans << "  " << lont_trans << '\n';

			// [latp =latp,  lonp = lonp + 180], Normaldirection2
			latp = 50, lonp = 180;
			TTransformedLongtitudeDirection lon_direction = NormalDirection2;
			double lat_trans2 = CartTransformation::latToLatTrans(lat, lon, latp, lonp);
			double lont_trans2 = CartTransformation::lonToLonTrans(lat, lon, lat_trans, latp, lonp, lon_direction);

			std::cout << lat_trans2 << "  " << lont_trans2 << "\n\n";
		}
		**/
		/*
		Container <Projection <double> *> projections;
		projections.load("projections_four2.txt");
		Projection <double> *proj = projections[0];


		//double xxx = CartTransformation::latLonToX(proj->getXEquat(), proj->getFThetaEquat(), proj->getTheta0Equat(), -80.0, 15.0, 1.0, 1.0, 1.0, 0.0, 1.0, 30.0, 0.0, 0.0, false);
		//double yyy = CartTransformation::latLonToX(proj->getYEquat(), proj->getFThetaEquat(), proj->getTheta0Equat(), -80.0, 15.0, 1.0, 1.0, 1.0, 0.0, 1.0, 30.0, 0.0, 0.0, false);


		double lat = 20, lon = 15, R = 1, RO = 180 / M_PI;;
		double y = 0.5*(sign(lat + atan(cos(0.5*lon) / tan(20/RO))*RO) - 1)*R*((1 + sin(20/RO) - cos(20)) / 2 + sin(-atan(cos(0.5*lon) / tan(20))*RO)*cos(20) - (1 + cos(-atan(cos(0.5*lon) / tan(20))*RO))*sin(20)*cos(lon / 2)) +
			0.5*(sign(lat + atan(cos(0.5*lon) / tan(20/RO))*RO) + 1)*R*((1 + sin(20/RO) - cos(20)) / 2 + sin(lat)*cos(20) - (1 + cos(lat))*sin(20)*cos(lon / 2));


		//0.0600282      42.2105868 - 0.0000298
		//Point3DGeographic <double> cart_pole(42.2105868, -0.0000298);
		Point3DGeographic <double> cart_pole(90, 0);
		proj->setCartPole(cart_pole);
		proj->setLat0(28);
		proj->setLat1(0.0);
		proj->setLat2(60.0);
		proj->setLon0(0);
		proj->setR(0.1);
		proj->setC(0.5);

		std::string proj4_string = ProjectionToProj4::ProjectionToProj4String(proj);

		//Create graticule
		analysis_parameters.lat_step = 5.0;
		analysis_parameters.lon_step = 5.0;

		TMeridiansListF <double> ::Type meridians_exp;
		TParallelsListF <double> ::Type parallels_exp;
		Container <Node3DCartesian <double> *> mer_par_points;

		const TMinMax <double> lat_interval(10, 40);
		const TMinMax <double> lon_interval(-30, 30);
		unsigned int index = 0;
		double alpha = 0;
		//const TMinMax <double> lat_interval ( -83.0, 77.0 );
		
		Graticule::computeGraticule(proj, lat_interval, lon_interval, analysis_parameters.lat_step, analysis_parameters.lon_step, 0.1 * analysis_parameters.lat_step, 0.1 * analysis_parameters.lon_step, alpha, TransformedGraticule, meridians_exp, parallels_exp, &mer_par_points, index);
		
		//Create file name
		char output_file_graticule[MAX_TEXT_LENGTH];
		strcpy(output_file_graticule, "graticule_");
		strcat(output_file_graticule, proj->getProjectionName());
		strcat(output_file_graticule, ".dxf");

		//Export graticule into DXF file
		//const double font_height = 0.05 * proj->getR() * std::min(analysis_parameters.lat_step, analysis_parameters.lon_step) * M_PI / 180;
		DXFExport::exportGraticuleToDXF(output_file_graticule, meridians_exp, parallels_exp, mer_par_points, 0.01, 10.0, 10.0);

	

		/*
		//Load prOjection
		Container <Projection <double> *> projections;
		projections.load("projections_hamm.txt");
		Projection <double> *proj = projections[0];

		//Load point
		Container < Point3DGeographic <double> *> pl_reference;
		pl_reference.load<Dim2D>("E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\AirRoutes\\reference2.txt");
		int nn = pl_reference.size();

		Point3DGeographic <double> cart_pole(45, 0);
		proj->setCartPole(cart_pole);
		proj->setLat0(0.0);
		proj->setLon0(0.0);
		proj->setR(0.07);

		for (unsigned int i = 0; i < nn; i++)
		{

				//Convert geographic point to oblique aspect
				const double lat_trans = CartTransformation::latToLatTrans((pl_reference)[i]->getLat(), (pl_reference)[i]->getLon(), proj->getCartPole().getLat(), proj->getCartPole().getLon());
				const double lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), (pl_reference)[i]->getLon(), lat_trans, proj->getCartPole().getLat(), proj->getCartPole().getLon(), ReversedDirection2);
				
				std::cout << lat_trans << "  " << lon_trans << '\n';

    				Point3DGeographic<double> p_temp(lat_trans, lon_trans);
				const double x_temp = CartTransformation::latLonToX(&p_temp, proj, analysis_parameters.print_exceptions) + 0.2415;
				const double y_temp = CartTransformation::latLonToY(&p_temp, proj, analysis_parameters.print_exceptions) - 0.1274;

				std::cout << x_temp  << "  " << y_temp << '\n';
				
		}


		Point3DGeographic<double> p1(45, 16), p(50, 15);

		double lat_transp1 = CartTransformation::latToLatTrans(p1.getLat(), p1.getLon(), p.getLat(), p.getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(p1.getLat(), p1.getLon(), lat_transp1, p.getLat(), p.getLon(), NormalDirection2);
		
		//analysis_parameters.lon0 = -70.3426043;
		/*
		TAnalyzedProjParametersList <double> ::Type analyzed_proj_parameters_list;
		TAnalyzedProjParameters <double> analyzed_proj_parameters;
		strcpy ( analyzed_proj_parameters.proj_name, "eqc" );
		analyzed_proj_parameters.latp = 0;
		analyzed_proj_parameters.lonp = 90;
		analyzed_proj_parameters.lat0 = 0;
		analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );

		/*
		TAnalysisParameters <double> analysis_parameters ( true);
		analysis_parameters.print_exceptions = false;
		analysis_parameters.remove_outliers = true;
		*/
		/*

		double nnn = -37;
		double kkk = nnn -  fmod (nnn, 10.0);

		//Get a sample proj sample
		Container <Projection <double> *> projections;
		projections.load ( "projections_api.txt" );
		Projection <double> *proj = projections[0];

		Point3DGeographic <double> cart_pole ( 90.0, 0.0 );
		proj->setCartPole ( cart_pole );
		proj->setLat0 ( 0.0 );
		proj->setLon0 ( 71.8 );

		//Get limits
		//const TMinMax <double> lon_interval ( -176.49, 168.13 );
		//const TMinMax <double> lon_interval ( -21.913862000000002, 120.98391300000000 );
		const TMinMax <double> lat_interval ( 50,70 );
		const TMinMax <double> lon_interval ( -20, 160);
		//const TMinMax <double> lat_interval ( -83.0, 77.0 );

		//const TMinMax <double> lon_interval ( 8.512969444, 8.626222222 );
		//const TMinMax <double> lat_interval ( 47.32697778, 47.38734444 );

		//Create data structures for graticule
		unsigned int index = 0;
		TMeridiansListF <double> ::Type meridians_exp;
		TParallelsListF <double> ::Type parallels_exp;
		Container <Node3DCartesian <double> *> mer_par_points;


		*/
		/*
		Matrix <double> A(4, 4);

		A(0, 0) = 3.8020;    A(0, 1) = 3.0397;    A(0, 2) = 2.1557;    A(0, 3) = 2.1239;
		A(1, 0) = 3.0397;    A(1, 1) = 3.8405;    A(1, 2) = 2.5496;    A(1, 3) = 2.4076;
		A(2, 0) = 2.1557;    A(2, 1) = 2.5496;    A(2, 2) = 3.2655;    A(2, 3) = 2.8074;
		A(3, 0) = 2.1239;    A(3, 1) = 2.4076;    A(3, 2) = 2.8074;    A(3, 3) = 3.8146;
		
		unsigned int mm = 3;
		//const Matrix <double> AC = ones(mm, mm, 1.0) * 10;
		Matrix <double> B = A - 10;
		B.print();
		double mmm = median(B, 1);
		double mmas = mad(B, 1);

		//Outliers detection
		Container < Node3DCartesian <double> *> PP, QQ, RR, SS;
		QQ.load<Dim2D>("E:\\tomas\\cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t25_26\\eqdc\\Terr1\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt");
		PP.load<Dim2D>("E:\\tomas\\cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t25_26\\eqdc\\Terr1\\t7_random_decr_prop_test_eqdc_proj#-80_-160_80_160_10.txt");

		//Convert to matrices
		unsigned int m = PP.size();

		Matrix <double> P(m, 2), Q(m, 2), W = eye(2 * m, 2 * m, 1.0), Eps(2 * m, 1);
		Matrix <unsigned int> I(2 * m, 1);

		for (int i = 0; i < m; i++)
		{
			P(i, 0) = PP[i]->getX(); P(i, 1) = PP[i]->getY();
			Q(i, 0) = QQ[i]->getX(); Q(i, 1) = QQ[i]->getY();
		}

		double eps_init = 0, eps = 0;
		unsigned int iter = 0;
		P = P * 1.0 / 100000;  Q = Q * 1.0 / 100000;
		TDevIndexPairs <double>::Type min_pairs;
		//Outliers::findOutliersME(Q, P, 2.0, 1.0e-10, SimilarityScheme, TukeyFunction, 20, W, I, Eps, eps_init, eps, iter);
		Outliers::findOutliersME(PP, QQ, RR, SS, 2.0, 1.0e-10, SimilarityScheme, DanishFunction2, 20, min_pairs, eps_init, eps, iter);
		*/
		
		//Matrix <double> A(7, 7);
		/*
		A(0, 0) = 8.9765;    A(0, 1) = 5.8414;    A(0, 2) = 8.1420;    A(0, 3) = 6.0109;    A(0, 4) = 6.7039;
		A(1, 0) = 5.8414;    A(1, 1) = 9.6482;    A(1, 2) = 8.3619;    A(1, 3) = 7.2455;    A(1, 4) = 5.3522;
		A(2, 0) = 8.1420;    A(2, 1) = 8.3619;    A(2, 2) = 14.1154;    A(2, 3) = 9.5547;    A(2, 4) = 8.4110;
		A(3, 0) = 6.0109;    A(3, 1) = 7.2455;    A(3, 2) = 9.5547;    A(3, 3) = 10.1250;    A(3, 4) = 6.2939;
		A(4, 0) = 6.7039;    A(4, 1) = 5.3522;    A(4, 2) = 8.4110;    A(4, 3) = 6.2939;    A(4, 4) = 8.2205;
		*/
		/*
		A(0, 0) = 72.2415224;	 A(0, 1) = -2258.9677161;     A(0, 2) = 861.3824291;	  A(0, 3) = -5657.5916356;    A(0, 4) = 0.0000000;       A(0, 5) = 0.0000000;  A(0, 6) = 0.0000007;
		A(1, 0) = -2258.9677161; A(1, 1) = 777774.4457143;    A(1, 2) = -120625.5584789;  A(1, 3) = 60657.8402929;    A(1, 4) = 0.0000000;       A(1, 5) = 0.0000000;  A(1, 6) = 240220.9126091;
		A(2, 0) = 861.3824291;	 A(2, 1) = -120625.5584789;   A(2, 2) = 571552.1902263;   A(2, 3) = 28304.6444900;    A(2, 4) = 0.0000000;       A(2, 5) = 0.0000000;  A(2, 6) = -403506.9759796;
		A(3, 0) = -5657.5916356; A(3, 1) = 60657.8402929;     A(3, 2) = 28304.6444900;    A(3, 3) = 639262.2044900;   A(3, 4) = 0.0000000;       A(3, 5) = 0.0000000;  A(3, 6) = -128585.5374394;
		A(4, 0) = 0.0000000;     A(4, 1) = 0.0000000;         A(4, 2) = 0.0000000;        A(4, 3) = 0.0000000;	      A(4, 4) = 0.0000000;       A(4, 5) = 0.0000000;  A(4, 6) = 0.0000000;
		A(5, 0) = 0.0000000;     A(5, 1) = 0.0000000;         A(5, 2) = 0.0000000;        A(5, 3) = 0.0000000;        A(5, 4) = 0.0000000;       A(5, 5) = 0.0000000;  A(5, 6) = 0.0000000;
		A(6, 0) = 0.0000007;	 A(6, 1) = 240220.9126091;    A(6, 2) = -403506.9759796;  A(6, 3) = -128585.5374394;  A(6, 4) = 0.0000000;       A(6, 5) = 0.0000000;  A(6, 6) = 855123.1402561;


		const unsigned int m = A.rows(), n = A.cols();

		Matrix <double> LL(m, 1, 1), RR(m, 1, 1);
		Matrix <double> T = tridiag(A, LL, RR);
		T.print();

		//Matrix <double> CH = chol(A, upper);
		//CH.print();

		Matrix <double> R(m, m), E(m, m);
		bool pos_def = true;
		gill(A, R, E, pos_def);

		Matrix <double> L2(m, m), d(m, 1), e(m, 1), neg(m, 1);
		gill(A, L2, d, e, neg);
		L2.print();
		d.print();
		e.print();
		neg.print();

		Matrix <double> BBB = L2*diag(d)*trans(L2);

		Matrix <double> V(m, m, 0, 1), L(m, 1);
		Matrix <unsigned int> IX(m,1);

		eig(A, V, L);
		V.print();
		L.print();
		sort(L, IX, 0);
		L.print();
		IX.print();
		//AT.print();
		
		/*

		//Create graticule
		analysis_parameters.lat_step = 10.0;
		analysis_parameters.lon_step = 10.0;

		Graticule::computeGraticule ( proj, lat_interval, lon_interval, analysis_parameters.lat_step, analysis_parameters.lon_step, 0.1 * analysis_parameters.lat_step, 0.1 * analysis_parameters.lon_step, TransformedGraticule, meridians_exp, parallels_exp, &mer_par_points, index );

		//Create file name
		char output_file_graticule[MAX_TEXT_LENGTH];
		strcpy ( output_file_graticule, "graticule_" );
		strcat ( output_file_graticule, proj->getProjectionName() );
		strcat ( output_file_graticule, ".dxf" );

		//Export graticule into DXF file
		DXFExport::exportGraticuleToDXF ( output_file_graticule, meridians_exp, parallels_exp, mer_par_points, 10.0, 10.0, 10.0 );

		/*
		double latp1 = 0;
		double lonp1 = 71.699;

		double lat_transp1 = CartTransformation::latToLatTrans(latp1, lonp1, cart_pole.getLat(), cart_pole.getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(latp1, lonp1, lat_transp1, cart_pole.getLat(), cart_pole.getLon(), NormalDirection );

		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';

		double xxp1 = CartTransformation::latLonToX ( projections[0]->getXEquat(), lat_transp1, lon_transp1, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );
		double yyp1 = CartTransformation::latLonToY ( projections[0]->getYEquat(), lat_transp1, lon_transp1, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );

		std::cout << xxp1 << "  " << yyp1 << '\n';
		*/

		/*
		Container <Projection <double> *> projections;
		projections.load ( "projections_nicol.txt" );
		Projection <double> *proj = projections[0];

		Point3DGeographic <double> cart_pole ( 90, 0 );
		proj->setCartPole ( cart_pole );
		proj->setLat0 ( 0.0 );
		proj->setLon0 ( 0.0 );

		double lat = -60, lon = -60;
		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';

		double xxp1 = CartTransformation::latLonToX ( projections[0]->getXEquat(), lat, lon, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );
		double yyp1 = CartTransformation::latLonToY ( projections[0]->getYEquat(), lat, lon, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );
		*/

		/*
		Sample <double> s1, s2;
		s1.setCrossNearestNeighbourDistanceRatio (1.03e5);
		s1.setHelmertTransformationRatio (1.55e5);
		s1.setHomotheticTransformationRatio (1.53e5);
		s1.setGNTurningFunctionRatio (1);
		s1.setVoronoiCellTurningFunctionRatio (3.73e0);
		s1.setCrossNearestNeighbourDistanceRatioPosition (22);
		s1.setHelmertTransformationRatioPosition (2);
		s1.setHomotheticTransformationRatioPosition (2);
		s1.setGNTurningFunctionRatioPosition (1);
		s1.setVoronoiCellTurningFunctionRatioPosition (27);

		s2.setCrossNearestNeighbourDistanceRatio (9.37e4);
		s2.setHelmertTransformationRatio (1.30e5);
		s2.setHomotheticTransformationRatio (1.30e5);
		s2.setGNTurningFunctionRatio (1);
		s2.setVoronoiCellTurningFunctionRatio (6.15e0);
		s2.setCrossNearestNeighbourDistanceRatioPosition (4);
		s2.setHelmertTransformationRatioPosition (1);
		s2.setHomotheticTransformationRatioPosition (1);
		s2.setGNTurningFunctionRatioPosition (1);
		s2.setVoronoiCellTurningFunctionRatioPosition (51);

		double c1 = s1.getSampleCost(analysis_parameters.analysis_type);
		double c2 = s2.getSampleCost(analysis_parameters.analysis_type);

		/*
		//Matrix <double> J ( 2, 2 );
		//J ( 0, 0 ) = 1; J ( 0, 1 ) = 2; J ( 1, 0 ) = 3; J ( 1, 1 ) = 4;
		Matrix <double> J ( 4, 4 );
		J ( 0, 0 ) = 1; J ( 0, 1 ) = 2;  J ( 0, 2 ) = 0; J ( 0, 3 ) = 0;
		J ( 1, 0 ) = 3; J ( 1, 1 ) = 4;  J ( 1, 2 ) = 0; J ( 1, 3 ) = 0;
		J ( 2, 0 ) = -11; J ( 2, 1 ) = 7;  J ( 2, 2 ) = 0; J ( 2, 3 ) = 0;
		J ( 3, 0 ) = 12; J ( 3, 1 ) = 17;  J ( 3, 2 ) = 0; J ( 3, 3 ) = 0;

		//Matrix <double> B = MatrixOperations::pinv ( A );

		/*
		//Differential evolution
		unsigned int population = 20;
		double eps = 0.01, max_iterations = 300,  F = 0.8, CR = 0.5;

		Matrix <double> J ( 1, 2 ), b ( 1, 2 ), arg_min ( 1, 2 );
		J ( 0, 0 ) = -500.0;
		J ( 0, 1 ) = -500.0;
		b ( 0, 0 ) = 500.0;
		b ( 0, 1 ) = 500.0;

		DifferentialEvolution::getMinimum ( FSchwefel<double>(), a, b, population, max_iterations, eps, F, CR, arg_min );

		double lll=-347;
		/*
		Container < Node3DCartesian <double> *> points, points_un;
		points.push_back( new Node3DCartesian <double> ( 555.8071, 109.3) );
		points.push_back( new Node3DCartesian <double> ( 555.0003, -0.0001) );
		points.push_back( new Node3DCartesian <double> ( 555.0001, 0.0001) );
		points.push_back( new Node3DCartesian <double> ( 554.9999, -0.0003) );
		points.push_back( new Node3DCartesian <double> ( 554.9997, -0.0004) );

		points.removeDuplicateElements ( points.begin(), points.end(),  sortPointsByX (), isEqualPointByPlanarCoordinates <Node3DCartesian <double> *> () );

		points.print();
		*/
		/*
		Sample <double> s;
		s.setCrossNearestNeighbourDistanceRatioPosition ( 1 );
		s.setHelmertTransformationRatioPosition ( 1 );
		s.setHelmertTransformationRatioPosition ( 1 );
		s.setGNTurningFunctionRatioPosition ( -1 );
		s.setVoronoiCellTurningFunctionRatioPosition ( 558 );

		int nn1 = 0;

		int ss = 30;
		bool uu = ss;
		double  tt = ss > 0 ? log ( ss * ( float ) ( ( bool ) ( ++nn1 ) ) ) : 0;

		float f = true;

		double s1 = analysis_parameters.analysis_type.a_cnd ? log ( s.getCrossNearestNeighbourDistanceRatioPosition() * ( s.getCrossNearestNeighbourDistanceRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0 ) ) : 0 ;
		double s2 = analysis_parameters.analysis_type.a_homt ? log ( s.getHomotheticTransformationRatioPosition() * ( s.getHomotheticTransformationRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0 ) ) : 0 ;
		double s3 = analysis_parameters.analysis_type.a_helt ? log ( s.getHelmertTransformationRatioPosition() * ( s.getHelmertTransformationRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0 ) ) : 0 ;
		double s4 = analysis_parameters.analysis_type.a_gn_tf ? log ( s.getGNTurningFunctionRatioPosition() * ( s.getGNTurningFunctionRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0 ) ) : 0 ;
		double s5 = analysis_parameters.analysis_type.a_vd_tf ? log ( s.getVoronoiCellTurningFunctionRatioPosition() * ( s.getVoronoiCellTurningFunctionRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0 ) ) : 0 ;

		int lh = s.getGNTurningFunctionRatioPosition() > 0 && ++ nn1 ? 1.0 : -1.0;

		int lll = 36;
		*/
		/*
		Container < Node3DCartesian <double> *> c1, c2, c3, c4;
		c1.load<Dim3D>("D:\\tomas\\cpp\\detectproj\\detectproj\\data\\pl1.txt");
		c2.load<Dim3D>("D:\\tomas\\cpp\\detectproj\\detectproj\\data\\pl2.txt");

		TTransformationKeyHelmert2D <double> kh;
		TTransformationKeyHomothetic2D <double> kho;
		TWeights <double> ::Type w;
		w.push_back(0.01);
		w.push_back(0.02);
		w.push_back(0.03);
		HelmertTransformation2D::transformPoints(&c1,&c2,&c3, kh);
		HomotheticTransformation2D::transformPoints(&c1,&c2,&c4, kho);
		TAccuracyCharacteristics <double> accu = Transformation2D::getAccuracyCharacteristics ( &c1, &c2, &c3, kh, w );

		/*
		//Check transformation
		Container < Node3DCartesian <double> *> points1, points2, points3, points4, points5, points6, points7, points8, points9;
		points1.load<Dim3D> ( "D:\\tomas\\cpp\\detectproj\\detectproj\\data\\outliers12.txt" );
		points2.load<Dim3D> ( "D:\\tomas\\cpp\\detectproj\\detectproj\\data\\outliers13.txt" );

		TTransformationKeyHelmert2D <double> key, key_w, key2, key3;
		TDevIndexPairs <double>::Type min_pairs, min_pairs2;
		Transformation2D::findOutliersIRLS ( &points1, &points2, &points3, &points4, key, min_pairs );
		Transformation2D::transformPoints ( &points3, &points4, &points5, key );
		TAccuracyCharacteristics <double> deviations = Transformation2D::getAccuracyCharacteristics ( &points3, &points4, &points5, key );

		points3.clear();
		Transformation2D::transformPoints ( &points1, &points2, &points3, key2 );
		TAccuracyCharacteristics <double> deviations2 = Transformation2D::getAccuracyCharacteristics ( &points1, &points2, &points3, key2 );
		int h = 777;

		TWeights <double> ::Type weights;
		Transformation2D::findOutliersIRLS ( &points1, &points2, &points7, &points8, key, min_pairs2 );
		Transformation2D::transformPoints ( &points7, &points8, &points9, key );
		TAccuracyCharacteristics <double> deviations3 = Transformation2D::getAccuracyCharacteristics ( &points7, &points8, &points9, key );
		int lllll = 333;
		lllll++;
		*/

		/*
		//Node3DCartesian <double> *n1 = new Node3DCartesian <double> (0.0,0.0);
		//Node3DCartesian <double> *n2 = new Node3DCartesian <double> (100.0,0.0);
		//Node3DCartesian <double> *n3 = new Node3DCartesian <double> (0.0,0);
		//Node3DCartesian <double> *n4 = new Node3DCartesian <double> (-10.0,0);
		//Node3DCartesian <double> *n1 = new Node3DCartesian <double> (-2548062.2382922294,   3530186.2840261245);
		//Node3DCartesian <double> *n2 = new Node3DCartesian <double> (-3516492.7637307635,   3618481.3385191644);
		//Node3DCartesian <double> *n3 = new Node3DCartesian <double> (-2548062.2382922294,   3530186.2840261245);
		//Node3DCartesian <double> *n4 = new Node3DCartesian <double> (-1506591.6135771051,   3435231.9190545371);
		/*
		Node3DCartesian <double> *n1 = new Node3DCartesian <double> (-3662804.3586703199,   1207639.0032263314 );
		Node3DCartesian <double> *n2 = new Node3DCartesian <double> (-2690803.6641360773,   1178171.2399238907 );
		Node3DCartesian <double> *n3 = new Node3DCartesian <double> (-2690803.6641360773,   1178171.2399238907 );
		Node3DCartesian <double> *n4 = new Node3DCartesian <double> (-1645493.6042309552,   1146480.9853837891 );
		*/
		/*
		Node3DCartesian <double> *n1 = new Node3DCartesian <double> ( -3662804.3586703199,   1207639.0032263314 );
		Node3DCartesian <double> *n2 = new Node3DCartesian <double> ( -2690803.6641360773,   1178171.2399238907 );
		Node3DCartesian <double> *n3 = new Node3DCartesian <double> ( -2548062.2382922294,   3530186.2840261245 );
		Node3DCartesian <double> *n4 = new Node3DCartesian <double> ( -2690803.6641360773,   1178171.2399238907 );
		*/

		/*
		Node3DCartesian <double> *n1 = new Node3DCartesian <double> (-1506591.6135771051,   3435231.9190545371);
		Node3DCartesian <double> *n2 = new Node3DCartesian <double> (-2548062.2382922294,   3530186.2840261245);
		Node3DCartesian <double> *n3 = new Node3DCartesian <double> (-2548062.2382922294,   3530186.2840261245);
		Node3DCartesian <double> *n4 = new Node3DCartesian <double> (-3516492.7637307635,   3618481.3385191644);
		*/
		/*
		double x_int, y_int;
		unsigned short pos = LineLinePosition::get2LineSegmentsPosition ( n1->getX(), n1->getY(), n2->getX(), n2->getY(), n3->getX(), n3->getY(), n4->getX(), n4->getY(), x_int, y_int );
		int kkk = 7689;
		*/
		/*
		double lat1= 10, lon1 = 10, lat2 = -10, lon2 = -10, latp = 0, lonp = 0;
		double lat_trans1 = CartTransformation::latToLatTrans(lat1, lon1, latp, lonp);
		double lon_trans1 = CartTransformation::lonToLonTrans(lat1, lon1, lat_trans1, latp, lonp);
		double lat_trans2 = CartTransformation::latToLatTrans(lat1, lon2, latp, lonp);
		double lon_trans2 = CartTransformation::lonToLonTrans(lat1, lon2, lat_trans2, latp, lonp);
		double lat_trans3 = CartTransformation::latToLatTrans(lat2, lon1, latp, lonp);
		double lon_trans3 = CartTransformation::lonToLonTrans(lat2, lon1, lat_trans3, latp, lonp);
		double lat_trans4 = CartTransformation::latToLatTrans(lat2, lon2, latp, lonp);
		double lon_trans4 = CartTransformation::lonToLonTrans(lat2, lon2, lat_trans4, latp, lonp);
		*/

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\vogt\\vogt_test5.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\vogt\\vogt_ref.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\vogt\\vogt_test.txt";
		//analysis_parameters.reference_file = "d:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\vogt\\vogt_ref.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\test\\Noise100\\Grid\\t5\\eck5\\t5_grid_mov_test_eck5#10_120_15_125.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\test\\Noise100\\Grid\\t5\\t5_grid_mov_ref#10_120_15_125.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\proj\\points_test_wer.txt";
		//analysis_parameters.reference_file = "d:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\proj\\points_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\france_test_12.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\france_reference_12.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\Europe\\M2\\NM\\europe_test_50.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\Europe\\M2\\NM\\europe_reference_50.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\europe_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\europe_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\france_test_25.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\france_reference_25.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Africa_Delisle\\africa_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Africa_Delisle\\africa_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_sbirka\\africa_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_sbirka\\africa_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\eqc_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\eqc_ref.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\mls_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\mls_ref.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\zurich\\zurich\\zurich_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\zurich\\zurich\\zurich_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\greece2\\greece_test2.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\greece2\\greece_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\europe\\europest\\europest_test2.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\europe\\europest\\europest_reference2.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\meyer\\meyer\\meyer_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\meyer\\meyer\\meyer_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\meyer\\meyer\\meyerup_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\meyer\\meyer\\meyerup_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_aeqd_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_aeqd_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_ortho_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_ortho_reference.txt";

		//           analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M2\\DE\\visualization\\hondius_world_east_test.txt";
		//           analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M2\\DE\\visualization\\hondius_world_east_reference.txt";

		//	analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_west\\M2\\BFGS\\visualization\\hondius_world_west_test.txt";
		//        analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_west\\M2\\BFGS\\visualization\\hondius_world_west_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\America_Suetter\\M2\\BFGS\\batch\\america_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\America_Suetter\\M2\\BFGS\\batch\\america_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\GreatBrittain\\M2\\DE\\batch\\great_brittain_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\GreatBrittain\\M2\\DE\\batch\\great_brittain_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\cassini_soldner_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\cassini_soldner_reference3.txt";
		/*
		analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\Africa\\M2\\NM\\africa_test.txt";
		analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\africa\\M2\\NM\\africa_reference.txt";

		analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M2\\NM\\hondius_world_east_test.txt";
		analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M2\\NM\\hondius_world_east_reference.txt";
		*/

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk6\\m2\\nm\\batch\\jtsk_test6.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk6\\m2\\nm\\batch\\jtsk_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk\\m2\\de\\jtsk_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk\\m2\\de\\jtsk_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\africa_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\africa_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\EastHemisphereGrid\\M2\\NM\\east_hemisphere_test_grid.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\EastHemisphereGrid\\M2\\NM\\east_hemisphere_reference_grid.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\WestHemisphere\\M2\\BFGS\\west_hemisphere_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\WestHemisphere\\M2\\BFGS\\west_hemisphere_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\europe_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\europe_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\danemark_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\danemark_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk_test6.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk_reference.txt";
		//analysis_parameters.test_file = "jtsk_test8.txt";
		//analysis_parameters.reference_file = "jtsk_reference.txt";

		//analysis_parameters.test_file = "jtsk_test6.txt";
		//analysis_parameters.reference_file = "jtsk_reference.txt";

		//analysis_parameters.test_file = "jtsk_test.txt";
		//analysis_parameters.reference_file = "jtsk_reference.txt";

		//analysis_parameters.test_file = "s42_test.txt";
		//analysis_parameters.reference_file = "s42_reference.txt";

		//analysis_parameters.test_file = "africa_test.txt";
		//analysis_parameters.reference_file = "africa_reference.txt";

		//analysis_parameters.test_file = "hondius_world_west_test.txt";
		//analysis_parameters.reference_file = "hondius_world_west_reference.txt";

		//analysis_parameters.test_file = "hondius_world_east_test.txt";
		//analysis_parameters.reference_file = "hondius_world_east_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\s42_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\s42_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\s42_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\s42_reference.txt";

		//analysis_parameters.test_file = "//home//tomas//cpp//detectproj//detectproj//tests//maps//other//jtsk_test.txt";
		//analysis_parameters.reference_file = "//home//tomas//cpp//detectproj//detectproj//tests//maps//other//jtsk_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\time_complexity\\merc_20_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\time_complexity\\merc_20_ref.txt";
		/*
		int m = 6;
		Matrix <double> V(m, 1), VX(0.5 * m, 1), VY(0.5 * m, 1);
		for (int i = 0; i < 6; i++) V(i, 0) = i;
		VX = V(0, 0.5 * m - 1, 0, 0); VY = V(0.5 * m, m - 1, 0, 0);
		Matrix <double> VXY1 = VX % VX;
		Matrix <double> VXY2 = VY % VY;
		VXY1.print();
		VXY2.print();
		Matrix <double> VXY = VX % VX + VY % VY;
		VXY.print();

		/*
		Projection <double> * proj = new ProjectionConicLimits<double> ();
		Container < Point3DGeographic <double> *> points;
		double latp_min, latp_max, lonp_min, lonp_max;

		points.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\i34.txt" );
		CartAnalysis::findLatPLonPIntervals ( &points, proj, latp_min, latp_max, lonp_min, lonp_max );
		double latp = 90;
		double lonp = 0.0;//

		lonp = lonp_min;
		double lonp0 = 0;

		//Set pointer to correct value: do not compute different lonp values for normal aspect
		double * p_lonp = ( latp == MAX_LAT ? &lonp0 : &lonp );

		for ( ; ( latp == MAX_LAT ? *p_lonp == 0.0 : *p_lonp <= ( lonp_min < lonp_max ? lonp_max : 180.0 ) );
		*p_lonp = lonp + analysis_parameters.lonp_step )
		{
		std::cout << *p_lonp << '\n';
		}

		//Internal cycle: only when latp_min > latp_max
		for ( lonp = ( lonp_min < lonp_max ? lonp : -180.0 ) ; ( latp == MAX_LAT || lonp_min <  lonp_max ) ? false : lonp <= lonp_max ;
		lonp =  lonp + analysis_parameters.lonp_step )
		{
		std::cout  << lonp << '\n';
		}

		delete proj;
		*/
		/*
		Container < Point3DGeographic <double> *> points1;
		Container < Point3DGeographic <double> *> *points2 = NULL;
		Container < Point3DGeographic <double> *, NonDestructable> points3;
		Container < Point3DGeographic <double> *, NonDestructable> *points4 = NULL;
		Container < Node3DCartesian <double> *, NonDestructable> *nodes2;
		points1.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\cassini_soldner_reference.txt" );
		points3.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\cassini_soldner_test.txt" );

		points3 = &points1;
		points3[0]->print();
		&points4 = points1;
		*points2 = points4;
		*/
		/*
		//Least squares
		Container <Projection <double> *> projections;
		projections.load ( "projections_eqc.txt" );
		Projection <double> *proj = projections[0];

		//Load reference points
		Container < Point3DGeographic <double> *> points_r;
		//points_r.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\mls_ref.txt" );
		points_r.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_aeqd_reference.txt" );

		//Load test points
		Container < Node3DCartesian <double> *> points_t;
		//points_t.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\mls_test.txt" );
		points_t.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\azimuthal_aeqd_test.txt" );

		//Create matrices
		unsigned int m = points_t.size();

		//Apply Minimum least squares algorithm
		Matrix <double> X ( 3, 1 ), W ( 2 * m, 2 * m, 0.0, 1 ), L ( 2 * m, 1 );
		TProjectionAspect aspect = ObliqueAspect;

		//Set initial values
		//X ( 0, 0 ) = 0.0; X ( 1, 0 ) = -110.0; X ( 2, 0 ) = 25.0;
		X ( 0, 0 ) = 90.0; X ( 1, 0 ) = 0.0; X ( 2, 0 ) = 20.0;
		double lat = 10.0, lon = 10.0;

		double y = ArithmeticParser::parseEq ( proj->getYEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() );

		Matrix <double> args ( 1, 5 ); args ( 0, 0 ) = 6380000; args ( 0, 1 ) = 90; args ( 0, 2 ) = 0; args ( 0, 3 ) = 20; args ( 0, 4 ) = 0;

		double dR = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getXEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX4, NUM_DERIV_STEP, false );

		double dx_dlatp = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getXEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX4, NUM_DERIV_STEP, false );
		double dx_dlonp = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getXEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX5, NUM_DERIV_STEP, false );
		double dy_dlatp = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getYEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX4, NUM_DERIV_STEP, false );
		double dy_dlonp = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getYEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX5, NUM_DERIV_STEP, false );
		double dx_dlat0 = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getXEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX1, NUM_DERIV_STEP, false );
		double dy_dlat0 = 180.0 / M_PI * NumDerivative::getDerivative ( FProjEquationDerivative6Var <double> ( proj->getYEquat(), lat, lon, proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getDy(), proj->getLat1(), proj->getLat2() ), args, VariableX1, NUM_DERIV_STEP, false );

		double fklfg = 37;
		*/
		/*
		Matrix <double> J ( 4, 3 ), L ( 3, 3 ), U ( 4, 4 ), B1 ( 4, 3 ), V1 ( 3, 3 ),  P2 ( 3 , 3 ), Q ( 3, 3 ), R ( 3, 3 ), B ( 3, 1 ), C ( 1 , 3 );
		Matrix <unsigned int> P ( 3, 3 );

		J ( 0, 0 ) = 1.0; J ( 0, 1 ) = 2.0; J ( 0, 2 ) = 3.0;
		J ( 1, 0 ) = 4.0; J ( 1, 1 ) = 5.0; J ( 1, 2 ) = 6.0;
		J ( 2, 0 ) = 1.0; J ( 2, 1 ) = 0.0; J ( 2, 2 ) = 1.0;
		J ( 3, 0 ) = 0.0; J ( 3, 1 ) = 1.0; J ( 3, 2 ) = -1.0;

		//MatrixOperations::bidiag ( A, U, B1, V1 );

		Matrix <double> LL ( 7, 7 ), W (7, 7, 0.0, 1.0 );
		for ( int i = 0; i < 7; i++ ) LL ( i, i ) = 100 * i;
		*/
		/*
		Matrix <double> AA = MatrixOperations::load <double> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\AA.txt" );
		const unsigned int m = AA.rows(), n = AA.cols();
		Matrix <double> UU (m, m ), BB ( m, n), VV ( n, n );

		MatrixOperations::bidiag ( AA, UU, BB, VV );
		/*
		Matrix <double> AA = MatrixOperations::load <double> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\A.txt" );
		Matrix <double> LL = MatrixOperations::load <double> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\LA.txt" );
		Matrix <double> W = MatrixOperations::load <double> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\WA.txt" );
		//Matrix <double> AA = trans (AAA);

		//AA = trans ( AA ) * AA;
		const unsigned int m = AA.rows(), n = AA.cols();
		Matrix <double> UU (m, m ), BB ( m, n), VV ( n, n );

		//MatrixOperations::svd ( AA, UU, BB, VV, 500 );
		MatrixOperations::bidiag ( AA, UU, BB, VV );

		Matrix <double> Back = UU * BB* trans (VV);

		Matrix <double> AA_inv = MatrixOperations::pinvs ( AA );
		Matrix <double> AA2_inv = MatrixOperations::pinv ( AA );
		Matrix <double> AA3_inv = MatrixOperations::pinv1 ( AA );

		Matrix <double> dx = MatrixOperations::pinvs ( MatrixOperations::trans ( AA ) * W * AA ) * MatrixOperations::trans ( AA ) * W * LL;
		Matrix <double> dx2 = MatrixOperations::mlsqr(AA, W, LL,);
		Matrix <double> CC1 = AA * AA_inv * AA;
		//Matrix <double> CC2 = AA * AA2_inv * AA;

		Matrix <double> DIF = CC1 - AA;
		//Matrix <double> DIF2 = CC2 - AA;

		AA_inv.print();
		//AA2_inv.print();
		//DIF.print();
		//DIF2.print();

		//double res = sum (DIF -DIF2);

		double l = 337;
		/*
		Matrix <double> AJ ( 2, 2 ), UU ( 2, 2 ), BB ( 2, 2 ), VV ( 2, 2 );
		AJ ( 0, 0 ) = 1.0; AJ ( 0, 1 ) = 2.0;
		AJ ( 1, 0 ) = 3.0; AJ ( 1, 1 ) = 4.0;
		MatrixOperations::bidiag ( AA, UU, BB, VV );
		*/

		//Matrix <double> A_inv = MatrixOperations::pinvs ( A );
		//A_inv.print();

		//J ( 1, 0 ) = 0.0; J ( 1, 1 ) = 0.0; J ( 1, 2 ) = 0.0;
		//J ( 2, 0 ) = 0.0; J ( 2, 1 ) = 0.0; J ( 2, 2 ) = 0.0;
		//B (2, 0 ) = 5;
		//Matrix <double> AI = MatrixOperations::pinv ( A );
		/*
		Matrix <double> A1 ( 4, 4 ), L1 ( 4, 4 ), U1 ( 4, 4 ), P1 ( 4 , 4 );

		A1 ( 0, 0 ) = 1.0; A1 ( 0, 1 ) = 2.0; A1 ( 0, 2 ) = 3.0; A1 ( 0, 3 ) = 4.0;
		A1 ( 1, 0 ) = 4.0; A1 ( 1, 1 ) = 5.0; A1 ( 1, 2 ) = 6.0; A1 ( 1, 3 ) = 19.0;
		A1 ( 2, 0 ) = 1.0; A1 ( 2, 1 ) = 0.0; A1 ( 2, 2 ) = 1.0; A1 ( 2, 3 ) = -6.0;
		A1 ( 3, 0 ) = 1.0; A1 ( 3, 1 ) = 8.0; A1 ( 3, 2 ) = 8.0; A1 ( 3, 3 ) = -7.0;

		Matrix <double> A1_inv = MatrixOperations::pinvs ( A1 );
		A1_inv.print();

		A1 ( 0, 0 ) = -21.494185260204677; A1 ( 0, 1 ) = -2.9775494732751060;	   A1 ( 0, 2 ) = -3.1171221048348770;	A1 ( 0, 3 ) = -2.1866378944364060;
		A1 ( 1, 0 ) = -3.0565012460326580e-016;	A1 ( 1, 1 ) = -10.056550061238653;  A1 ( 1, 2 ) = -2.6568370421181817;	A1 ( 1, 3 ) = -9.2963468505483213;
		A1 ( 2, 0 ) = 1.4948707788404388e-015;	A1 ( 2, 1 ) = 0.00000000000000000; A1 ( 2, 2 ) = -1.4915651897180007;	A1 ( 2, 3 ) =  1.0156418206697397;
		A1 ( 3, 0 ) = 3.3290501264859602e-015;	A1 ( 3, 1 ) = 0.00000000000000000; A1 ( 3, 2 ) = 0.00000000000000000;	A1 ( 3, 3 ) = -0.87465515777456271;

		MatrixOperations::lu ( A1, L1, U1, P1 );
		Matrix <double> AII = MatrixOperations::pinv ( A1 );

		Matrix <double> A2 ( 8, 4 ), L2 ( 8, 1 );
		A2 ( 0, 0 ) = -46.706595565192401;   A2 ( 0, 1 ) = 5201303.6890143296 ;  A2 ( 0, 2 ) = -6935071.5601962805;   A2 ( 0, 3 ) = 1689866.1668453738;
		A2 ( 1, 0 ) = -57.756419708020985;   A2 ( 1, 1 ) = 7739727.0237239264 ;  A2 ( 1, 2 ) = -6791805.2074869163;   A2 ( 1, 3 ) = 2089653.8166015269   ;
		A2 ( 2, 0 ) = -62.461120393127203;   A2 ( 2, 1 ) = 14668719.092789855 ;  A2 ( 2, 2 ) = -8777768.1241737865;   A2 ( 2, 3 ) = 2259872.0143301692   ;
		A2 ( 3, 0 ) = -90.393512835726142;   A2 ( 3, 1 ) = 19409418.300810147 ;  A2 ( 3, 2 ) = -1552149.1542892903;   A2 ( 3, 3 ) = 3270479.0484867431   ;
		A2 ( 4, 0 ) = 48.593955868855119 ;   A2 ( 4, 1 ) = -4176998.0210933369;  A2 ( 4, 2 ) = -4176998.0210933369;   A2 ( 4, 3 ) = 0.00000000000000000   ;
		A2 ( 5, 0 ) = 55.636196834966540 ;   A2 ( 5, 1 ) = -3120137.5439594686;  A2 ( 5, 2 ) = -4819891.3836225914;   A2 ( 5, 3 ) = 0.00000000000000000   ;
		A2 ( 6, 0 ) = 69.351466507650912 ;   A2 ( 6, 1 ) = -2628631.2487878883;  A2 ( 6, 2 ) = -5034929.8781290278;   A2 ( 6, 3 ) = 0.00000000000000000   ;
		A2 ( 7, 0 ) = 72.707025418058038 ;   A2 ( 7, 1 ) = 560350.87873502634 ;  A2 ( 7, 2 ) = -5504298.5704466701;   A2 ( 7, 3 ) = 0.00000000000000000   ;

		L2 ( 0, 0 ) = 4686655.1327758105;
		L2 ( 1, 0 ) = 5795421.5792721715;
		L2 ( 2, 0 ) = 6267502.7400934119;
		L2 ( 3, 0 ) = 9070308.4931811783;
		L2 ( 4, 0 ) = 4869572.4639018495;
		L2 ( 5, 0 ) = 5575271.4803221179;
		L2 ( 6, 0 ) = 6949670.9729914768;
		L2 ( 7, 0 ) = 7285929.1907577384;

		Matrix<double> A22 = MatrixOperations::trans ( A2 ) * A2;

		//MatrixOperations::qr ( A, Q, R, P );
		Matrix <double> AI2 = MatrixOperations::pinv ( A22 );
		Matrix <double> AI3 = MatrixOperations::pinvs ( A22 );
		Matrix <double> DF = AI2 - AI3;
		Matrix <double> RR = MatrixOperations::trans ( A2 ) * L2;
		Matrix <double> XX = AI2 * RR;
		double nnn = 35;

		/*
		//Load cartographic projections from file
		Container <Projection <double> *> projections;
		projections.load ( "projections.txt" );
		Matrix <double> args(1,2);
		args(0,0) = 50.0; args(0,1) = 15.0;
		Projection <double> *proj = projections[0];

		double derx = NumDerivative::getDerivative(FProjEquationDerivative2Var <double> (proj->getXEquat(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
		proj->getDy(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0()) , args, VariableX2, 0.01, false );

		//double derx2 = NumDerivative::lonDerivative(0.01, proj->getXEquat(), 50.0, 15.0, proj->getR(), proj->getA(), proj->getB(), proj->getDx(),
		//	proj->getDy(), proj->getLat0(), proj->getLat1(), proj->getLat2(), proj->getLon0(), false );

		double kkk=37;

		*/
		/*
		//Load cartographic projections from file
		Container <Projection <double> *> projections;
		projections.load ( "projections_hamm.txt" );
		Container < Point3DGeographic <double> *> points;

		//points.load<Dim2D> ( "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\greece\\greece_test.txt" );
		//points.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\greece\\greece_reference.txt" );

		for (int i = 0; i < points.size(); i++)
		{
			double lat_transp1 = CartTransformation::latToLatTrans(points[i]->getLat(), points[i]->getLon(), 0.0, 90.0);
			double lon_transp1 = CartTransformation::lonToLonTrans(points[i]->getLat(), points[i]->getLon(), lat_transp1, 0.0, 90.0, NormalDirection );

			//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';

			double xxp1 = CartTransformation::latLonToX ( projections[0]->getXEquat(), lat_transp1, lon_transp1, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );
			double yyp1 = CartTransformation::latLonToY ( projections[0]->getYEquat(), lat_transp1, lon_transp1, projections[0]->getR(), 0.0, 0.0, 0.0, projections[0]->getLat0(), 0.0, 0.0, 0.0 );

			std::cout << xxp1 << "  " << yyp1 << '\n';
		}

		double kkk=37;
		kkk++;
		*/
		/*
		Container < Point3DCartesian <double> *> ppoints1, ppoints2, ppoints3;
		ppoints1.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk_test6.txt" );

		Container <Projection <double> *> proj_list2;
		proj_list2.load ( "projections_lccjtsk.txt" );
		Container < Point3DGeographic <double> *> ppoints;
		ppoints.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk_reference.txt" );

		double lat0=78.5, lon0 = -17.6666;
		double lat1 = 50.09044, lon1 = 	14.39908 - lon0;
		double latp = (59+42.0/60+42.6969/3600), lonp = 42+31.0/60+31.41725/3600, R=6380.703;

		/*
		double lat0=78.5887, lon0 = -17.6777;
		double latp = 59.091, lonp = 42.548;
		double dx = -1.0034, dy = -1.927, R = 6.3863;
		*/
		//double lat0=80.0709, lon0 = 0;
		//double latp = 49.391, lonp = 25.3351;
		//double dx = 36.603, dy = -1170.551, R = 6390.289;
		/*
		Point3DGeographic <double> pole (latp, lonp );
		proj_list2[0]->setR(R);
		proj_list2[0]->setCartPole(pole);
		proj_list2[0]->setLat0(lat0);
		proj_list2[0]->setLon0(lon0);
		*/
		/*

		Container <Projection <double> *> proj_list2;
		proj_list2.load ( "projections_lccjtsk.txt" );
		Container < Point3DGeographic <double> *> ppoints;
		ppoints.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk_reference.txt" );

		double lat0=78.5, lon0 = -17.6666;
		double lat1 = 50.09044, lon1 = 	14.39908 - lon0, dx =0, dy = 0;
		double latp = (59+42.0/60+42.6969/3600), lonp = 42+31.0/60+31.41725/3600, R=6380.703;

		Point3DGeographic <double> pole (latp, lonp );
		proj_list2[0]->setR(R);
		proj_list2[0]->setCartPole(pole);
		proj_list2[0]->setLat0(lat0);
		proj_list2[0]->setLon0(lon0);

		double res = 0;
		for (int i = 0; i < ppoints.size(); i++)
		{
		TTransformedLongtitudeDirection trans_lon_dir = proj_list2[0]-> getLonDir();
		double lat_transp1 = CartTransformation::latToLatTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, lat_transp1, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon(), trans_lon_dir );
		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';

		double R_check =  proj_list2[0]->getR();
		double lat0_check =  proj_list2[0]->getLat0();

		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';
		double xxp1 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 ) + dx;
		double yyp1 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 ) + dy;

		//Point3DCartesian <double> *p = new Point3DCartesian <double> (xxp1, yyp1);
		//ppoints2.push_back (p);

		std::cout << lat_transp1 << "  " << lon_transp1 << '\n';
		std::cout << xxp1 << "  " << yyp1 << '\n';

		//res += ( xxp1 - ppoints1[i]->getX() ) * ( xxp1 - ppoints1[i]->getX() ) + ( yyp1 - ppoints1[i]->getY() ) * ( yyp1 - ppoints1[i]->getY() );
		}*/

		/*
		lat0=83.9910, lon0 = -17.6666;
		latp = 50.8799, lonp = 42.8266;
		dx = 24.6760, dy = -1003.704, R = 6385.869;

		Point3DGeographic <double> pole2(latp, lonp );
		proj_list2[0]->setR(R);
		proj_list2[0]->setCartPole(pole2);
		proj_list2[0]->setLat0(lat0);
		proj_list2[0]->setLon0(lon0);

		res = 0;

		for (int i = 0; i < ppoints.size(); i++)
		{
		TTransformedLongtitudeDirection trans_lon_dir = proj_list2[0]-> getLonDir();
		double lat_transp1 = CartTransformation::latToLatTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, lat_transp1, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon(), trans_lon_dir );
		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';

		double R_check =  proj_list2[0]->getR();
		double lat0_check =  proj_list2[0]->getLat0();

		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';
		double xxp1 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 ) + dx;
		double yyp1 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 ) + dy;

		Point3DCartesian <double> *p = new Point3DCartesian <double> (xxp1, yyp1);
		ppoints3.push_back (p);

		std::cout << xxp1 << "  " << yyp1 << '\n';

		res += ( xxp1 - ppoints1[i]->getX() ) * ( xxp1 - ppoints1[i]->getX() ) + ( yyp1 - ppoints1[i]->getY() ) * ( yyp1 - ppoints1[i]->getY() );
		}

		double bear1 = Bearing::getBearing (ppoints2[0], ppoints2[1]);
		double bear2 = Bearing::getBearing (ppoints3[0], ppoints3[1]);

		std::cout << res;
		int lll = 29;

		*/
		/*
		Container <Projection <double> *> proj_list2;
		proj_list2.load ( "projections_lccjtsk.txt" );
		Container < Point3DGeographic <double> *> ppoints;
		Container < Point3DCartesian <double> *> ppoints1;
		ppoints.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk_reference_jtsk.txt" );

		double lat0=78.5, lon0 = -17.6666;
		double lat1 = 50.09044, lon1 = 	14.39908 - lon0;
		double latp = 59.7, lonp = 42.5;

		for (int i = 0; i < ppoints.size(); i++)
		{
		double lat_transp1 = CartTransformation::latToLatTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, lat_transp1, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon(), NormalDirection );
		//std::cout << lat_transp1 << "  " << lon_transp1 << '\n';
		double xxp1 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 );
		double yyp1 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 );
		std::cout << xxp1 << "  " << yyp1 << '\n';
		}
		int lll = 29;
		/*
		//Load cartographic projections from file
		Container <Projection <double> *> proj_list2;
		proj_list2.load ( "projections.txt" );
		Container < Point3DGeographic <double> *> ppoints;
		Container < Point3DCartesian <double> *> ppoints1, ppoints2, ppoints3;
		ppoints.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk_reference_jtsk.txt" );
		ppoints1.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk_test.txt" );
		ppoints2.load<Dim2D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk_lcc.txt" );

		double lat0=78.5, lon0 = -17.6666;

		//double lat1 = 50.05037516, lon1 = 32.08599652;
		double lat1 = 50.09044, lon1 = 	14.39908 - lon0;

		//double latp = (59+42.0/60+42.6969/3600), lonp = 42+31.0/60+31.41725/3600 ;
		double latp = 59.7, lonp = 42.5 ;

		//Compute Helmert transformation
		TTransformationKeyHelmert2D <double> key_helmert;
		HelmertTransformation2D::transformPoints ( &ppoints2, &ppoints1, &ppoints3, key_helmert );

		//Compute ratio and percentage match using uncertainty regions
		TIndexList matched_points;
		TAccuracyCharacteristics <double> deviations = Transformation2D::getAccuracyCharacteristics ( &ppoints2, &ppoints1, &ppoints3, key_helmert );
		double helmert_transformation_ratio = deviations.std_dev;
		//double helmert_transformation_perc_match = Transformation2D::getMatchRatioTissotIndicatrix ( &ppoints2, &ppoints3, matched_points, CollectOn, 0.5 );

		for (int i = 0; i < ppoints.size(); i++)
		{
		double lat_transp1 = CartTransformation::latToLatTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon());
		double lon_transp1 = CartTransformation::lonToLonTrans(ppoints[i]->getLat(), ppoints[i]->getLon() -lon0, lat_transp1, proj_list2[0]->getCartPole().getLat(), proj_list2[0]->getCartPole().getLon(), NormalDirection );
		std::cout << lat_transp1 << "  " << lon_transp1 << '\n';
		double xxp1 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 );
		double yyp1 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_transp1, lon_transp1, proj_list2[0]->getR(), 0.0, 0.0, 0.0, proj_list2[0]->getLat0(), 0.0, 0.0, 0.0 );
		std::cout << xxp1 << "  " << yyp1 << '\n';
		}

		double lat_trans1 = CartTransformation::latToLatTrans(lat1, lon1, latp, lonp);
		double lon_trans1 = CartTransformation::lonToLonTrans(lat1, lon1, lat_trans1, latp, lonp, NormalDirection );

		double xx1 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_trans1, lon_trans1, 6378.0, 0.0, 0.0, 0.0, lat0, 0.0, 0.0, lon0 );
		double yy1 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_trans1, lon_trans1, 6378.0, 0.0, 0.0, 0.0, lat0, 0.0, 0.0, lon0 );

		//lon0 = 10;
		//lon = 30;
		double lat2 = -10, lon2 = -180;

		double lat_trans2 = CartTransformation::latToLatTrans(lat2, lon2, latp, lonp);
		double  lon_trans2 = CartTransformation::lonToLonTrans(lat2, lon2, lat_trans2, latp, lonp, NormalDirection );

		double xx2 = CartTransformation::latLonToX ( proj_list2[0]->getXEquat(), lat_trans2, lon_trans2, 6378.0, 0.0, 0.0, 0.0, 45.0, 0.0, 0.0, lon0 );
		double yy2 = CartTransformation::latLonToY ( proj_list2[0]->getYEquat(), lat_trans2, lon_trans2, 6378.0, 0.0, 0.0, 0.0, 45.0, 0.0, 0.0, lon0 );

		double kkk=37;
		*/
		/*
		//Load cartographic projections from file
		Container <Projection <double> *> proj_list2;
		proj_list2.load ( "projections.txt" );
		double latp = 90, lonp = 0, lat = -80, lon = 0;
		double xx = CartTransformation::latLonToX ( proj_list2[19]->getXEquat(), lat, lon, 6378.0, 0.0, 0.0, 0.0, 45.0, 0.0, 0.0, 0.0 );
		double yy = CartTransformation::latLonToY ( proj_list2[19]->getYEquat(), lat, lon, 6378.0, 0.0, 0.0, 0.0, 45.0, 0.0, 0.0, 0.0 );
		double sofsd = 37;

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\test2\\Noise100\\Grid\\t2\\merc\\t2_grid_decr_prop_test_merc#-30_0_30_60.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\test2\\Noise100\\Grid\\t2\\t2_grid_decr_prop_ref#-30_0_30_60.txt";
		*/
		/*
		Container < Node3DCartesian <double> *> ppoints1, ppoints2, ppoints3, ppoints4, ppoints5;
		ppoints1.load<Dim3D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\turf1.txt" );
		ppoints2.load<Dim3D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\turf2.txt" );
		ppoints3.load<Dim3D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\f9_inner.txt" );
		ppoints4.load<Dim3D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\f4_inner.txt" );
		ppoints5.load<Dim3D> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\f5_inner.txt" );
		Container <Face <double> *> faces;
		for ( unsigned int i = 0; i < ppoints1.size(); i++ )
		{
		Container <HalfEdge<double> *> hl1, hl2, hl3, hl4, hl5;
		Face <double> * f1 = new Face <double> ( &ppoints1, &hl1 );
		Face <double> * f2 = new Face <double> ( &ppoints2, &hl2 );
		Face <double> * f3 = new Face <double> ( &ppoints3, &hl3 );
		Face <double> * f4 = new Face <double> ( &ppoints4, &hl4 );
		Face <double> * f5 = new Face <double> ( &ppoints5, &hl5 );

		double f1_area = FaceArea::getFaceAreJ ( f1 );
		double f2_area = FaceArea::getFaceAreJ ( f5 );

		double res1 = TurningFunction::compare2FacesUsingTurningFunction<double> ( f1, f2, RotationInvariant, ScaleInvariant );
		double res2 = TurningFunction::compare2FacesUsingTurningFunction<double> ( f2, f1, RotationInvariant, ScaleInvariant );
		double res3 = TurningFunction::compare2FacesUsingTurningFunction<double>( f1, f2,RotationDependent, ScaleInvariant );
		double res4 = TurningFunction::compare2FacesUsingTurningFunction<double>( f3, f4,RotationDependent, ScaleInvariant );
		//double lkll = 37;
		std::cout << res1 << " " << res2 << '\n';

		ppoints5.push_back ( ppoints5[0] );
		ppoints5.erase ( ppoints5.begin() );

		faces.push_back (f1);
		DXFExport::exportFacesToDXF <double> ( "D:\\Tomas\\Cpp\\detectproj\\detectproj\\out\\f1_tf.dxf", &faces );
		}

		/*
		const char * eq = "R*ln(tg(u/2+ 45))";
		char operand[32];
		while (*eq != 0)
		{
		ArithmeticParser::findSequence(&eq, operand);
		// += strlen(operand);
		std::cout << operand << '\n' ;
		}
		*/
		/*
		TTurningFunction <double>::Type tff1, tff2;
		tff1[0.0]= 5;
		tff1[0.1] = 10;
		tff1[0.2] = 20;
		tff1[0.6] = 25;
		tff1[0.8] = 15;
		tff1[1.0] = 5;

		tff2[0.0]= 15;
		tff2[0.3] = 25;
		tff2[0.5] = 15;
		tff2[0.6] = 30;
		tff2[0.7] = 15;
		tff2[0.9] = 20;
		tff2[1.0] = 15;

		double dd1 = TurningFunction::compareTurningFunctions2<double>(tff1, tff2, RotationInvariant);
		*/

		/*
		TTurningFunction <double>::Type tf1, tf2, tf3, tf4, tf5, tf6;

		tf1[0.00000] =   291.29336;
		tf1[0.39120] =   202.19054;
		tf1[0.41608] =   184.70378;
		tf1[0.51914] =   114.57959;
		tf1[0.84750] =   54.71410;
		tf1[0.90812] =   21.88901;
		tf1[1.00000] =   291.29336;

		tf2[0.00000] =   288.59129;
		tf2[0.38346] =   205.42527;
		tf2[0.43423] =   187.40149;
		tf2[0.49437] =   118.79139;
		tf2[0.81925] =   44.18262;
		tf2[0.95925] =   15.21895;
		tf2[1.00000] =  288.59129;

		tf3[0.00000] =    21.88901;
		tf3[0.09188] =    291.29336;
		tf3[0.48308] =    202.19054;
		tf3[0.50796] =    184.70378;
		tf3[0.61102] =    114.57959;
		tf3[0.93938] =    54.71410;
		tf3[1.00000] =    21.88901;

		tf4[0.00000] =    15.21895;
		tf4[0.04075] =    288.59129;
		tf4[0.42421] =    205.42527;
		tf4[0.47498] =    187.40149;
		tf4[0.53512] =    118.79139;
		tf4[0.86000] =    44.18262;
		tf4[1.00000] =    15.21895;

		tf5[0.00000] = 250/100.0;
		tf5[0.1] =    20/100.0;
		tf5[0.2] =    300/100.0;
		tf5[0.4] =    20/100.0;
		tf5[0.5] =    300/100.0;
		tf5[0.85] =   200/100.0;
		tf5[1.00000] = 610/100.0;

		tf6[0.00000] = 180/100.0;
		tf6[0.3] =     250/100.0;
		tf6[0.8] =     180/100.0;
		tf6[0.9] =     150/100.0;
		tf6[1.00000] = 640/100.0;

		/*
		double d1 = TurningFunction::compareTurningFunctions<double>(tf1, tf4, RotationDependent );
		double d2 = TurningFunction::compareTurningFunctions<double>(tf3, tf2, RotationDependent );
		double d3 = TurningFunction::compareTurningFunctions2<double>(tf1, tf4, RotationDependent );
		double d4 = TurningFunction::compareTurningFunctions2<double>(tf3, tf2, RotationDependent );
		*/
		//double d5 = TurningFunction::compareTurningFunctions<double>(tf5, tf6, RotationDependent );
		//double d6 = TurningFunction::compareTurningFunctions<double>(tf5, tf6, RotationDependent );

		//std::cout << d2;

		//PointEllipsePosition::getPointEllipsePosition(&Point3DCartesian <double>(7,-7), 0, 0, 10, 5, 45);

		/*
		Container <Projection <double> *> proj_list;
		std::cout << ">> Loading cartographic projections...";
		proj_list.load ( "sinu.txt" );
		std::cout << " Completed." << std::endl;
		std::cout << proj_list.size() << " cartographic projection(s) have been loaded." << std::endl << std::endl;
		*output << proj_list.size() << " cartographic projection(s) have been loaded." << std::endl << std::endl;

		double fi = 60, la = -60;
		double h = CartDistortion::H (0.01, proj_list[0]->getXEquat(), proj_list[0]->getYEquat(), fi, la, proj_list[0]->getR(), proj_list[0]->getA(), proj_list[0]->getB(),
		proj_list[0]->getDx(), proj_list[0]->getDy(), proj_list[0]->getLat0(), proj_list[0]->getLat1(), proj_list[0]->getLat2(), proj_list[0]->getLon0() );
		double k = CartDistortion::K (0.01, proj_list[0]->getXEquat(), proj_list[0]->getYEquat(), fi, la, proj_list[0]->getR(), proj_list[0]->getA(), proj_list[0]->getB(),
		proj_list[0]->getDx(), proj_list[0]->getDy(), proj_list[0]->getLat0(), proj_list[0]->getLat1(), proj_list[0]->getLat2(), proj_list[0]->getLon0() );
		double s = CartDistortion::S (0.01, proj_list[0]->getXEquat(), proj_list[0]->getYEquat(), fi, la, proj_list[0]->getR(), proj_list[0]->getA(), proj_list[0]->getB(),
		proj_list[0]->getDx(), proj_list[0]->getDy(), proj_list[0]->getLat0(), proj_list[0]->getLat1(), proj_list[0]->getLat2(), proj_list[0]->getLon0() );
		double theta = CartDistortion::Theta (0.01, proj_list[0]->getXEquat(), proj_list[0]->getYEquat(), fi, la, proj_list[0]->getR(), proj_list[0]->getA(), proj_list[0]->getB(),
		proj_list[0]->getDx(), proj_list[0]->getDy(), proj_list[0]->getLat0(), proj_list[0]->getLat1(), proj_list[0]->getLat2(), proj_list[0]->getLon0() );
		TTissotIndicatrix <double> tiss = CartDistortion::Tiss (0.01, proj_list[0]->getXEquat(), proj_list[0]->getYEquat(), fi, la, proj_list[0]->getR(), proj_list[0]->getA(), proj_list[0]->getB(),
		proj_list[0]->getDx(), proj_list[0]->getDy(), proj_list[0]->getLat0(), proj_list[0]->getLat1(), proj_list[0]->getLat2(), proj_list[0]->getLon0() );

		short tt1 = PointEllipsePosition::getPointEllipsePosition(&Point3DCartesian <double>(3339510, -6679030), 4469420, -6122440,
		1000000*tiss.a_tiss, 1000000 * tiss.b_tiss, tiss.Ae_proj + tiss.b_mer - 90 );
		short tt2 = PointEllipsePosition::getPointEllipsePosition(&Point3DCartesian <double>(3339510, -6679030), 4469420, -6782217,
		1000000*tiss.a_tiss, 1000000 * tiss.b_tiss, tiss.Ae_proj + tiss.b_mer - 90 );
		short tt3 = PointEllipsePosition::getPointEllipsePosition(&Point3DCartesian <double>(3339510, -6679030), 3339510, -5826014,
		1000000*tiss.a_tiss, 1000000 * tiss.b_tiss, tiss.Ae_proj + tiss.b_mer - 90 );

		short tt4 = PointEllipsePosition::getPointEllipsePosition(&Point3DCartesian <double>(3339510, -6679030), 3339510, -6114730,
		1000000*tiss.a_tiss, 1000000 * tiss.b_tiss, tiss.Ae_proj + tiss.b_mer - 90 );

		/*
		double lll = 3;
		*/
		/*

		Container < Point3DCartesian <double> > points1, points2;
		points1.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tfl3.txt");
		points2.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tfl4.txt");

		double res = TangentFunction::compare2PolyLinesUsingTangentFunction( &points1, &points2, RotationDependent );

		//Point3DCartesian <double> p1 (10,0);
		//Point3DCartesian <double> p2 (0,0);
		//Point3DCartesian <double> p3 (0,10);
		//double angle = Angle3Points::getAngle3Points(&p3, &p2, &p1);

		Container < Node3DCartesian <double> *> points1, points2;
		points1.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tf23.txt");
		points2.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tf24.txt");

		Container <HalfEdge<double> *> hl;
		Face <double> * f1 = new Face <double> (&points1, &hl);
		Face <double> * f2 = new Face <double> (&points2, &hl);

		double res = TangentFunction::compare2FacesUsingTangentFunction <double>( f1, f2,RotationDependent );
		double res2 = InnerDistance::compare2FacesUsingInnerDistances( f1, f2);
		double res3 = TARCriterion::compare2FacesUsingTARCriterion( f1, f2);

		Container < Point3DCartesian <double> > points1, points2;

		points1.load<Dim2D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tf17.txt");
		points2.load<Dim2D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\tf18.txt");

		TTangentFunction <double>::Type tf1, tf2;

		for (int i = 0; i < points1.size(); i++)
		{
		tf1[points1[i].getX()] = points1[i].getY();
		tf2[points2[i].getX()] = points2[i].getY();
		}

		double res = TangentFunction::compareTangentFunctions <double>( tf1, tf2 );
		*/
		/*
		Container < Node3DCartesian <double> *> points3, points4;

		points3.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\id_vor1.txt");
		points4.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\id_vor2.txt");

		//double res2 = TangentFunction::compare2PolyLinesUsingTangentFunction( &points3, &points4, RotationDependent );
		//double res2 = TARCriterion::compare2FacesUsingTARCriterion(f1, f2);

		Container <HalfEdge<double> *> hl;
		ProjectionAzimuthal <double> pp;
		Face <double> * f1 = new Face <double> (points3, hl);
		Face <double> * f2 = new Face <double> (points4, hl);
		double res3= InnerDistance::compare2FacesUsingInnerDistances<double>(f1, f2);

		/*
		Container <Node3DCartesian <double> *> points1, points2;
		Container <HalfEdge<double> *> hl;

		points1.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\id7.txt");
		points2.load<Dim3D> ("D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\id8.txt");

		Face <double> * f1 = new Face <double> (&points1, &hl);
		Face <double> * f2 = new Face <double> (&points2, &hl);

		double res = InnerDistance::compare2FacesUsingInnerDistances(f1, f2);
		double res2 = TARCriterion::compare2FacesUsingTARCriterion(f1, f2);
		*/
		/*
		Container <Node3DCartesian <double> *> global_points;
		global_points.push_back (new Node3DCartesian <double> ( 0, 0 ) );
		global_points.push_back (new Node3DCartesian <double> ( 200, 0 ) );
		global_points.push_back (new Node3DCartesian <double> ( 200, 100 ) );
		global_points.push_back (new Node3DCartesian <double> ( 0, 100 ) );

		Container <Node3DCartesian <double> *> local_points;
		local_points.push_back (new Node3DCartesian <double> ( -50, 100 ) );
		local_points.push_back (new Node3DCartesian <double> ( -50, -100 ) );
		local_points.push_back (new Node3DCartesian <double> ( 50, -100 ) );
		local_points.push_back (new Node3DCartesian <double> ( 50, 100 ) );

		//local_points.push_back (new Node3DCartesian <double> ( 0, 0 ) );

		Container <Node3DCartesian <double> *> transformed_points;

		HelmertTransformation2D::transformPoints( &global_points, &local_points, &transformed_points );
		transformed_points.print();
		double test = Transformation2D::getStandardDeviation( &global_points, &transformed_points );

		/*
		Container <Node3DCartesian <double> *> points, points2;
		points.push_back (new Node3DCartesian <double> ( -2.5981, 1.500 ) );
		points.push_back (new Node3DCartesian <double> ( 0, 0 ) );
		points.pFush_back (new Node3DCartesian <double> (1.5, 2.5981 ) );
		points.push_back (new Node3DCartesian <double> (-1.0981, 4.0981) );
		points.push_back (new Node3DCartesian <double> (-0.5490, 2.049 ) );

		points2.push_back (new Node3DCartesian <double> ( 0, 0 ) );
		points2.push_back (new Node3DCartesian <double> (1.5, 2.5981 ) );
		points2.push_back (new Node3DCartesian <double> (-1.0981, 4.0981) );
		points2.push_back (new Node3DCartesian <double> (-0.5490, 2.049 ) );
		points2.push_back (new Node3DCartesian <double> ( -2.5981, 1.500 ) );

		Container <HalfEdge <double> *> edges1, edges2;
		Face <double> *f1 = new Face <double> (&points, &edges1);
		Face <double> *f2 = new Face <double> (&points2, &edges2);

		double res = TangentFunction::compare2FacesUsingTangentFunction(f1, f2, RotationDependent);
		*/

		//Test number of arguments
		//if ( argc < 3 ) throw Error ( "Error: not enough arguments in command line. " );

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat50\\m50000000\\Noise0.001\\Grid\\t1\\sinu\\t1_grid_test_sinu#100.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat50\\m50000000\\Noise0.001\\Grid\\t1\\t1_grid_ref#100.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tiha\\reichard\\reichard_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tiha\\reichard\\reichard_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\merc_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\merc_ref.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Portulan\\M2\\BFGS\\portulan_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Portulan\\M2\\BFGS\\portulan_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tiha\\albrizzi_africa\\albrizzi_test.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tiha\\albrizzi_africa\\albrizzi_reference.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat40\\m100000\\Noise0.001\\Grid\\t1\\sinu\\t1_grid_test_sinu#1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat40\\m100000\\Noise0.001\\Grid\\t1\\t1_grid_ref#1.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat10\\m5000000\\Noise0.002\\Grid\\t1\\sinu\\t1_grid_test_sinu#7.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat10\\m5000000\\Noise0.002\\Grid\\t1\\t1_grid_ref#7.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat30\\m50000000\\Noise0.003\\Meridians\\t1\\sinu\\t1_circle_test_sinu#9.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\lat30\\m50000000\\Noise0.003\\Meridians\\t1\\t1_circle_ref#9.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\m5\\Noise0.001\\Grid\\t1\\sinu\\t1_grid_test_sinu#1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\geoinformatica\\m5\\Noise0.001\\Grid\\t1\\t1_grid_ref#1.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\points_test_merc3.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\data\\points_reference_merc3.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\eck5\\t1_grid_test_eck5#-80_0_80_160.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t4\\aeqd\\t4_grid_decr_prime_mer_test_aeqd#-80_0_80_0.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t4\\t4_grid_decr_prime_mer_ref#-80_0_80_0.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t3\\aeqd\\t3_grid_decr_equat_test_aeqd#-50_-160_50_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t3\\t3_grid_decr_equat_ref#-50_-160_50_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\aeqd\\t2_grid_decr_prop_test_aeqd#-40_-60_40_60.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\t2_grid_decr_prop_ref#-40_-60_40_60.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t5\\aeqd\\t5_grid_mov_test_aeqd#-10_-100_-5_-95.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t5\\t5_grid_mov_ref#-10_-100_-5_-95.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t5\\eck5\\t5_grid_mov_test_eck5#-10_-100_-5_-95.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t5\\t5_grid_mov_ref#-10_-100_-5_-95.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1\\Grid\\t1\\eck5\\t1_grid_test_eck5#-80_80_0_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1\\Grid\\t1\\t1_grid_ref#-80_80_0_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1\\Grid\\t5\\eck5\\t5_grid_mov_test_eck5#-10_100_0_110.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1\\Grid\\t5\\t5_grid_mov_ref#-10_100_0_110.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\aeqd\\t2_grid_decr_prop_test_aeqd#-10_0_10_20.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\t2_grid_decr_prop_ref#-10_0_10_20.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\aeqd\\t2_grid_decr_prop_test_aeqd#-30_0_30_60.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t2\\t2_grid_decr_prop_ref#-30_0_30_60.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t1\\aeqd\\t1_grid_test_aeqd#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t1\\wer\\t1_grid_test_wer#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\aeqd\\t1_grid_test_aeqd#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\eck5\\t1_grid_test_eck5#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\sinu\\t1_grid_test_sinu#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t3\\merc\\t3_grid_decr_equat_test_merc#-20_0_20_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t3\\t3_grid_decr_equat_ref#-20_0_20_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t1\\laea\\t1_grid_test_laea#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\Noise0.25\\Random\\t1\\laea\\t1_random_test_laea#-80_0_80_160.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\Noise0.25\\Random\\t1\\laea\\t1_random_ref#-80_0_80_160.txt";


		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\-5_5\\merc\\t5_grid_mov_test_merc#-5_5_5_15_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\-5_5\\t5_grid_mov_ref#-5_5_5_15_1.txt";

		//.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\45_55\\merc\\t5_grid_mov_test_merc#45_55_55_65_20.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\45_55\\t5_grid_mov_ref#45_55_55_65_20.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\25_25\\merc\\t5_grid_mov_test_merc#25_25_35_35_19.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise50\\Grid\\t5\\25_25\\t5_grid_mov_ref#25_25_35_35_19.txt";


		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic2\\habilitation\\t5_minus2\\Noise50\\Grid\\t5\\65_65\\eck5\\t5_grid_mov_test_eck5#65_65_75_75_10.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic2\\habilitation\\t5_minus2\\Noise50\\Grid\\t5\\65_65\\t5_grid_mov_ref#65_65_75_75_10.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic2\\habilitation\\t5_minus2\\Noise50\\Grid\\t5\\15_15\\eck5\\t5_grid_mov_test_eck5#15_15_25_25_10.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic2\\habilitation\\t5_minus2\\Noise50\\Grid\\t5\\15_15\\t5_grid_mov_ref#15_15_25_25_10.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Random\\t6\\eck5\\t6_random_test_eck5#-80_-160_80_160_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Random\\t6\\t6_random_ref#-80_-160_80_160_1.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise0\\Grid\\t5\\45_55\\merc\\t5_grid_mov_test_merc#45_55_55_65_9.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t5_merc_new\\Noise0\\Grid\\t5\\45_55\\t5_grid_mov_ref#45_55_55_65_9.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Random\\t6\\laea\\t6_random_test_laea#-80_-160_80_160_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Random\\t6\\t6_random_ref#-80_-160_80_160_1.txt";


		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Grid\\t6\\merc\\t6_grid_test_merc#-80_-160_80_160_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6\\Noise90\\Grid\\t6\\t6_grid_ref#-80_-160_80_160_1.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_merc\\Noise90\\Random\\t6\\merc\\t6_random_test_merc#-80_-160_80_160_11.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_merc\\Noise90\\Random\\t6\\t6_random_ref#-80_-160_80_160_11.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_laea\\Noise90\\Random\\t6\\laea\\t6_random_test_laea#-80_-160_80_160_15.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_laea\\Noise90\\Random\\t6\\t6_random_ref#-80_-160_80_160_15.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_minus2\\Noise2.5\\Random\\t6\\wer\\t6_random_test_wer#-80_-160_80_160_15.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_minus2\\Noise2.5\\Random\\t6\\t6_random_ref#-80_-160_80_160_15.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_laea\\Noise90\\Random\\t6\\laea\\t6_random_test_laea#-80_-160_80_160_15.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_laea\\Noise90\\Random\\t6\\t6_random_ref#-80_-160_80_160_15.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t1\\Noise100\\Random\\t1\\eck5\\t1_random_test_eck5#-80_0_80_160_15.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t1\\Noise100\\Random\\t1\\t1_random_ref#-80_0_80_160_15.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_minus2\\Noise2.5\\Random\\t6\\eck5\\t6_random_test_eck5#-80_-160_80_160_15.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t6_minus2\\Noise2.5\\Random\\t6\\t6_random_ref#-80_-160_80_160_15.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-10\\merc\\t7_random_decr_prop_test_merc#-10_-160_10_20_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-10\\t7_random_decr_prop_ref#-15_-160_10_20_1.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-80\\merc\\t7_random_decr_prop_test_merc#-80_-160_80_160_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-80\\t7_random_decr_prop_ref#-80_-160_80_160_1.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-80\\merc\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-80\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file =      "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk6\\m6\\bfgs\\jtsk_test6.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk6\\m6\\bfgs\\jtsk_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk8\\m2\\bfgs\\jtsk_test8.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\jtsk\\jtsk8\\m2\\bfgs\\jtsk_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\M6\\BFGS\\s42_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\other\\s42\\M6\\BFGS\\s42_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M6\\BFGS\\hondius_world_east_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_east\\M6\\BFGS\\hondius_world_east_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_west\\M6\\BFGS\\hondius_world_west_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Hondius_world_west\\M6\\BFGS\\hondius_world_west_reference.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\America_Seutter\\M6\\BFGS\\america_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\America_Seutter\\M6\\BFGS\\america_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Safarik\\M6\\BFGS\\safarik_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Safarik\\M6\\BFGS\\safarik_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Africa_Delisle\\M6\\BFGS\\africa_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\map_collection_uk\\Africa_Delisle\\M6\\BFGS\\africa_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\GreatBrittain\\M6\\BFGS\\brittain_test.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\davidrumsay\\GreatBrittain\\M6\\BFGS\\brittain_reference.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t8_minus2\\Noise90\\Random\\t8\\-80\\sinu\\t8_random_decr_equat_test_sinu#-80_-160_80_160_1.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t8_minus2\\Noise90\\Random\\t8\\-80\\t8_random_decr_equat_ref#-80_0_80_160_15.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t9_minus2\\Noise10\\Random\\t9\\0\\wer\\t9_random_decr_central_mer_test_wer#-80_0_80_0_1.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t9_minus2\\Noise10\\Random\\t9\\0\\t9_random_decr_central_mer_ref#-80_0_80_0_1.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t9_minus2\\Noise90\\Random\\t9\\20\\eck5\\t9_random_decr_central_mer_test_eck5#-80_0_80_20_1.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t9_minus2\\Noise90\\Random\\t9\\20\\t9_random_decr_central_mer_ref#-80_0_80_20_1.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\eqdc\\Terr2\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\eqdc\\Terr2\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\eqdc\\t13_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\eqdc\\t13_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\merc\\t13_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\merc\\t13_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\laea\\t13_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80\\laea\\t13_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80_1\\eqdc\\GND500_2\\t13_random_decr_prop_test_eqdc#-80_-160_80_160_1.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t13\\Noise500\\Random\\t13\\-80_1\\eqdc\\GND500_2\\t13_random_decr_prop_ref#-80_-160_80_160_1.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\laea\\Terr2\\M8\\Minus1\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\laea\\Terr2\\M8\\Minus1\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr2\\M8\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr2\\M8\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr2\\M8\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr2\\M8\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr1\\M7s\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr1\\M7s\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr1\\M7\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr1\\M7\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr2\\M7\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\eqdc\\Terr2\\M7\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

/*
analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\eqdc\\Terr2\\M6\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\eqdc\\Terr2\\M6\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";


analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\laea\\Terr2\\M6\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\laea\\Terr2\\M6\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";


analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\merc\\Terr2\\M6\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\merc\\Terr2\\M6\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";
*/
//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\eqdc\\Terr1\\M7\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t17_18\\eqdc\\Terr1\\M7\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr1\\M7\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t15_16\\merc\\Terr1\\M7\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t19_20\\laea\\Terr1\\M7\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t19_20\\laea\\Terr1\\M7\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t23_24\\eqdc\\Terr1\\M7\\target1\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t23_24\\eqdc\\Terr1\\M7\\target1\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\DE\\M7S\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\DE\\M7S\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr1\\DE\\M7S\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr1\\DE\\M7S\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr1\\DE\\M7S\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr1\\DE\\M7S\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\DE\\M7S\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\DE\\M7S\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\BFGSH\\Terr3\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\BFGSH\\Terr3\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\BFGSH\\Terr3\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\BFGSH\\Terr3\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\Nelder\\Terr3\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\Nelder\\Terr3\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\DE\\M7S\\t7_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr1\\DE\\M7S\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

 //analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\Nelder\\Terr4\\t7_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr1\\Nelder\\Terr4\\t7_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\BFGSH\\Terr4\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\BFGSH\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\Nelder\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\BFGSH\\Terr4\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\BFGSH\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr2\\DE\\Terr4\\t2_random_decr_prop_test_laea#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\laea\\Terr2\\DE\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\DE\\Terr4\\t2_random_decr_prop_test_eqdc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\eqdc\\Terr2\\DE\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\DE\\Terr4\\t2_random_decr_prop_test_merc#-80_-160_80_160_10.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t27_28\\merc\\Terr2\\DE\\Terr4\\t2_random_decr_prop_ref#-80_-160_80_160_10.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tobler\\projA\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tobler\\projA\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tobler\\projB\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tobler\\projB\\reference.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\AmericaSanson\\Graticule\\BFGSH\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\AmericaSanson\\Graticule\\BFGSH\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\Italy\\Graticule\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\Italy\\Graticule\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\NorthernRegions\\Graticule\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\ContentVsGraticule\\NorthernRegions\\Graticule\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\AirRoutes\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\AirRoutes\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Mercator\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Mercator\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Ortelius\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Ortelius\\reference.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Tectonico\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\WorldPolitical\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\WorldPolitical\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Agriculture\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\WorldMaps\\Agriculture\\reference.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Hondius_world\\testw.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Hemisphere\\Hondius_world\\referencew.txt";

//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Asia\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Asia\\reference.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Europe\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\Continents\\Europe\\reference.txt";


//analysis_parameters.test_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\CentralAmerica\\test.txt";
//analysis_parameters.reference_file = "E:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\Habilitation\\LargeTerritories\\CentralAmerica\\reference.txt";		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-70\\merc\\t7_random_decr_prop_test_merc#-70_-160_70_140_1.txt";
		
//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-70\\t7_random_decr_prop_ref#-70_-160_70_140_1.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-60\\merc\\t7_random_decr_prop_test_merc#-60_-160_60_120_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-60\\t7_random_decr_prop_ref#-60_-160_60_120_1.txt";


		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-50\\merc\\t7_random_decr_prop_test_merc#-50_-160_50_100_1.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\synthetic\\habilitation\\t7_merc\\Noise90\\Random\\t7\\-50\\t7_random_decr_prop_ref#-50_-160_50_100_1.txt";


		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\t1\\merc\\t1_grid_test_merc#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\merc\\t2_grid_decr_prop_test_merc#-70_0_70_140.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\t2_grid_decr_prop_ref#-70_0_70_140.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\merc\\t2_grid_decr_prop_test_merc#-30_0_30_60.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\t2_grid_decr_prop_ref#-30_0_30_60.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\merc\\t2_grid_decr_prop_test_merc#-60_0_60_120.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t2\\t2_grid_decr_prop_ref#-60_0_60_120.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\t1\\merc\\t1_grid_test_merc#-30_0_30_60.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\t1\\t1_grid_ref#-30_0_30_60.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.5\\Grid\\t1\\eck5\\t1_grid_test_eck5#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.5\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t1\\laea\\t1_grid_test_laea#-80_0_80_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.25\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\mercator2.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\testik\\mercator1.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.75\\Grid\\t4\\wer\\t4_grid_decr_prime_mer_test_wer#-80_0_80_20.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.75\\Grid\\t4\\t4_grid_decr_prime_mer_ref#-80_0_80_20.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t3\\eqdc\\t3_grid_decr_equat_test_eqdc#-10_0_10_160.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise1.0\\Grid\\t3\\t3_grid_decr_equat_ref#-10_0_10_160.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.75\\Grid\\t5\\sinu\\t5_grid_mov_test_sinu#60_40_70_50.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.75\\Grid\\t5\\t5_grid_mov_ref#60_40_70_50.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.53\\Grid\\t5\\sinu\\t5_grid_mov_test_sinu#60_60_65_65.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.53\\Grid\\t5\\t5_grid_mov_ref#60_60_65_65.txt";

		//test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.5\\Grid\\t5\\sinu\\t5_grid_mov_test_sinu#70_10_80_20.txt";
		//reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0.5\\Grid\\t5\\t5_grid_mov_ref#70_10_80_20.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise50\\Grid\\t5\\sinu\\t5_grid_mov_test_sinu#40_110_45_115.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise50\\Grid\\t5\\t5_grid_mov_ref#40_110_45_115.txt";

		//analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\eck5\\t1_grid_test_eck5#-80_0_80_160.txt";
		//analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		/*
		analysis_parameters.test_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\eqdc\\t1_grid_test_eqdc#-80_0_80_160.txt";
		analysis_parameters.reference_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\test\\Noise0\\Grid\\t1\\t1_grid_ref#-80_0_80_160.txt";

		analysis_parameters.match_method = MatchTissotIndicatrix;
		analysis_parameters.analyze_oblique_aspect = true;
		analysis_parameters.analysis_type.a_cnd = true;
		analysis_parameters.analysis_type.a_and = true;
		analysis_parameters.analysis_type.a_homt = true;
		analysis_parameters.analysis_type.a_helt = true;
		analysis_parameters.analysis_type.a_gn_ad = false;
		analysis_parameters.analysis_type.a_gn_tf = false;
		analysis_parameters.analysis_type.a_nnng = true;
		analysis_parameters.analysis_type.a_sig = true;
		analysis_parameters.analysis_type.a_vd_al = true;
		analysis_parameters.analysis_type.a_vd_tf = true;
		analysis_parameters.analysis_type.a_vd_tar = true;
		analysis_parameters.analysis_type.a_vd_id = true;
		analysis_parameters.perform_heuristic = true;
		*/

		//Process all parameters and arguments of the command line
		{
			bool proj_set = false, lat0_set = false, latp_set = false, lonp_set = false, projections_file_set = false;

			while ( --argc > 0 )
			{
				//Get - (follows some parameter)
				if ( *argv[argc] == '-' )
				{
					//Process parameter after -
					for ( char * parameter = argv[argc]; ; )
					{
						switch ( * ( ++ parameter ) )
						{
							//Help
						case '?':
							{
								//Print help
								printHelp ( &std::cout );

								//Stop analysis_type
								return 0;
							}

							//Terminate character \0 of the argument
						case '\0':
							{
								break;
							}

							//Perform heuristic: throw nonperspective samples
						case 'h':
							{
								analysis_parameters.perform_heuristic = true;
								continue;
							}

							//Detect projection in normal aspect
						case 'n':
							{
								analysis_parameters.analyze_normal_aspect = true;
								continue;
							}

							//Detect projection in transverse aspect
						case 't':
							{
								analysis_parameters.analyze_transverse_aspect = true;
								continue;
							}

							//Detect projection in oblique aspect
						case 'o':
							{
								analysis_parameters.analyze_oblique_aspect = true;
								continue;
							}

							//Remove outliers
						case 'r':
							{
								analysis_parameters.remove_outliers = true;
								continue;
							}

							//Correct rotation (inappropriate rotated samples)
						case 'c':
							{
								analysis_parameters.correct_rotation = true;
								continue;
							}

							//Print warnings and exceptions
						case 'w':
							{
								analysis_parameters.print_exceptions = true;
								continue;
							}

							//Throw exception
						default:
							{
								std::cout << parameter << '\n';
								throw Error ( "Error: Invalid parameter in command line!" );
							}
						}

						//Stop processing of the argument
						break;
					}
				}

				//Get + ( follows some type of cartometric analysis_type and analyzed projection definition )
				else if ( *argv[argc] == '+' )
				{
					//Get new command line parameter
					char * attribute = const_cast <char *> ( argv[argc] + 1 ),  *value = NULL;
					char * command = attribute;

					//Find splitter: =
					for ( ; ( *command != '=' ) && ( *command != '\0' ) ; command++ );

					//We found splitter, trim command and copy to value
					if ( ( *command == '=' ) && ( *command != '\0' ) )
					{
						*command = '\0';
						value = command + 1;
					}

					//Throw exception
					if ( attribute == NULL || value == NULL ) throw Error ( "Error: Invalid value in command line!" );

					//Set analyzed projection
					if ( !strcmp ( "proj", attribute ) )
					{
						//Current proj, set proj name
						if ( !proj_set ) proj_set = true;

						//Skip current projection and start next projection
						else
						{
							//Add old projection to the list
							analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );

							//Reset other parameters
							lat0_set = false; latp_set = false; lonp_set = false;
						}

						//Set proj name
						strcpy ( analyzed_proj_parameters.proj_name, value );
					}

					//Set detection method
					else if ( !strcmp ( "met", attribute ) )
					{
						if ( !strcmp ( "si", value ) || !strcmp ( "m1", value ) ) analysis_parameters.analysis_method = SimplexMethod;
						else if ( !strcmp ( "de", value ) || !strcmp ( "m2", value ) ) analysis_parameters.analysis_method  = DifferentialEvolutionMethod;
						else if ( !strcmp ( "ls", value ) || !strcmp ( "m3", value ) ) analysis_parameters.analysis_method  = NonLinearLeastSquaresMethod;
						else if ( !strcmp ( "sir", value ) || !strcmp ( "m4", value ) ) analysis_parameters.analysis_method  = SimplexRotMethod;
						else if ( !strcmp ( "der", value ) || !strcmp ( "m5", value ) ) analysis_parameters.analysis_method  = DifferentialEvolutionRotMethod;
						else if ( !strcmp ( "lsr", value ) || !strcmp ( "m6", value ) ) analysis_parameters.analysis_method  = NonLinearLeastSquaresRotMethod;
						else if (!strcmp("sir2", value) || !strcmp("m7", value)) analysis_parameters.analysis_method = SimplexRot2Method;
						else if (!strcmp("der2", value) || !strcmp("m8", value)) analysis_parameters.analysis_method = DifferentialEvolutionRot2Method;
						else if (!strcmp("lsr2", value) || !strcmp("m9", value)) analysis_parameters.analysis_method = NonLinearLeastSquaresRot2Method;
						else if (!strcmp("sis", value) || !strcmp("m10", value)) analysis_parameters.analysis_method = SimplexShiftsMethod;
						else if (!strcmp("des", value) || !strcmp("m11", value)) analysis_parameters.analysis_method = DifferentialEvolutionShiftsMethod;
						else if (!strcmp("lss", value) || !strcmp("m12", value)) analysis_parameters.analysis_method = NonLinearLeastSquaresShiftsMethod;
						else throw Error ( "Error: Invalid analysis_method type in command line!" );
					}

					//Set type of analysis
					else if ( !strcmp ( "an", attribute ) )
					{
						if ( !strcmp ( "cnd", value ) ) analysis_parameters.analysis_type.a_cnd = true;
						else if ( !strcmp ( "homt", value ) ) analysis_parameters.analysis_type.a_homt = true;
						else if ( !strcmp ( "helt", value ) ) analysis_parameters.analysis_type.a_helt = true;
						else if ( !strcmp ( "gntf", value ) ) analysis_parameters.analysis_type.a_gn_tf = true;
						else if ( !strcmp ( "vdtf", value ) ) analysis_parameters.analysis_type.a_vd_tf = true;
						else if ( !strcmp ( "all", value ) )
						{
							analysis_parameters.analysis_type.a_cnd = true; analysis_parameters.analysis_type.a_homt = true; analysis_parameters.analysis_type.a_helt = true;
							analysis_parameters.analysis_type.a_gn_tf = true; analysis_parameters.analysis_type.a_vd_tf = true; analysis_parameters.analysis_type.a_vd_id = true;
						}

						//Switch both homt and helt analyzes for the rotation correction
						else if ( analysis_parameters.correct_rotation )
						{
							analysis_parameters.analysis_type.a_homt = true;
							analysis_parameters.analysis_type.a_helt = true;
						}
						else throw Error ( "Error: Invalid analysis_type type in command line!" );
					}

					//Export graticule of the first k samples into dxf file
					else if ( !strcmp ( "gr", attribute ) )
					{
						//Set amount of samples exported to DXF
						analysis_parameters.exported_graticule = atof ( value );

						//Enable homothetic transformation for DE method used for estiamtion of R, dx, dy
						if ( ( analysis_parameters.analysis_method == DifferentialEvolutionMethod ) && ( !analysis_parameters.analysis_type.a_homt ) )
						{
							analysis_parameters.analysis_type.a_homt = true;
						}
					}

					//Latitude of the cartographic pole in analyzed projection
					else if ( !strcmp ( "latp", attribute ) )
					{
						//Current proj, set latp
						if ( !latp_set ) latp_set = true;

						//Skip current projection and start next projection
						else
						{
							//Projection not set: bad parameter in command line
							if ( !proj_set ) throw Error ( "Error: Invalid projection ID in command line!" );

							//Add old projection to the list
							analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );

							//Reset old values
							proj_set = false; lat0_set = false; lonp_set = false;
						}

						//Set latp
						analyzed_proj_parameters.latp = atof ( value );

						//Correct value
						if ( fabs ( analyzed_proj_parameters.latp ) > MAX_LAT )
							analyzed_proj_parameters.latp = 90.0;
					}

					//Longitude of the cartographic pole in analyzed projection
					else if ( !strcmp ( "lonp", attribute ) )
					{
						//Current proj, set lonp
						if ( !lonp_set )  lonp_set = true;

						//Skip current projection and start next projection
						else
						{
							//Projection not set: bad parameter in command line
							if ( !proj_set ) throw Error ( "Error: Invalid proj type in command line!" );

							//Add old projection to the list
							analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );

							//Reset old values
							proj_set = false; latp_set = false; lat0_set = false;
						}

						//Set lonp
						analyzed_proj_parameters.lonp = atof ( value );

						//Correct value
						if ( fabs ( analyzed_proj_parameters.lonp ) > MAX_LON )
							analyzed_proj_parameters.lonp = 0.0;
					}

					//Latitude of the undistorted parallel in analyzed projection
					else if ( !strcmp ( "lat0", attribute ) )
					{
						//Current proj, set lat0
						if ( !lat0_set ) lat0_set = true;

						//Skip previous projection and start next projection
						else
						{
							//Projection not set: bad parameter in command line
							if ( !proj_set ) throw Error ( "Error: Invalid proj type in command line!" );

							//Add old projection to the list
							analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );

							//Reset old values
							proj_set = false; latp_set = false; lonp_set = false;
						}

						//Set lat0
						analyzed_proj_parameters.lat0 = atof ( value );

						//Correct value
						if ( fabs ( analyzed_proj_parameters.lat0 ) > MAX_LAT0 )
							analyzed_proj_parameters.lat0 = 45.0;
					}

					//Set ä new central meridian shifted by lon0
					else if ( !strcmp ( "lon0", attribute ) )
					{
						analysis_parameters.lon0 = atof ( value );
					}

					//Set the latitude increment of adjacent parallels in the graticule
					else if ( !strcmp ( "dlat", attribute ) )
					{
						analysis_parameters.lat_step = atof ( value );

						//Correct value
						if ( analysis_parameters.lat_step  > MAX_LAT_STEP )
							analysis_parameters.lat_step  = MAX_LAT_STEP;
						else if ( analysis_parameters.lat_step  < MIN_LAT_STEP )
							analysis_parameters.lat_step  = MIN_LAT_STEP;
					}

					//Set the longitude increment of adjacent meridians in the graticule
					else if ( !strcmp ( "dlon", attribute ) )
					{
						analysis_parameters.lon_step = atof ( value );

						//Correct value
						if ( analysis_parameters.lon_step  > MAX_LON_STEP )
							analysis_parameters.lon_step  = MAX_LON_STEP;
						else if ( analysis_parameters.lon_step  < MIN_LON_STEP )
							analysis_parameters.lon_step  = MIN_LON_STEP;
					}

					//Set increment value for the latitude of the cartographic pole
					else if ( !strcmp ( "dlatp", attribute ) )
					{
						analysis_parameters.latp_step = atof ( value );
					}

					//Set increment value for the longitude of the cartographic pole
					else if ( !strcmp ( "dlonp", attribute ) )
					{
						analysis_parameters.lonp_step = atof ( value );
					}

					//Set increment value for the latitude of the undistorted parallel
					else if ( !strcmp ( "dlat0", attribute ) )
					{
						analysis_parameters.lat0_step = atof ( value );
					}

					//Set printed results
					else if ( !strcmp ( "res", attribute ) )
					{
						analysis_parameters.printed_results = atoi ( value );
					}

					//Set sensitivity of the heuristic
					else if ( !strcmp ( "sens", attribute ) )
					{
						analysis_parameters.heuristic_sensitivity_ratio  = atof ( value );

						//Correct value
						if ( analysis_parameters.heuristic_sensitivity_ratio > MAX_HEURISTIC_SENSTIVITY_RATIO )
							analysis_parameters.heuristic_sensitivity_ratio = MAX_HEURISTIC_SENSTIVITY_RATIO;
						else if ( analysis_parameters.heuristic_sensitivity_ratio < MIN_HEURISTIC_SENSTIVITY_RATIO )
							analysis_parameters.heuristic_sensitivity_ratio = MIN_HEURISTIC_SENSTIVITY_RATIO;
					}

					//Set number of repetition
					else if ( !strcmp ( "rep", attribute ) )
					{
						analysis_parameters.analysis_repeat  = ( int ) atof ( value );

						//Correct value
						if ( analysis_parameters.analysis_repeat > MAX_ANALYSIS_REPEAT )
							analysis_parameters.analysis_repeat = MAX_ANALYSIS_REPEAT;
						else if ( analysis_parameters.analysis_repeat < MIN_ANALYSIS_REPEAT )
							analysis_parameters.analysis_repeat = MIN_ANALYSIS_REPEAT;
					}

					//Set heuristic increment
					else if ( !strcmp ( "dsens", attribute ) )
					{
						analysis_parameters.heuristic_sensitivity_increment  = atof ( value );

						//Correct value
						if ( analysis_parameters.heuristic_sensitivity_increment > MAX_HEURISTIC_SENSTIVITY_ICREMENT )
							analysis_parameters.heuristic_sensitivity_increment = MAX_HEURISTIC_SENSTIVITY_ICREMENT;
						else if ( analysis_parameters.heuristic_sensitivity_increment < MIN_HEURISTIC_SENSTIVITY_INCREMENT )
							analysis_parameters.heuristic_sensitivity_increment = MIN_HEURISTIC_SENSTIVITY_INCREMENT;
					}

					//Set match type: circle or Tissot Indicatrix
					else if ( !strcmp ( "match", attribute ) )
					{
						if ( !strcmp ( "circ", value ) ) analysis_parameters.match_method = MatchCircle;
						else if ( !strcmp ( "tiss", value ) ) analysis_parameters.match_method = MatchTissotIndicatrix;

						else throw Error ( "Error: Invalid result match type type in command line!" );
					}

					//Set max residuals
					else if ( !strcmp ( "maxerror", attribute ) )
					{
						analysis_parameters.max_error = atof ( value );
					}

					//Bad argument
					else
					{
						std::cout << attribute << '\n';
						throw Error ( "Error: Invalid atribute in command line!" );
					}
				}

				//Process input / output file
				else
				{
					//Set reference file
					if ( analysis_parameters.reference_file == NULL )
					{
						analysis_parameters.reference_file = argv[argc];
					}

					//Set test file
					else if ( analysis_parameters.test_file == NULL )
					{
						analysis_parameters.test_file = argv[argc];
					}

					//Set projections file
					else if ( ( !strcmp ( "projections.txt", analysis_parameters.projections_file ) ) && ( !projections_file_set ) )
					{
						projections_file_set = true;
						analysis_parameters.projections_file = argv[argc];
					}

					//Throw exception
					else throw Error ( "Error: Too many input files in command line!" );
				}
			}

			//Add analyzed projection properties from command line, if necessary
			if ( proj_set || latp_set || lonp_set || lat0_set )
			{
				//Add analyzed projection to the list
				if ( proj_set ) analyzed_proj_parameters_list.push_back ( analyzed_proj_parameters );
				else throw Error ( "Error: Invalid proj type in command line!" );

				//Reset all values
				proj_set = false; lat0_set = false; latp_set = false; lonp_set = false;
			}
		}

		//Test input / output file
		if ( ( analysis_parameters.test_file == NULL ) || ( analysis_parameters.reference_file == NULL ) ) throw ErrorFileRead ( "ErrorFileRead, no input file set. ", "Analysis cancelled." );

		/*
		//Delete after batch tests
		//analysis_parameters.projections_file= "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\maps\\tiha\\reichard\\projections.txt";
		analysis_parameters.projections_file = "D:\\Tomas\\Cpp\\detectproj\\detectproj\\projects\\vs2010\\projections_hire5.txt";
		analysis_parameters.projections_file= "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\Geoinformatica\\projections_red.txt";
		//analysis_parameters.projections_file= "D:\\Tomas\\Cpp\\detectproj\\detectproj\\tests\\Geoinformatica\\projections_four.txt";
		analysis_parameters.analysis_type.a_cnd = false;
		analysis_parameters.analysis_type.a_vd_tf = false;
		analysis_parameters.analysis_type.a_gn_tf = false;
		analysis_parameters.match_method = MatchCircle;
		*/

		//Create output protocol for results of analysis
		char output_file_text[MAX_TEXT_FILE_LENGTH];
		strcpy ( output_file_text, analysis_parameters.test_file );
		strcat ( output_file_text, ".log" );

		//There is no protocol file, no analysis_type has been performed (jump files with created protocol)
		FILE * output_file_test = fopen ( output_file_text, "r" );

		//Log file exist, test if successfully computed (the last character is 1)
		if ( output_file_test )
		{
			//Go to end of the file
			if ( fseek ( output_file_test, -1, SEEK_END ) != 0 )
			{
				fclose ( output_file_test );
				throw ErrorFileRead ( "ErrorFileRead: can not process the file.", output_file_text );
			}

			//Read last char
			char last_char[1];

			if ( fread ( last_char, 1, 1, output_file_test ) != 1 )
			{
				fclose ( output_file_test );
				throw ErrorFileRead ( "ErrorFileRead: can not read the file.", output_file_text );
			}

			//Last protocol successfully computed and closed
			if ( *last_char == '1' ) return 0;
		}

		//Otherwise create output log file
		output_file.open ( output_file_text );

		//Create output file
		output = &output_file;

		//Print header
		//printHeader ( &std::cout );
		printHeader ( output );

		//All analysis diabled
		if ( !analysis_parameters.analysis_type.a_cnd && !analysis_parameters.analysis_type.a_homt &&  !analysis_parameters.analysis_type.a_helt &&
			!analysis_parameters.analysis_type.a_gn_tf && !analysis_parameters.analysis_type.a_vd_tf )
		{
			throw Error ( "Error: All analyzes diabled, can not compute cartometric analysis." );
		}

		*output << "Prerelease version, only for testing...\n";

		//At least one analysis is enabled, continue
		//Get start time
		time_t start, end;
		struct tm * start_ascii, *end_ascii;
		time ( &start );
		start_ascii = localtime ( &start );
		std::cout << ">> Starting...: " << asctime ( start_ascii ) << std::endl;
		*output << ">> Starting...: " << asctime ( start_ascii ) << std::endl;

		//Create new lists of test and reference points
		Container <Node3DCartesian <double> *> nl_test, nl_transformed;
		Container <Point3DGeographic <double> *> nl_reference;
		Container <Node3DCartesianProjected <double> *> nl_projected;

		std::cout << ">> Loading points...";
		nl_test.load <Dim2D> ( analysis_parameters.test_file );
		nl_reference.load <Dim2D> ( analysis_parameters.reference_file );
		std::cout << "Completed." << std::endl;

		//Get test / reference points count
		std::cout.flush();
		const unsigned int n_test = nl_test.size(), n_reference = nl_reference.size();
		std::cout << "Test points: " << analysis_parameters.test_file << std::endl;
		std::cout << "Reference points: " << analysis_parameters.reference_file << std::endl;
		std::cout << n_test << " test point(s) have been loaded." << std::endl;
		std::cout << n_reference << " reference point(s) have been loaded." << std::endl << std::endl;
		*output << "Test points: " << analysis_parameters.test_file << std::endl;
		*output << "Reference points: " << analysis_parameters.reference_file << std::endl;
		*output << n_test << " test point(s) have been loaded." << std::endl;
		*output << n_reference << " reference point(s) have been loaded." << std::endl;

		//Different number of test and reference points
		if ( n_test != n_reference )
		{
			throw ErrorBadData ( "ErroBadData: different number of test and reference points, ", "can not perform cartometric analysis." );
		}

		//Find and remove duplicate points
		//nl_test.removeDuplicateElements ( nl_test.begin(), nl_test.end(),  sortPointsByX (), isEqualPointByPlanarCoordinates <Node3DCartesian <double> *> () );
		//nl_reference.removeDuplicateElements ( nl_reference.begin(), nl_reference.end(), sortPointsByLat (), isEqualPointByPlanarCoordinates <Point3DGeographic <double> *> () );

		//Total unique points
		const unsigned int n_test_unique = nl_test.size(), n_reference_unique = nl_reference.size();

		//Different number of test and reference points
		if ( n_test_unique != n_reference_unique )
			throw ErrorBadData ( "ErroBadData: number of unique points in both files different, ", "can not perform cartometric analysis." );

		//Not enough unique points
		if ( n_test_unique < 3 )
			throw ErrorBadData ( "EErroBadData: insufficient number of unique points ( at least 3 ), ", "can not perform cartometric analysis." );

		//Total duplicate points
		if ( n_test_unique < n_test )
			*output << n_test - n_test_unique << " duplicate point(s) have been found. Removed... " << std::endl << std::endl;

		//Lat/lon values outside the limits
		for ( unsigned int i = 0; i < n_test_unique; i++ )
		{
			if ( ( fabs ( nl_reference[i]->getLat() ) > MAX_LAT ) || ( fabs ( nl_reference[i]->getLon() ) > MAX_LON ) )
				throw ErrorBadData ( "EErroBadData: latitude or longitude outside intervals, ", "can not perform cartometric analysis." );
		}

		//Load cartographic projections from file
		Container <Projection <double> *> proj_list;
		std::cout << ">> Loading cartographic projections...";
		proj_list.load ( analysis_parameters.projections_file );
		std::cout << " Completed." << std::endl;
		std::cout << proj_list.size() << " cartographic projection(s) have been loaded." << std::endl << std::endl;
		*output << proj_list.size() << " cartographic projection(s) have been loaded." << std::endl << std::endl;

		/*
		//Use a shifted central meridian: Subtract lon0 from lon
		if ( analysis_parameters.lon0 != 0.0 )
		{
		//Correct all points, set new central meridian
		for ( unsigned int i = 0; i < n_test_unique; i++ )
		{
		nl_reference [i]->setLon ( CartTransformation::redLon0 ( ( nl_reference ) [i]->getLon(), analysis_parameters.lon0 ) );
		}
		}
		*/

		//Check if analyzed projection is present in the list of cartographic projections
		//Process all extracted projections from command line
		for ( TAnalyzedProjParametersList <double> ::Type ::const_iterator i_analyzed_projections = analyzed_proj_parameters_list.begin(); i_analyzed_projections != analyzed_proj_parameters_list.end(); ++ i_analyzed_projections )
		{
			//Compare with all projections stored in the list of projections
			for ( TItemsList <Projection <double> *> ::Type ::const_iterator i_projections = proj_list.begin(); i_projections != proj_list.end(); ++ i_projections )
			{
				//Compare analyzed projection name with names of projections stored in the list
				if ( !strcmp ( ( * i_analyzed_projections ).proj_name, ( *i_projections )->getProjectionName() ) )
				{
					//Clone analyzed projection
					Projection <double> *analyzed_proj = ( *i_projections )->clone();

					//Set analyzed projections properties: if not set by user in command line, use default projection values from configuration file
					//We set lat0 in command line and analyzed projection does not define own lat0 in configuration file
					if ( ( fabs ( ( * i_analyzed_projections ).lat0 ) <= MAX_LAT0 ) && ( analyzed_proj->getLat0() == 45.0 ) )
						analyzed_proj->setLat0 ( ( * i_analyzed_projections ).lat0 );

					//We set latp, lonp in command and analyzed projection does not define own latp,lonp in configuration file (have a normal aspect)
					if ( ( fabs ( ( * i_analyzed_projections ).latp ) < MAX_LAT ) && ( fabs ( ( * i_analyzed_projections ).lonp ) <= MAX_LON ) && ( analyzed_proj->getCartPole().getLat() == MAX_LAT ) )
					{
						Point3DGeographic <double> cart_pole ( ( * i_analyzed_projections ).latp, ( * i_analyzed_projections ).lonp );
						analyzed_proj->setCartPole ( cart_pole );
					}

					//Add projection to the list of analyzed projections
					analysis_parameters.analyzed_projections.push_back ( analyzed_proj );
				}
			}
		}

		//Find meridians and parallels using RANSAC
		TMeridiansList <double> ::Type meridians;
		TParallelsList <double> ::Type parallels;

		if ( analysis_parameters.analysis_type.a_gn_tf )
		{
			//Perform RANSAC
			//Ransac::ransacFitMeridiansAndParallels ( nl_reference, meridians, parallels, MAX_RANSAC_ERROR_GRATICULE, 1.0, true, true, output );

			//Disable meridians/parallels analysis: not enough meridians / parallels have been found
			analysis_parameters.analysis_type.a_gn_tf = ( ( meridians.size() > 0 ) || ( parallels.size() > 0 ) );
		}

		//Precompute geometric structures for test dataset: Create Voronoi diagram and merge Voronoi cells
		Container <HalfEdge <double> *> hl_dt_test, hl_vor_test, hl_merge_test;
		Container <Node3DCartesian <double> *> nl_vor_test, intersections_test;
		Container <VoronoiCell <double> *> vor_cells_test;
		Container <Face <double> *> faces_test;

		//Precompute Voronoi diagram for the test dataset and merge Voronoi cells: only if outliers are not removed
		if ( analysis_parameters.analysis_type.a_vd_tf )
		{
			try
			{
				//Compute Voronoi cells with full topology
				Voronoi2D::VD ( ( Container <Node3DCartesian <double> *> & ) nl_test, nl_vor_test, hl_dt_test, hl_vor_test, vor_cells_test, AppropriateBoundedCells, TopologicApproach, 0, analysis_parameters.print_exceptions, output );

				//Are there enough Voronoi cells for analysis? If not, disable analysis
				analysis_parameters.analysis_type.a_vd_tf = ( faces_test.size() < MIN_BOUNDED_VORONOI_CELLS );

				//There is enough Voronoi cells for analysis
				if ( analysis_parameters.analysis_type.a_vd_tf )
				{
					//Merge Voronoi cells with all adjacent cells
					for ( TItemsList <Node3DCartesian <double> *> ::Type ::iterator i_points_test = nl_test.begin(); i_points_test != nl_test.end() ; i_points_test ++ )
					{
						//Get Voronoi cell
						VoronoiCell <double> *vor_cell_test = dynamic_cast < VoronoiCell <double> * > ( ( *i_points_test ) -> getFace() );

						//Merge Voronoi cell only if exists and it is bounded
						Face <double> * face_test = NULL;

						if ( ( vor_cell_test != NULL )  && ( vor_cell_test->getBounded() ) )
						{
							Voronoi2D::mergeVoronoiCellAndAdjacentCells ( vor_cell_test, &face_test, intersections_test, hl_merge_test );
						}

						//Add face to the list
						faces_test.push_back ( face_test );
					}
				}
			}

			//Throw exception
			catch ( Error & error )
			{
				//faces_test.clear();
				//nl_vor_test.clear();
				//hl_dt_test.clear();
				//hl_vor_test.clear();
				//vor_cells_test.clear();

				//Print exception
				if ( analysis_parameters.print_exceptions )
				{
					error.printException ( output );
				}
			}
		}

		//Compute cartometric analysis_type and create samples
		//Performed + thrown samples = total samples
		unsigned int total_created_or_thrown_samples = 0;
		Container <Sample <double> > sl;
		std::cout << std::endl << ">> Computing analysis & creating samples, please wait: " << std::endl;
		*output << std::endl << ">> Computing analysis & creating samples, please wait: " << std::endl;

		//Repeat analysis with the different sensitivity, if neccessary
		for ( unsigned int i = 0; ( i < analysis_parameters.analysis_repeat + 1 ) && ( analysis_parameters.heuristic_sensitivity_ratio < MAX_HEURISTIC_SENSTIVITY_RATIO ) &&
			( sl.size() == 0 ); i++, analysis_parameters.heuristic_sensitivity_ratio += 2 )
		{
			//Print analysis parameters
			std::cout << "Heuristic senstivity: " << analysis_parameters.heuristic_sensitivity_ratio << std::endl;
			*output << "Heuristic senstivity: " << analysis_parameters.heuristic_sensitivity_ratio << std::endl;

			//Perform analysis: simplex method
			if (analysis_parameters.analysis_method == SimplexMethod || analysis_parameters.analysis_method == SimplexRotMethod || analysis_parameters.analysis_method == SimplexRot2Method || analysis_parameters.analysis_method == SimplexShiftsMethod)
			{
				CartAnalysis::computeAnalysisForAllSamplesSim ( sl, proj_list, nl_test, nl_reference, meridians, parallels, faces_test, analysis_parameters,
					total_created_or_thrown_samples, output );
			}

			//Perform analysis: differential evolution
			else if (analysis_parameters.analysis_method == DifferentialEvolutionMethod || analysis_parameters.analysis_method == DifferentialEvolutionRotMethod || analysis_parameters.analysis_method == DifferentialEvolutionRot2Method || analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
			{
				CartAnalysis::computeAnalysisForAllSamplesDE ( sl, proj_list, nl_test, nl_reference, meridians, parallels, faces_test, analysis_parameters,
					total_created_or_thrown_samples, output );
			}

			//Perform analysis: minimum least squares
			else
			{
				CartAnalysis::computeAnalysisForAllSamplesMLS ( sl, proj_list, nl_test, nl_reference, meridians, parallels, faces_test, analysis_parameters,
					total_created_or_thrown_samples, output );
			}

			//Print total tested samples
			std::cout << std::endl << sl.size() << " sample(s) have been added for testing." << std::endl;
			*output << std::endl << sl.size() << " sample(s) have been added for testing." << std::endl;

			//Print thrown samples
			if ( total_created_or_thrown_samples > 0 )
			{
				std::cout  << total_created_or_thrown_samples - sl.size() << " sample(s) have been thrown [" << fixed << std::setprecision ( 2 ) << ( 100.0 * ( total_created_or_thrown_samples - sl.size() ) / ( total_created_or_thrown_samples ) ) << " %]."  << std::endl << std::endl;
				*output << total_created_or_thrown_samples - sl.size() << " sample(s) have been thrown [" << fixed << std::setprecision ( 2 ) << ( 100.0 * ( total_created_or_thrown_samples - sl.size() ) / ( total_created_or_thrown_samples ) ) << " %]."  << std::endl << std::endl;
			}
		}

		//No sample had been computed and added
		if ( sl.size() == 0 )
		{
			if ( analysis_parameters.analysis_repeat == 0 ) throw ErrorBadData ( "ErrorBadData: no sample had been computed. ", "Try to increase heuristic sensitivity ratio \"sens\", disable heuristic or check projection equations." );
			else throw ErrorBadData ( "ErrorBadData: no sample had been computed. ", "Disable heuristic or check projection equations." );
		}

		//Sort computed results
		output->flush();
		std::cout << ">> Sorting all samples...";
		CartAnalysis::sortSamplesByComputedRatios ( sl, analysis_parameters.analysis_type );
		std::cout << " Completed."  << std::endl << std::endl;

		//Print results into output
		//Enable after analysis
		std::cout << "Print results";
		//CartAnalysis::printResults(sl, nl_test, nl_reference, analysis_parameters, &std::cout);
		CartAnalysis::printResults2 ( sl, nl_test, nl_reference, analysis_parameters, output );

		//Compute solution diversity
		/*
		const double diversity = CartAnalysis::solutionDiversity ( sl, proj_list, nl_test, nl_reference, analysis_parameters.printed_results );
		std::cout << std::setprecision ( 3 ) << "Solution diversity (first " << analysis_parameters.printed_results << " candidates): "  << '\n' << diversity << " km \n\n";
		*output << std::setprecision ( 3 ) << "Solution diversity (first " << analysis_parameters.printed_results << " candidates): "  << '\n' << diversity << " km \n\n";
		*/

		//Create graticule for the results
		if ( analysis_parameters.exported_graticule > 0 )
		{
			std::cout << ">> Exporting points and graticules, please wait... \n";

			//Create graticules
			for ( unsigned int i = 0; i < std::min ( ( unsigned int ) analysis_parameters.exported_graticule, sl.size() ); i++ )
			{
				//Get a sample projection
				Projection <double> *proj = sl[i].getProj();

				//Set properties for a projection
				proj->setR ( sl[i].getR() );
				Point3DGeographic <double> cart_pole ( sl[i].getLatP(), sl[i].getLonP() );
				proj->setCartPole ( cart_pole );
				proj->setLat0 ( sl[i].getLat0() );
				//proj->setLat1(sl[i].getLat1());
				//proj->setLat2(sl[i].getLat2());
				proj->setLon0 ( sl[i].getLon0() );
				proj->setDx ( sl[i].getDx() );
				proj->setDy ( sl[i].getDy() );
				proj->setC ( sl[i].getC() );
				
				//Get limits
				TMinMax <double> lon_interval ( ( * std::min_element ( nl_reference.begin(), nl_reference.end(), sortPointsByLon () ) )->getLon(),
					( * std::max_element ( nl_reference.begin(), nl_reference.end(), sortPointsByLon () ) )->getLon() );
				TMinMax <double> lat_interval ( ( * std::min_element ( nl_reference.begin(), nl_reference.end(), sortPointsByLat () ) )->getLat(),
					( * std::max_element ( nl_reference.begin(), nl_reference.end(), sortPointsByLat () ) )->getLat() );

				//Stretch over the whole planishere
				if ((lon_interval.min_val < MIN_LON + 20) && (lon_interval.max_val > MAX_LON - 20))
				{
					lon_interval.min_val = MIN_LON;
					lon_interval.max_val = MAX_LON;
				}

				//std::cout << sl[i].getR() << "  " << sl[i].getDx() << "   " << sl[i].getDy() << "   " << sl[i].getAlpha() << '\n';
				//std::cout << lon_interval.min_val << "  " << lon_interval.max_val << "   " << lat_interval.min_val << "   " << lat_interval.max_val << "   " << analysis_parameters.lat_step << "   " << analysis_parameters.lon_step << '\n';
				//lat_interval.min_val=20; lat_interval.max_val=50;
				//lon_interval.min_val = 0; lon_interval.max_val = 50;
				//proj->setLon0(0);

				//Create data structures for the graticule representation
				unsigned int index = 0;
				TMeridiansListF <double> ::Type meridians_exp;
				TParallelsListF <double> ::Type parallels_exp;
				Container <Node3DCartesian <double> *> mer_par_points;

				//Set font height
				const double font_height = 0.05 * proj->getR() * std::min ( analysis_parameters.lat_step, analysis_parameters.lon_step ) * M_PI / 180;

				//Create graticule
				double alpha = sl[i].getAlpha();
				//alpha = 26;
				//std::cout << "grat = " << analysis_parameters.latp_step;
				Graticule::computeGraticule ( proj, lat_interval, lon_interval, analysis_parameters.lat_step, analysis_parameters.lon_step, 0.1 * analysis_parameters.lat_step, 0.1 * analysis_parameters.lon_step, alpha, TransformedGraticule, meridians_exp, parallels_exp, &mer_par_points, index );

				//Compute all points
				Container <Node3DCartesianProjected <double> *> nl_projected_temp;

				for (unsigned int j = 0; j < n_test_unique; j++)
				{
					//Get type of the direction
					TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

					//Reduce lon
					const double lon_red = CartTransformation::redLon0(nl_reference[j]->getLon(), proj->getLon0());

					double lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0, x_temp = 0.0, y_temp = 0.0;
					try
					{
						//Convert geographic point to oblique aspect
						lat_trans = CartTransformation::latToLatTrans(nl_reference[j]->getLat(), lon_red, cart_pole.getLat(), cart_pole.getLon());
						lon_trans = CartTransformation::lonToLonTrans(nl_reference[j]->getLat(), lon_red, cart_pole.getLat(), cart_pole.getLon(), trans_lon_dir);

						//Change coordinates of the point
						Point3DGeographic<double> point_geo(lat_trans, lon_trans);

						for (unsigned int k = 0; k < 3; k++)
						{
							try
							{
								//Compute equations
								x_temp = CartTransformation::latLonToX(&point_geo, proj, false);
								y_temp = CartTransformation::latLonToY(&point_geo, proj, false);

							}

							//2 attempt to avoid the singularity
							catch (Error &error)
							{
								lat_trans = point_geo.getLat();
								lon_trans = point_geo.getLon();

								//Move in latitude direction
								if (k == 0)
								{
									if (lat_trans == MAX_LAT) lat_trans -= 2*GRATICULE_ANGLE_SHIFT;
									else lat_trans += 2*GRATICULE_ANGLE_SHIFT;
								}

								//Move in longitude direction
								else if (k == 1)
								{
									if (lon_trans == MAX_LON) lon_trans -= 2*GRATICULE_ANGLE_SHIFT;
									else lon_trans += 2*GRATICULE_ANGLE_SHIFT;
								}

								//Neither first nor the second shhifts do not bring improvement
								else if (k == 2)
								{
									throw;
								}

								point_geo.setLat(lat_trans);
								point_geo.setLon(lon_trans);
							}
						}

						//Get shifts: need to be additionally subtracted
						const double dx = proj->getDx(), dy = proj->getDy();

						//Non-rotated projection
						if (alpha == 0.0)
						{
							x = x_temp; 
							y = y_temp;
						}

						//Rotated projection
						else
						{
							x = (x_temp - dx) * cos(alpha * M_PI / 180) - (y_temp - dy) * sin(alpha * M_PI / 180) + dx;
							y = (x_temp - dx) * sin(alpha * M_PI / 180) + (y_temp - dy) * cos(alpha * M_PI / 180) + dy;
						}
					}

					//Get error argument
					catch (ErrorMath <double> &error)
					{
						//Throw new math error
						//throw ErrorMathInvalidArgument <double>("ErrorMathInvalidArgument: error in coordinate functions (lat/lon).", "Can not compute exported points.", lat_trans);
					}

					//Create new cartographic point
					Node3DCartesianProjected <double> *n_projected = new Node3DCartesianProjected <double>(x, y);

					//Add point to the list
					nl_projected_temp.push_back(n_projected);
				}

				//nl_projected_temp.print();
				//nl_reference.print();
				
				//Convert position to char
				char posit_text [5];
				sprintf ( posit_text, "%d", i + 1 );

				//Create file name
				char output_file_graticule[MAX_TEXT_LENGTH], output_file_points_ref[MAX_TEXT_LENGTH], output_file_points_test[MAX_TEXT_LENGTH], output_file_proj4[MAX_TEXT_LENGTH];
				strcpy ( output_file_graticule, analysis_parameters.test_file );
				strcpy( output_file_points_test, output_file_graticule);  //Copy string without ID and projection name
				strcat ( output_file_graticule, "_" );
				
				strcat ( output_file_graticule, posit_text );
				strcat ( output_file_graticule, "_" );
				strcat ( output_file_graticule, proj->getProjectionName() );

				strcpy ( output_file_points_ref, output_file_graticule ); //Copy full string with ID and projection name
				strcpy (output_file_proj4, output_file_graticule); //Copy full string with ID and projection name
				strcat (output_file_graticule, "_grat.dxf");
				strcat (output_file_points_ref, "_points_ref.dxf");
				strcat (output_file_points_test, "_points_test.dxf");
				strcat(output_file_proj4, "_proj4.bat");

				//Export graticule into DXF file
				DXFExport::exportGraticuleToDXF ( output_file_graticule, meridians_exp, parallels_exp, mer_par_points, font_height, analysis_parameters.lat_step, analysis_parameters.lon_step );

				//Export estimated points into DXF file
				DXFExport::exportPointsToDXF(output_file_points_ref, nl_projected_temp, font_height);

				//Export test points into DXF file
				if (i == 0)
					DXFExport::exportPointsToDXF(output_file_points_test, nl_test, font_height);

				//Create Proj.4 definition string
				std::string proj4_string = ProjectionToProj4::ProjectionToProj4String(proj);
				std::ofstream output_file(output_file_proj4);
				output_file << proj4_string;
				output_file.close();

				//Print hash
				std::cout.flush();
				std::cout << ".";
			}
		}

		std::cout << " Completed."  << std::endl << std::endl;

		//Get time difference
		time ( &end );
		end_ascii = localtime ( &end );
		float time_difference_all = difftime ( end, start );

		//End
		std::cout << ">> End Analysis: " << asctime ( end_ascii );
		std::cout << std::setprecision ( 1 ) << ">> Processing takes " << time_difference_all << "sec." << std::endl;
		*output << ">> End Analysis: " << asctime ( end_ascii );
		*output << std::setprecision ( 1 ) << ">> Processing takes " << time_difference_all << "sec." << std::endl;
		std::cout << std::endl  << "Detection has been successfully completed..." << std::endl << std::endl;
		std::cout << "**************************************************************************" << std::endl << std::endl << std::endl;
		*output << "1";

		//Close file
		output_file.close();
	}

	//Throw exception
	catch ( Error & error )
	{
		std::cout << std::endl  << "Detection has not been successfully completed..." << std::endl;
		*output << std::endl  << "Detection has not been successfully completed..." << std::endl;
		error.printException();
		error.printException ( output );
		*output << "0";

		//Close file
		output_file.close();
	}

	//Memory leak analysis for VS 2010
	//_CrtDumpMemoryLeaks();
}


