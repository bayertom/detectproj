// Description: Definiton of the global constants

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


//Global constants definition
//Please, do not chnage the values


#ifndef M_PI					//Pi definition
#define M_PI					3.141592653589793238462643383279502884
#endif

#ifndef RO					//Ro definition
#define RO					( 180.0 / M_PI )
#endif

#ifndef ANGLE_ROUND_ERROR 							//Angle accuracy (cartographic distortion calculation)
#define ANGLE_ROUND_ERROR   			1.0e-8
#endif

#ifndef ARGUMENT_ROUND_ERROR
#define ARGUMENT_ROUND_ERROR 			1.0e-5				//Round error of the function argument (i.e. asin(arg))
#endif

#ifndef POSITION_ROUND_ERROR							//Round error of the point position to the line
#define POSITION_ROUND_ERROR 			1.0e-7
#endif

#ifndef MAX_FLOAT_OPER_ERROR							//Maximum floating operation error
#define MAX_FLOAT_OPER_ERROR 			1.0e-15
#endif

#ifndef GRATICULE_ANGLE_SHIFT							//Angle shift used to avoid singular points when a graticule is computed
#define GRATICULE_ANGLE_SHIFT 			1.0e-3
#endif

#ifndef AREA_ROUND_ERROR							//Round error of the point position to the line
#define AREA_ROUND_ERROR 			1.0e-5
#endif

#ifndef MAX_EDGES								//Maximum Delaunay/Voronoi half edges stored in vector
#define MAX_EDGES				10000000
#endif

#ifndef BUFF									//Length of the text buffer
#define BUFF					1024
#endif

#ifndef MIN_ANGLE_DIFF								//Minimum angle difference between 2 cartographic points. Points having angle difference < MIN_ANGLE_DIF are unique
#define MIN_ANGLE_DIFF				1.0e-6
#endif

#ifndef MIN_POSITION_DIFF							//Minimum position difference between two 2D points. Points having angle difference < MIN_ANGLE_DIF are unique
#define MIN_POSITION_DIFF			5.0e-4
#endif

#ifndef MAX_INT 								//Maximum int value
#define MAX_INT 				2147483647
#endif

#ifndef MAX_DOUBLE 								//Maximum double value
#define MAX_DOUBLE				1.7e+300
#endif

#ifndef MIN_DOUBLE								//Minimum double value
#define MIN_DOUBLE				1.7e-300
#endif

#ifndef MAX_FLOAT								//Maximum float value
#define MAX_FLOAT				1.0e37
#endif

#ifndef MIN_FLOAT								//Minimum float value
#define MIN_FLOAT				1.0e-37
#endif

#ifndef MAX_FLOAT_EXPONENT							//Maximum exponent for the floating value
#define MAX_FLOAT_EXPONENT 			log(MAX_FLOAT)
#endif

#ifndef MAX_NODES								//Maximum of nodes stored in the list
#define MAX_NODES				5000000
#endif

#ifndef MAX_POINT_COORDINATE							//Maximum value of the cartographic point coordinate
#define MAX_POINT_COORDINATE			1.0e+9
#endif

#ifndef MAX_ROWS_MATRIX								//Maximum rows of the matrix
#define MAX_ROWS_MATRIX				1000
#endif

#ifndef MAX_COLUMNS_MATRIX							//Maximum cols of the matrix
#define MAX_COLUMNS_MATRIX			1000
#endif


#ifndef MAX_TEXT_LENGTH								//Maximum length of the char
#define MAX_TEXT_LENGTH				4096
#endif

#ifndef MAX_TEXT_FILE_LENGTH							//Maximum length of the text file
#define MAX_TEXT_FILE_LENGTH			800000
#endif

#ifndef MAX_SAMPLES								//Maximum cartographic samples stored in vector
#define MAX_SAMPLES				1000000
#endif

#ifndef MIN_BOUNDED_VORONOI_CELLS						//Minimum amount of unbounded Voronoi cells usable for analysis
#define MIN_BOUNDED_VORONOI_CELLS		3
#endif


#ifndef MIN_POINTS								//Minimum amount of points usable for analysis
#define MIN_POINTS				7
#endif

#ifndef MIN_ANALYSIS_REPEAT
#define MIN_ANALYSIS_REPEAT			0				//Minimum steps of repetiotions of cartometric analysis
#endif

#ifndef MAX_ANALYSIS_REPEAT
#define MAX_ANALYSIS_REPEAT			5				//Maximum steps of repetiotions of cartometric analysis
#endif

#ifndef MIN_HEURISTIC_SENSTIVITY_RATIO
#define MIN_HEURISTIC_SENSTIVITY_RATIO		0.1				//Set minimum heuristic sensitivity ratio for throwing non-perspecticve samples
#endif

#ifndef MAX_HEURISTIC_SENSTIVITY_RATIO
#define MAX_HEURISTIC_SENSTIVITY_RATIO		15.0				//Set maximum heuristic sensitivity ratio for throwing non-perspecticve samples
#endif

#ifndef MIN_HEURISTIC_SENSTIVITY_INCREMENT
#define MIN_HEURISTIC_SENSTIVITY_INCREMENT	0.2				//Set minimum heuristic sensitivity increasement
#endif

#ifndef MAX_HEURISTIC_SENSTIVITY_ICREMENT
#define MAX_HEURISTIC_SENSTIVITY_ICREMENT	5.0				//Set maximum heuristic sensitivity increasement
#endif

#ifndef RANSAC_MIN_LINE_POINTS
#define RANSAC_MIN_LINE_POINTS			4				//Minimum points on the line detected by the RANSAC
#endif

#ifndef NUM_DERIV_STEP
#define NUM_DERIV_STEP				0.001				//Set step for numeric derivative using Stirling method
#endif

#ifndef MAX_RANSAC_ERROR_GRATICULE
#define MAX_RANSAC_ERROR_GRATICULE		0.05				//Maximum acceptable error in meridian/parallels points detected by RANSAC
#endif

#ifndef MIN_LAT
#define MIN_LAT					-90.0				//Minimum latitude
#endif

#ifndef MAX_LAT
#define MAX_LAT					90.0				//Maximum latitude
#endif

#ifndef MIN_LON
#define MIN_LON					-180.0				//Minimum longitude
#endif

#ifndef MAX_LON
#define MAX_LON					180.0				//Maximum longitude
#endif

#ifndef MIN_LAT0
#define MIN_LAT0				10.0				//Minimum latitude of the undistorted parallel
#endif

#ifndef MAX_LAT0
#define MAX_LAT0				80.0				//Maximum latitude of the undistorted parallel
#endif

#ifndef MAX_C
#define MAX_C					10000000.0			//Maximum constant of the projection
#endif

#ifndef MIN_LAT_STEP
#define MIN_LAT_STEP				0.01				//Minimum latitude offset of two parallels
#endif

#ifndef MAX_LAT_STEP
#define MAX_LAT_STEP				50.0				//Maximum latitude offset of two parallels
#endif

#ifndef MIN_LON_STEP
#define MIN_LON_STEP				0.01				//Minimum longitude offset of two meridians
#endif

#ifndef MAX_LON_STEP
#define MAX_LON_STEP				50.0				//Maximum longitude offset of two meridians
#endif

#ifndef LAT_LIMIT_DIST
#define LAT_LIMIT_DIST				80.0				//Latitude limit for comutation of cartographic distortions
#endif

#ifndef IMPROVE_RATIO_STD_DEV
#define IMPROVE_RATIO_STD_DEV			10.0
#endif

#ifndef REM_DIV_ROT_ANGLE
#define REM_DIV_ROT_ANGLE			2
#endif

#ifndef TURNING_FUNCTION_MAX_DIFFERENCE
#define TURNING_FUNCTION_MAX_DIFFERENCE		1.0
#endif

#ifndef MATCHING_FACTOR
#define MATCHING_FACTOR				0.75
#endif

#ifndef EPS							//Maximum floating operation error
#define EPS 			1.0e-15
#endif
