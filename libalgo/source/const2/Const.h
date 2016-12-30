// Description: Definiton of the global constants

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

#ifndef CONST2_H
#define CONST2_H


//Global constants definition
//Please, do not chnage the values


#ifndef M_PI					//Pi definition
#define M_PI					3.141592653589793238462643383279502884
#endif

#ifndef RO					//Ro definition
#define RO					( 180.0 / M_PI )
#endif

#ifndef R0					//Initial Earth radius
#define R0					6380.0
#endif

#ifndef MAX_ANGULAR_DIFF			//Numeric threshold for the difference of two angles
#define MAX_ANGULAR_DIFF			1.0e-9
#endif

#ifndef BUFF					//Length of the text buffer
#define BUFF					1024
#endif

#ifndef MIN_FLOAT				//Minimum float value
#define MIN_FLOAT				1.0e-37
#endif

#ifndef MAX_FLOAT				//Maximum float value
#define MAX_FLOAT				1.0e37
#endif

#ifndef MAX_DOUBLE 				//Maximum double value
#define MAX_DOUBLE				1.7e+300
#endif

#ifndef EPS					//Maximum floating operation error: numeric threshold
#define EPS 					1.0e-15
#endif

#ifndef GRATICULE_LAT_LON_SHIFT			//Shift of the meridian/parallel point, when graticule constructed
#define GRATICULE_LAT_LON_SHIFT 		1.0e-3
#endif

#ifndef MAX_NLS_STEP_LENGTH			//Maximum length of the NLS step
#define MAX_NLS_STEP_LENGTH 			1.0e3
#endif

#ifndef MIN_LAT1				//Minimum latitude of the  true parallel
#define MIN_LAT1 				-85.0
#endif

#ifndef MAX_LAT1				//Minimum latitude of the true parallel
#define MAX_LAT1 				85.0
#endif

#ifndef MIN_LAT					//Minimum latitude
#define MIN_LAT					-90.0
#endif

#ifndef MAX_LAT					//Maximum latitude
#define MAX_LAT					90.0
#endif

#ifndef MIN_LON					//Minimum longitude
#define MIN_LON					-180.0
#endif

#ifndef MAX_LON					//Maximum longitude
#define MAX_LON					180.0
#endif

#ifndef MIN_POSITION_DIFF			//Minimum position difference between two 2D points. Points having angle difference < MIN_ANGLE_DIF are unique
#define MIN_POSITION_DIFF			5.0e-4
#endif

#ifndef ANGLE_ROUND_ERROR 			//Angle accuracy (cartographic distortion calculation)
#define ANGLE_ROUND_ERROR   			1.0e-8
#endif

#ifndef ARGUMENT_ROUND_ERROR			//Round error of the function argument (i.e. asin(arg))
#define ARGUMENT_ROUND_ERROR 			1.0e-5				
#endif

#ifndef NUM_DERIV_STEP				//Set step for numeric derivative using Stirling method
#define NUM_DERIV_STEP				0.001
#endif

#ifndef MAX_FLOAT_OPER_ERROR			//Maximum floating operation error
#define MAX_FLOAT_OPER_ERROR 			1.0e-15
#endif

#endif