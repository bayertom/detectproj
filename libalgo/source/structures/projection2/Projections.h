// Description: List of all implemented map projections
// Supported equations in the non-closed form,

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

#ifndef Projections_H
#define Projections_H

#include <memory>
#include <list>

#include "libalgo/source/types/TListS.h"

#include "ProjectionAzimuthal.h"
#include "ProjectionConic.h"
#include "ProjectionCylindrical.h"
#include "ProjectionEllipsoidal.h"
#include "ProjectionMiscellaneous.h"
#include "ProjectionPolyConic.h"
#include "ProjectionPseudoAzimuthal.h"
#include "ProjectionPseudoConic.h"
#include "ProjectionPseudoCylindrical.h"


//Definition of all projections
class Projections
{
public:

	template <typename T>
	static void init(TListS <Projection<T> > &projections);		//Initialize all projections

	//*******************List of all supported projection equations  **************************************
	template <typename T>
	static T X_def(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_def(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	//Albers equal area
	template <typename T>
	static T X_aea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_aea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Azimuthal Equidistant
	template <typename T>
	static T X_aeqd(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_aeqd(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Aitoff
	template <typename T>
	static T X_aitoff(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_aitoff(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Apianus
	template <typename T>
	static T X_api(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_api(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Apianus elliptical
	template <typename T>
	static T X_apiel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_apiel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Armadillo
	template <typename T>
	static T X_armad(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_armad(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//August Epicycloidal
	template <typename T>
	static T X_august(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_august(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Bacon Globular
	template <typename T>
	static T X_bacon(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_bacon(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Boggs Eumorphic
	template <typename T>
	static T X_boggs(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_boggs(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Bonne
	template <typename T>
	static T X_bonne(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_bonne(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Breusign
	template <typename T>
	static T X_breus(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_breus(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Central Cylindrical
	template <typename T>
	static T X_cc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_cc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Cylindrical Equal Area
	template <typename T>
	static T X_cea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_cea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Clark Perspective
	template <typename T>
	static T X_clar(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_clar(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Collignon
	template <typename T>
	static T X_collg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_collg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Craster Parabolic
	template <typename T>
	static T X_crast(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_crast(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Denoyer Semi-Elliptical
	template <typename T>
	static T X_denoy(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_denoy(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Eckert I
	template <typename T>
	static T X_eck1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Eckert II
	template <typename T>
	static T X_eck2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Eckert III
	template <typename T>
	static T X_eck3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Eckert IV
	template <typename T>
	static T X_eck4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_eck4(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_eck4(const T lat, const T theta);

	
	//Eckert V
	template <typename T>
	static T X_eck5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Eckert VI
	template <typename T>
	static T X_eck6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eck6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_eck6(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_eck6(const T lat, const T theta);

	
	//Eisenlohr
	template <typename T>
	static T X_eisen(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eisen(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Equidistant Cylindrical
	template <typename T>
	static T X_eqc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eqc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Equidistant Conic (true parallel lat1)
	template <typename T>
	static T X_eqdc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eqdc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Equidistant Conic (true parallels lat1, lat2)
	template <typename T>
	static T X_eqdc2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eqdc2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Equidistant Conic (true parallel lat1, pole = point)
	template <typename T>
	static T X_eqdc3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_eqdc3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Fahey
	template <typename T>
	static T X_fahey(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_fahey(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Foucaut Sine-Tangent
	template <typename T>
	static T X_fouc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_fouc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Foucaut Sinusoidal
	template <typename T>
	static T X_fouc_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_fouc_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Fournier I Globular
	template <typename T>
	static T X_fourn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_fourn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Fournier II Elliptical
	template <typename T>
	static T X_fourn2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_fourn2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Gall Stereographic
	template <typename T>
	static T X_gall(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_gall(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Ginsburg VIII (TsNIIGAiK)
	template <typename T>
	static T X_gins8(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_gins8(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Gnomonic
	template <typename T>
	static T X_gnom(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_gnom(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Goode Homolosine
	template <typename T>
	static T X_goode(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_goode(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Hammer
	template <typename T>
	static T X_hammer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_hammer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Hatano Asymmetrical Equal Area
	template <typename T>
	static T X_hataea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_hataea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_hataea(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_hataea(const T lat, const T theta);


	//La Hire
	template <typename T>
	static T X_hire(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_hire(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Interrupted Goode Homolosine
	template <typename T>
	static T X_igh(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_igh(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//James Perspective
	template <typename T>
	static T X_jam(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_jam(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Kavrayskiy V
	template <typename T>
	static T X_kav5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_kav5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Kavrayskiy VII
	template <typename T>
	static T X_kav7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_kav7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	//template <typename T>
	//T X_krovak(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	//template <typename T>
	//T Y_krovak(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Lambert Azimuthal Equal Area
	template <typename T>
	static T X_laea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_laea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Lagrange Conformal
	template <typename T>
	static T X_lagrng(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_lagrng(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Larrivee
	template <typename T>
	static T X_larr(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_larr(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Laskowski
	template <typename T>
	static T X_lask(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_lask(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Lambert Conformal Conic
	template <typename T>
	static T X_lcc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_lcc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Lambert Equal Area Conic (standard parallel lat1)
	template <typename T>
	static T X_leac(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_leac(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Lambert Equal Area Conic (standard parallels lat1, lat2)
	template <typename T>
	static T X_leac2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_leac2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Loximuthal
	template <typename T>
	static T X_loxim(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_loxim(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//McBryde-Thomas Sine I.
	template <typename T>
	static T X_mbt_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_mbt_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//McBryde-Thomas Flat-Pole Sine III
	template <typename T>
	static T X_mbt_s3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_mbt_s3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//McBryde-Thomas Flat-Pole Quartic IV
	template <typename T>
	static T X_mbtfpq(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_mbtfpq(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_mbtfpq(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_mbtfpq(const T lat, const T theta);

	
	//McBryde-Thomas Flat-Pole Sine II
	template <typename T>
	static T X_mbtfps(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_mbtfps(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_mbtfps(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_mbtfps(const T lat, const T theta);


	//Mercator
	template <typename T>
	static T X_merc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_merc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Miller Cylindrical
	template <typename T>
	static T X_mill(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_mill(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Mollweide
	template <typename T>
	static T X_moll(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_moll(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_moll(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_moll(const T lat, const T theta);


	//Nell
	template <typename T>
	static T X_nell(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_nell(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_nell(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_nell(const T lat, const T theta);


	//Nell-Hammer
	template <typename T>
	static T X_nell_h(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_nell_h(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Nicolosi Globular
	template <typename T>
	static T X_nicol(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_nicol(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Ortelius
	template <typename T>
	static T X_ortel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_ortel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Orthographic
	template <typename T>
	static T X_ortho(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_ortho(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Parabolic
	template <typename T>
	static T X_parab(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_parab(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Perspective azimuthal
	template <typename T>
	static T X_pers(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_pers(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);
	

	//Far-side Perspective Azimuthal
	template <typename T>
	static T X_persf(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_persf(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Near-side Perspective Azimuthal
	template <typename T>
	static T X_persn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_persn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Hassler Polyconic, American
	template <typename T>
	static T X_poly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_poly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Putnins P1
	template <typename T>
	static T X_putp1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Putnins P2
	template <typename T>
	static T X_putp2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_putp2(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_putp2(const T lat, const T theta);


	//Putnins P3
	template <typename T>
	static T X_putp3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Putnins P3P
	template <typename T>
	static T X_putp3p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp3p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Putnins P4P
	template <typename T>
	static T X_putp4p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp4p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Putnins P5
	template <typename T>
	static T X_putp5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Putnins P5P
	template <typename T>
	static T X_putp5p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp5p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Putnins P6
	template <typename T>
	static T X_putp6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_putp6(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_putp6(const T lat, const T theta);


	//Putnins P6P
	template <typename T>
	static T X_putp6p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_putp6p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_putp6p(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_putp6p(const T lat, const T theta);


	//Quartic Authalic
	template <typename T>
	static T X_qua_aut(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_qua_aut(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Rectangular Polyconic
	template <typename T>
	static T X_rpoly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_rpoly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Sinusoidal
	template <typename T>
	static T X_sinu(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_sinu(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Solovyev Azimuthal
	template <typename T>
	static T X_solo(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_solo(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Stereographic
	template <typename T>
	static T X_stere(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_stere(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Twilight General Vertical Perspective
	template <typename T>
	static T X_twi(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_twi(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Urmaev V
	template <typename T>
	static T X_urm5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_urm5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Van der Grinten I
	template <typename T>
	static T X_vandg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_vandg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Van der Grinten II
	template <typename T>
	static T X_vandg2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_vandg2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Van der Grinten III
	template <typename T>
	static T X_vandg3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);
	
	template <typename T>
	static T Y_vandg3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Van der Grinten IV
	template <typename T>
	static T X_vandg4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);
	
	template <typename T>
	static T Y_vandg4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);
	

	//Wagner I
	template <typename T>
	static T X_wag1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Wagner II
	template <typename T>
	static T X_wag2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Wagner III
	template <typename T>
	static T X_wag3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	
	//Wagner IV
	template <typename T>
	static T X_wag4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_wag4(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_wag4(const T lat, const T theta);


	//Wagner VI
	template <typename T>
	static T X_wag6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Wagner VII
	template <typename T>
	static T X_wag7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wag7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Werner-Staab
	template <typename T>
	static T X_wer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Werenskiold I
	template <typename T>
	static T X_weren(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_weren(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	///Winkel I
	template <typename T>
	static T X_wink1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wink1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);


	//Winkel II
	template <typename T>
	static T X_wink2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wink2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T FTheta_wink2(const T lat, const T theta);

	template <typename T>
	static T FThetaDer_wink2(const T lat, const T theta);


	//Winkel Tripel
	template <typename T>
	static T X_wintri(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T>
	static T Y_wintri(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c);

	template <typename T, typename FTheta, typename FThetaDer>
	static T NewtonRaphson(FTheta ftheta, FThetaDer fthetader, const T lat, const T theta0, const unsigned int MAX_ITERATIONS = 20, const T MAX_DIFF = 1.0e-5);

	template <typename T>
	static int sign(T value)
	{
		return (T(0) < value) - (value < T(0));
	}
};

//Constructor of the class Projection
template<typename T>
Projection <T>::Projection() : R(1.0), lon0(0.0), dx(0.0), dy(0.0), c(0.5), X(&Projections::X_def), Y(&Projections::Y_def), name("default") {}

#include "Projections.hpp"

#endif