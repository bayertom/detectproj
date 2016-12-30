// Description: List of all implemented map projections
// Supported equations in the non-closed form,

// Copyright (c) 2015 - 2016
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3.0 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef Projections_HPP
#define Projections_HPP

#include <cmath>

#include "libalgo/source/exceptions/MathOverflowException.h"
#include "libalgo/source/exceptions/MathZeroDevisionException.h"

template <typename T>
void Projections::init(TListS <Projection<T> > &projections)
{
	// Initialize all supported projections and add them to the list
	auto aea = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_aea<T>, Y_aea<T>, "aea");
	auto aeqd = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_aeqd<T>, Y_aeqd<T>, "aeqd");
	auto aitoff = std::make_shared <ProjectionPseudoAzimuthal < T >> (R0, 90.0, 0.0, 10, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_aitoff<T>, Y_aitoff<T>, "aitoff");
	auto api = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_api<T>, Y_api<T>, "api");
	auto apiel = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_apiel<T>, Y_apiel<T>, "apiel");
	auto armad = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_armad<T>, Y_armad<T>, "armad");
	auto august = std::make_shared <ProjectionMiscellaneous < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_august<T>, Y_august<T>, "august");
	auto bacon = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_bacon<T>, Y_bacon<T>, "bacon");
	auto boggs = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_boggs<T>, Y_boggs<T>, "boggs");
	auto bonne = std::make_shared <ProjectionPseudoConic < T >> (R0, 90.0, 0.0, 40.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_bonne<T>, Y_bonne<T>, "bonne");
	
	auto breus = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_breus<T>, Y_breus<T>, "breus");
	auto cc = std::make_shared <ProjectionCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_cc<T>, Y_cc<T>, "cc");
	auto cea = std::make_shared <ProjectionCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_cea<T>, Y_cea<T>, "cea");
	auto clar = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_clar<T>, Y_clar<T>, "breus");
	auto collg = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection, 0.0, 0.0, 0.0, 1.0, X_collg<T>, Y_collg<T>, "collg");
	auto crast = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_crast<T>, Y_crast<T>, "crast");
	auto denoy = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_denoy<T>, Y_denoy<T>, "denoy");
	auto eck1 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck1<T>, Y_eck1<T>, "eck1");
	auto eck2 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck2<T>, Y_eck2<T>, "eck2");
	auto eck3 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck3<T>, Y_eck3<T>, "eck3");
	
	auto eck4 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck4<T>, Y_eck4<T>, "eck4");
	auto eck5 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck5<T>, Y_eck5<T>, "eck5");
	auto eck6 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eck6<T>, Y_eck6<T>, "eck6");
	auto eisen = std::make_shared <ProjectionMiscellaneous < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eisen<T>, Y_eisen<T>, "eisen");
	auto eqc = std::make_shared <ProjectionCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eqc<T>, Y_eqc<T>, "eqc");
	auto eqdc = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eqdc<T>, Y_eqdc<T>, "eqdc");
	auto eqdc2 = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eqdc2<T>, Y_eqdc2<T>, "eqdc2");
	auto eqdc3 = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_eqdc3<T>, Y_eqdc3<T>, "eqdc3");
	auto fahey = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_fahey<T>, Y_fahey<T>, "fahey");
	auto fouc = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_fouc<T>, Y_fouc<T>, "fouc");
	
	auto fouc_s = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_fouc_s<T>, Y_fouc_s<T>, "fouc_s");
	auto fourn = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_fourn<T>, Y_fourn<T>, "fourn");
	auto fourn2 = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_fourn2<T>, Y_fourn2<T>, "fourn2");
	auto gall = std::make_shared <ProjectionCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_gall<T>, Y_gall<T>, "gall");
	auto gins8 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_gins8<T>, Y_gins8<T>, "gins8");
	auto gnom = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_gnom<T>, Y_gnom<T>, "gnom");
	auto goode = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_goode<T>, Y_goode<T>, "goode");
	auto hammer = std::make_shared <ProjectionPseudoAzimuthal < T >> (R0, 90.0, 0.0, 10, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_hammer<T>, Y_hammer<T>, "hammer");
	auto hataea = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_hataea<T>, Y_hataea<T>, "hataea");
	auto hire = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_hire<T>, Y_hire<T>, "hire");
	
	auto igh = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_igh<T>, Y_igh<T>, "igh");
	auto jam = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_jam<T>, Y_jam<T>, "jam");
	auto kav5 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_kav5<T>, Y_kav5<T>, "kav5");
	auto kav7 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_kav7<T>, Y_kav7<T>, "kav7");
	auto laea = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_laea<T>, Y_laea<T>, "laea");
	auto lagrng = std::make_shared <ProjectionMiscellaneous < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_lagrng<T>, Y_lagrng <T>, "lagrng");
	auto larr = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_larr<T>, Y_larr <T>, "larr");
	auto lask = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_lask<T>, Y_lask<T>, "lask");
	auto lcc = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_lcc<T>, Y_lcc<T>, "lcc");
	auto leac = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_leac<T>, Y_leac<T>, "leac");
	
	auto leac2 = std::make_shared <ProjectionConic < T >> (R0, 90.0, 0.0, 40.0, 50.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_leac2<T>, Y_leac2<T>, "leac2");
	auto loxim = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 10.0, 0.0, 0.0, 1.0, X_loxim<T>, Y_loxim<T>, "loxim");
	auto mbt_s = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_mbt_s<T>, Y_mbt_s<T>, "mbt_s");
	auto mbt_s3 = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_mbt_s3<T>, Y_mbt_s3<T>, "mbt_s3");
	auto mbtfpq = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_mbtfpq<T>, Y_mbtfpq<T>, "mbtfpq");
	auto mbtfps = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_mbtfps<T>, Y_mbtfps<T>, "mbtfps");
	auto merc = std::make_shared <ProjectionCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_merc<T>, Y_merc<T>, "merc");
	auto mill = std::make_shared <ProjectionCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_mill<T>, Y_mill<T>, "mill");
	auto moll = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_moll<T>, Y_moll<T>, "moll");
	auto nell = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_nell<T>, Y_nell<T>, "nell");
	
	auto nell_h = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_nell_h<T>, Y_nell_h<T>, "nell_h");
	auto nicol = std::make_shared <ProjectionMiscellaneous < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_nicol<T>, Y_nicol<T>, "nicol");
	auto ortel = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_ortel<T>, Y_ortel<T>, "ortel");
	auto ortho = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_ortho<T>, Y_ortho<T>, "ortho");
	auto parab = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_parab<T>, Y_parab<T>, "parab");
	auto pers = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_pers<T>, Y_pers<T>, "pers");
	auto persf = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_persf<T>, Y_persf<T>, "persf");
	auto persn = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_persn<T>, Y_persn<T>, "persn");
	auto poly = std::make_shared <ProjectionPolyConic < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_poly<T>, Y_poly<T>, "poly");
	auto putp1 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp1<T>, Y_putp1<T>, "putp1");
	
	auto putp2 = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp2<T>, Y_putp2<T>, "putp2");
	auto putp3 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp3<T>, Y_putp3<T>, "putp3");
	auto putp3p = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp3p<T>, Y_putp3p<T>, "putp3p");
	auto putp4p = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp4p<T>, Y_putp4p<T>, "putp4p");
	auto putp5 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp5<T>, Y_putp5<T>, "putp5");
	auto putp5p = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp5p<T>, Y_putp5p<T>, "putp5p");
	auto putp6 = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp6<T>, Y_putp6<T>, "putp6");
	auto putp6p = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_putp6p<T>, Y_putp6p<T>, "putp6p");
	auto qua_aut = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_qua_aut<T>, Y_qua_aut<T>, "qua_aut");
	auto rpoly = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_rpoly<T>, Y_rpoly<T>, "rpoly");
	
	auto sinu = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_sinu<T>, Y_sinu<T>, "sinu");
	auto solo = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_solo<T>, Y_solo<T>, "solo");
	auto stere = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_stere<T>, Y_stere<T>, "stere");
	auto twi = std::make_shared <ProjectionAzimuthal < T >> (R0, 90.0, 0.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_twi<T>, Y_twi<T>, "twi");
	auto urm5 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_urm5<T>, Y_urm5<T>, "urm5");
	auto vandg = std::make_shared <ProjectionPolyConic < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_vandg<T>, Y_vandg<T>, "vandg");
	auto vandg2 = std::make_shared <ProjectionPolyConic < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_vandg2<T>, Y_vandg2<T>, "vandg2");
	auto vandg3 = std::make_shared <ProjectionPolyConic < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_vandg3<T>, Y_vandg3<T>, "vandg3");
	auto vandg4 = std::make_shared <ProjectionPolyConic < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_vandg4<T>, Y_vandg4<T>, "vandg4");
	auto wag1 = std::make_shared <ProjectionPseudoCylindrical < T >>(R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag1<T>, Y_wag1<T>, "wag1");
	
	auto wag2 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag2<T>, Y_wag2<T>, "wag2");
	auto wag3 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag3<T>, Y_wag3<T>, "wag3");
	auto wag4 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag4<T>, Y_wag4<T>, "wag4");
	auto wag6 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag6<T>, Y_wag6<T>, "wag6");
	auto wag7 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wag7<T>, Y_wag7<T>, "wag7");
	auto wer = std::make_shared <ProjectionPseudoAzimuthal < T >> (R0, 90.0, 0.0, 10, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wer<T>, Y_wer<T>, "wer");
	auto weren = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_weren<T>, Y_weren<T>, "weren");
	auto wink1 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wink1<T>, Y_wink1<T>, "wink1");
	auto wink2 = std::make_shared <ProjectionPseudoCylindrical < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wink2<T>, Y_wink2<T>, "wink2");
	auto wintri = std::make_shared <ProjectionPseudoAzimuthal < T >> (R0, 90.0, 0.0, 10.0, NormalDirection2, 0.0, 0.0, 0.0, 1.0, X_wintri<T>, Y_wintri<T>, "wintri");
	
	//Add projections to the list
	
	projections.push_back(aea);
	projections.push_back(aeqd);
	projections.push_back(aitoff);
	projections.push_back(api);
	projections.push_back(apiel);
	projections.push_back(armad);
	projections.push_back(august);
	projections.push_back(bacon);
	projections.push_back(boggs);
	projections.push_back(bonne);
	
	projections.push_back(breus);
	projections.push_back(cc);
	projections.push_back(cea);
	projections.push_back(clar);
	projections.push_back(collg);
	projections.push_back(crast);
	projections.push_back(denoy);
	projections.push_back(eck1);
	projections.push_back(eck2);
	projections.push_back(eck3);

	projections.push_back(eck4);
	projections.push_back(eck5);
	projections.push_back(eck6);
	projections.push_back(eisen);
	projections.push_back(eqc);
	projections.push_back(eqdc);
	projections.push_back(eqdc2);
	projections.push_back(eqdc3);
	projections.push_back(fahey);
	projections.push_back(fouc);

	projections.push_back(fouc_s);
	projections.push_back(fourn);
	projections.push_back(fourn2);
	projections.push_back(gall);
	projections.push_back(gins8);
	
	projections.push_back(gnom);
	
	projections.push_back(goode);
	projections.push_back(hammer);
	projections.push_back(hataea);
	projections.push_back(hire);

	//projections.push_back(igh);
	projections.push_back(jam);
	projections.push_back(kav5);
	projections.push_back(kav7);
	projections.push_back(laea);
	projections.push_back(lagrng);
	projections.push_back(larr);
	projections.push_back(lask);
	projections.push_back(lcc);
	projections.push_back(leac);
	
	projections.push_back(leac2);
	projections.push_back(loxim);
	projections.push_back(mbt_s);
	projections.push_back(mbt_s3);
	projections.push_back(mbtfpq);
	projections.push_back(mbtfps);
	projections.push_back(merc);
	projections.push_back(mill);
	projections.push_back(moll);
	projections.push_back(nell);
	
	projections.push_back(nell_h);
	projections.push_back(nicol);
	projections.push_back(ortel);
	projections.push_back(ortho);
	projections.push_back(parab);
	projections.push_back(pers);
	projections.push_back(persf);
	projections.push_back(persn);
	projections.push_back(poly);
	projections.push_back(putp1);
	
	projections.push_back(putp2);
	projections.push_back(putp3);
	projections.push_back(putp3p);
	projections.push_back(putp4p);
	projections.push_back(putp5);
	projections.push_back(putp5p);
	projections.push_back(putp6);
	projections.push_back(putp6p);
	projections.push_back(qua_aut);
	projections.push_back(rpoly);
	
	projections.push_back(sinu);
	
	projections.push_back(solo);
	projections.push_back(stere);
	projections.push_back(twi);
	projections.push_back(urm5);
	projections.push_back(vandg);
	projections.push_back(vandg2);
	projections.push_back(vandg3);
	projections.push_back(vandg4);
	projections.push_back(wag1);
	
	projections.push_back(wag2);
	projections.push_back(wag3);
	projections.push_back(wag4);
	projections.push_back(wag6);
	projections.push_back(wag7);
	projections.push_back(wer);
	projections.push_back(weren);
	projections.push_back(wink1);
	projections.push_back(wink2);
	projections.push_back(wintri);
	
}

//Default projection
template <typename T>
T Projections::X_def(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	return lat / RO;
}

template <typename T>
T Projections::Y_def(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	return lon / RO;
}


//Albers equal area
template <typename T>
T Projections::X_aea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{ 
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = sin(lat1 / RO) + sin(lat2 / RO);
	const T B = cos(lat1 / RO) * cos(lat1 / RO);
	const T C = B + A * (sin(lat1 / RO) - sin(lat / RO));
	
	//Throw exception
	if (C < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(C) in X_aea coordinate function, ", "C < 0: ", C);
	
	const T X = 2.0 * R / (A) * sqrt(C) * sin(A / 2.0 * lonr / RO) + dx;

	//Throw exception
	if (fabs (X) > MAX_FLOAT )  
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_aea coordinate function, ", "X_aea > MAX_FLOAT: ", X);
	
	return X;
}

template <typename T>
T Projections::Y_aea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T A = sin(lat1 / RO) + sin(lat2 / RO);
	const T B = cos(lat1 / RO) * cos(lat1 / RO);
	const T C = B + A * (sin(lat1 / RO) - sin((lat1 + lat2) / 2.0 / RO));
	const T D = B + A * (sin(lat1 / RO) - sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT  || C <= 0)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_aea coordinate function, ", "1.0 / A, A = 0:", A);

	//Throw exception
	if (D < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(D) in Y_aea coordinate function, ", "D < 0: ", D);

	const T Y = 2.0 * R / A * sqrt(C) - 2.0 * R / A * sqrt(D) * cos(A / 2.0 * lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_aea coordinate function, ", "Y_AEA > MAX_FLOAT: ", Y);

	return Y;
}


//Azimuthal equidistant
template <typename T>
T Projections::X_aeqd(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * (90 - lat) / RO * sin(lonr/RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_aeqd coordinate function, ", "X_aeqd > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_aeqd(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = -R * (90 - lat) / RO * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_aeqd coordinate function, ", "Y_aeqd > MAX_FLOAT: ", Y);

	return Y;
}


//Aitoff
template <typename T>
T Projections::X_aitoff(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	T X;
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat / RO) * cos(lonr / 2.0 / RO);

	//Throw exception
	if (fabs(A) > 1)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate acos(A) in X_aitoff coordinate function, ", "A > 1: ", A);
	
	const T B = acos(A);

	if (fabs(B) < MIN_FLOAT)
		X = dx;

	else
	{
		const T C = sin(lat / RO) / sin(B);
		const T D = 1.0 - C * C;
		
		if (D < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(D) in X_aitoff coordinate function, ", "D < 0: ", D);


		X = 2.0 * R * B * sign(lonr) * sqrt(D) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_aitoff coordinate function, ", "X_aitoff > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_aitoff(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	T Y;

	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat / RO) * cos(lonr / 2.0 / RO);

	//Throw exception
	if (A > 1)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate acos(A) in Y_aitoff coordinate function, ", "A > 1: ", A);

	const T B = acos(A);

	if (fabs(B) < MIN_FLOAT)
		Y = dy;

	else
	{
		const T C = sin(lat / RO) / sin(B);

		Y = R * B * C + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_aitoff coordinate function, ", "Y_aitoff > MAX_FLOAT: ", Y);

	return Y;
}


//Apianus
template <typename T>
T Projections::X_api(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs(lonr) < MAX_ANGULAR_DIFF)
	{
		X = dx;
	}

	else
	{
		//Analogous to Bacon globular, but different Y
		const T F = ((M_PI / 2) *(M_PI / 2) * RO / abs(lonr) + abs(lonr) / RO) / 2;
		const T Y = R * lat / RO;
		const T G = F * F - Y * Y / (R * R);

		//Throw exception
		if (G < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(G) in X_bacon coordinate function, ", "G < 0: ", G);

		X = R * sign(lonr) * (abs(lonr) / RO - F + sqrt(G)) + dx;
	}
	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_api coordinate function, ", "X_api > MAX_FLOAT: ", X);
	
	return X;

}

template <typename T>
T Projections::Y_api(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_api coordinate function, ", "Y_api > MAX_FLOAT: ", Y);

	return Y;
}


//Apianus Elliptical
template <typename T>
T Projections::X_apiel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = 2.0 * lat / (RO * M_PI);

	//Throw exception
	if (A > 1)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate acos(A) in X_apiel coordinate function, ", "A > 1: ", A);


	const T X = R * lonr / RO * cos(asin(A)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_apiel coordinate function, ", "X_apiel > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_apiel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_apiel coordinate function, ", "Y_apiel > MAX_FLOAT: ", Y);

	return Y;

}


//Armadillo
template <typename T>
T Projections::X_armad(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T lat_int = -atan(cos(0.5 * lonr / RO) / tan(20 / RO)) * RO;
	const T X =  (lat < lat_int ? R * (1.0 + cos(lat_int / RO)) * sin(lonr / 2.0 / RO) + dx : R * (1.0 + cos(lat / RO)) * sin(lonr / 2.0 / RO) + dx);

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_armad coordinate function, ", "X_armad > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_armad(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T lat_int = -atan(cos(0.5 * lonr / RO) / tan(20 / RO)) * RO;
	const T Y = (lat < lat_int ? R * ((1.0 + sin(20 / RO) - cos(20 / RO)) / 2.0 + sin(lat_int / RO) * cos(20 / RO) - (1.0 + cos(lat_int / RO)) * sin(20 / RO) * cos(lonr / 2.0 / RO)) + dy :
		                  R * ((1.0 + sin(20 / RO) - cos(20 / RO)) / 2.0 + sin(lat / RO) * cos(20 / RO) - (1.0 + cos(lat / RO)) * sin(20 / RO) * cos(lonr / 2.0 / RO)) + dy);

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_armad coordinate function, ", "Y_armad > MAX_FLOAT: ", Y);

	return Y;

}


//August Epicycloidal
template <typename T>
T Projections::X_august(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T C1 = sqrt(1.0 - tan(0.5 * lat / RO) * tan(0.5 * lat / RO));
	const T C = 1.0 + C1 * cos(0.5 * lonr / RO);

	//Throw exception
	if (fabs(C) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_august coordinate function, ", "1.0 / C, C = 0:", C);

	const T X1 = sin(0.5 * lonr / RO) * C1 / C;
	const T Y1 = tan(0.5 * lat / RO) / C;

	const T X = 4.0 * R * X1 * (3.0 + X1 * X1 - 3.0 * Y1 * Y1) / 3.0 + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_august coordinate function, ", "X_armad > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_august(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T C1 = sqrt(1.0 - tan(0.5 * lat / RO) * tan(0.5 * lat / RO));
	const T C = 1.0 + C1 * cos(0.5 * lonr / RO);

	//Throw exception
	if (fabs(C) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_august coordinate function, ", "1.0 / C, C = 0:", C);

	const T X1 = sin(0.5 * lonr / RO) * C1 / C;
	const T Y1 = tan(0.5 * lat / RO) / C;

	const T Y = 4.0 * R * Y1 * (3.0 + 3.0 * X1 * X1 - Y1 * Y1) / 3.0 + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_august coordinate function, ", "Y_armad > MAX_FLOAT: ", Y);

	return Y;
}


//Bacon Globular
template <typename T>
T Projections::X_bacon(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs(lonr) < MAX_ANGULAR_DIFF)
	{
		X = dx;
	}

	else
	{
		const T F = ((M_PI / 2) *(M_PI / 2) * RO / abs(lonr) + abs(lonr) / RO) / 2;
		const T Y = R * M_PI / 2.0 * sin(lat / RO) + dy;
		const T G = F * F - Y * Y / (R * R);

		//Throw exception
		if (G < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(G) in X_bacon coordinate function, ", "G < 0: ", G);

		X = R * sign(lonr) * (abs(lonr) / RO - F + sqrt(G)) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_bacon coordinate function, ", "X_bacon > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_bacon(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * (M_PI / 2) * sin(lat / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_bacon coordinate function, ", "Y_bacon > MAX_FLOAT: ", Y);

	return Y;
}


//Boggs Eumorphic
template <typename T>
T Projections::X_boggs(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF)
	{
		X = dx;
	}

	else
	{
		const T theta0 = lat;
		const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);
		
		const T A = 1.0 / cos(lat / RO) + 1.11072 / cos(theta / RO);

		//Throw exception
		if (fabs(A) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_boggs coordinate function, ", "1.0 / A, A = 0:", A);

		X = 2.00276 * R * lonr / RO / A + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_boggs coordinate function, ", "X_boggs > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_boggs(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	if (fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF)
	{
		Y = 0.49931 * R * (lat / RO + sqrt(2)) + dy;
	}

	else
	{
		const T theta0 = lat;
		const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);

		Y = 0.49931 * R * (lat / RO + sqrt(2) * sin(theta / RO)) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_boggs coordinate function, ", "Y_boggs > MAX_FLOAT: ", Y);

	return Y;
}


//Bonne
template <typename T>
T Projections::X_bonne(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(lat1) == MAX_LAT)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat1) in X_bonne coordinate function, ", "lat1 = +-90: ", lat1);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_bonne coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = 1.0 / A + (lat1 - lat) / RO;
	const T X = R * B * sin(lonr / RO * cos(lat / RO) / B) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_bonne coordinate function, ", "X_bonne > MAX_FLOAT: ", X);
	
	return X;
}

template <typename T>
T Projections::Y_bonne(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(lat1) == MAX_LAT)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat1) in Y_bonne coordinate function, ", "lat1 = +-90: ", lat1);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_bonne coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = 1.0 / A + (lat1 - lat) / RO;

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_bonne coordinate function, ", "1.0 / B, B = 0:", B);

	const T Y = R * 1.0 / A - R * B * cos(lonr / RO * cos(lat / RO) / B) + dy;
	
	return Y;
}


//Breusign
template <typename T>
T Projections::X_breus(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 2.0 / RO;
	const T B = cos(A);

	//Throw exception
	if (B < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(B) in X_breus coordinate function, ", "B < 0: ", B);

	const T X = 2.0 * R * sin(A) * sqrt(B) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_breus coordinate function, ", "X_breus > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_breus(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 2.0 / RO;
	const T B = cos(A);

	//Throw exception
	if (B < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(B) in Y_breus coordinate function, ", "B < 0: ", B);

	const T Y = -2.0 * R * sin(A) * sqrt(B) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_breus coordinate function, ", "Y_breus > MAX_FLOAT: ", Y);

	return Y;
}


//Central Cylindrical
template <typename T>
T Projections::X_cc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr / RO * cos (lat1 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_cc coordinate function, ", "X_cc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_cc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	//Throw exception
	if (fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_cc coordinate function, ", "lat = +-90: ", lat);

	const T Y = R * tan(lat / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_cc coordinate function, ", "Y_cc > MAX_FLOAT: ", Y);

	return Y;
}


//Cylindrical Equal Area
template <typename T>
T Projections::X_cea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr / RO * cos(lat1 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_cea coordinate function, ", "X_cea > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_cea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * sin(lat / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_cea coordinate function, ", "Y_cea > MAX_FLOAT: ", Y);

	return Y;
}


//Clark Perspective
template <typename T>
T Projections::X_clar(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.368 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_clar coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = R * 2.368 * R * cos(lat / RO) / A * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_clar coordinate function, ", "X_clar > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_clar(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.368 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_clar coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = -R * 2.368 * R * cos(lat / RO) / A * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_clar coordinate function, ", "Y_clar > MAX_FLOAT: ", Y);

	return Y;
}


//Collignon
template <typename T>
T Projections::X_collg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * 2.0 / sqrt(M_PI) * lonr / RO * sqrt(1.0 - sin(lat / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_collg coordinate function, ", "X_collg > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_collg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * sqrt(M_PI) * (1.0 - sqrt(1.0 - sin(lat / RO))) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_collg coordinate function, ", "Y_collg > MAX_FLOAT: ", Y);

	return Y;
}


//Craster Parabolic
template <typename T>
T Projections::X_crast(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * sqrt(3.0 / M_PI) * lonr / RO * (2.0 * cos(2.0 * lat / 3.0 / RO) - 1) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_crast coordinate function, ", "X_crast > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_crast(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * sqrt(3.0 * M_PI) * sin(lat / 3.0 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_crast coordinate function, ", "Y_crast > MAX_FLOAT: ", Y);

	return Y;
}


//Denoyer Semi-Elliptical
template <typename T>
T Projections::X_denoy(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X =  R * lonr / RO * cos((0.95 - 1.0 / 12 * abs(lonr) / RO + 1.0 / 600 * pow((abs(lonr) / RO), 3)) * lat / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_denoy coordinate function, ", "X_denoy > MAX_FLOAT: ", X);

	return X;
}


template <typename T>
T Projections::Y_denoy(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_denoy coordinate function, ", "Y_denoy > MAX_FLOAT: ", Y);

	return Y;
}


//Eckert I
template <typename T>
T Projections::X_eck1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 2.0 * sqrt(2.0 / (3.0 * M_PI)) * R * lonr / RO * (1.0 - fabs(lat / (M_PI * RO))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck1 coordinate function, ", "X_eck1 > MAX_FLOAT: ", X);

	return X;

}

template <typename T>
T Projections::Y_eck1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * sqrt(2.0 / (3.0 * M_PI)) * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck1 coordinate function, ", "Y_eck1 > MAX_FLOAT: ", Y);

	return Y;
}


//Eckert II
template <typename T>
T Projections::X_eck2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 2.0 * R / sqrt(6.0 * M_PI) * lonr / RO * sqrt(4.0 - 3.0 * sin(fabs(lat / RO))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck2 coordinate function, ", "X_eck2 > MAX_FLOAT: ", X);

	return X;

}

template <typename T>
T Projections::Y_eck2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = sqrt(2.0 * M_PI / 3) * R* (2.0 - sqrt(4.0 - 3.0 * sin(fabs(lat / RO)))) * sign(lat) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck2 coordinate function, ", "Y_eck2 > MAX_FLOAT: ", Y);

	return Y;
}


//Eckert III
template <typename T>
T Projections::X_eck3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = 1.0 - 4.0 * pow((lat / (M_PI * RO)), 2);

	//Throw exception
	if (A < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(A) in X_eck3 coordinate function, ", "A < 0: ", A);

	const T X = 2.0 / sqrt(M_PI * (4.0 + M_PI)) * R * (1.0 + sqrt(A)) * lonr / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck3 coordinate function, ", "X_eck3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eck3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 4.0 / sqrt(M_PI * (4.0 + M_PI)) * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck3 coordinate function, ", "Y_eck3 > MAX_FLOAT: ", Y);

	return Y;
}


//Eckert IV
template <typename T>
T Projections::X_eck4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T latr = lat / RO;
	const T theta0 = (0.895168 * latr + 0.0218849 * latr * latr * latr + 0.00806809 * latr * latr * latr * latr * latr) * RO;
	const T theta = NewtonRaphson(FTheta_eck4<T>, FThetaDer_eck4<T>, lat, theta0);

	const T X = 2.0 * R * lonr / RO * (1.0 + cos(theta / RO)) / sqrt(M_PI * (4.0 + M_PI)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck4 coordinate function, ", "X_eck4 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eck4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T latr = lat / RO;
	const T theta0 = (0.895168 * latr + 0.0218849 * latr * latr * latr + 0.00806809 * latr * latr * latr * latr * latr) * RO;
	const T theta = NewtonRaphson(FTheta_eck4<T>, FThetaDer_eck4<T>, lat, theta0);

	const T Y = 2.0 * R * sqrt(M_PI / (4.0 + M_PI)) * sin(theta / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck4 coordinate function, ", "Y_eck4 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_eck4(const T lat, const T theta)
{
	return theta / RO + sin(theta / RO) * (cos(theta / RO) + 2) - (2.0 + M_PI / 2) * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_eck4(const T lat, const T theta)
{
	return (1.0 + cos(theta / RO) * (cos( theta / RO) + 2) - sin(theta / RO) * sin(theta / RO))/RO;
}


//Eckert V
template <typename T>
T Projections::X_eck5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr / RO * (1.0 + cos(lat / RO)) / sqrt(2.0 + M_PI) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck5 coordinate function, ", "X_eck5 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eck5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * R * lat / RO / sqrt(2.0 + M_PI) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck5 coordinate function, ", "Y_eck5 > MAX_FLOAT: ", Y);

	return Y;
}


//Eckert VI
template <typename T>
T Projections::X_eck6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_eck6<T>, FThetaDer_eck6<T>, lat, theta0);

	const T X = R * (1.0 + cos(theta / RO)) * lonr / RO / sqrt(2.0 + M_PI) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eck6 coordinate function, ", "X_eck6 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eck6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_eck6<T>, FThetaDer_eck6<T>, lat, theta0);

	const T Y = 2.0 * R * theta / sqrt(2.0 + M_PI) / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eck6 coordinate function, ", "Y_eck6 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_eck6(const T lat, const T theta)
{
	return theta / RO + sin(theta / RO) - (1.0 + M_PI / 2) * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_eck6(const T lat, const T theta)
{
	return (1.0 + cos(theta / RO)) / RO;
}


//Eisenlohr
template <typename T>
T Projections::X_eisen(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T S1 = sin(0.5 * lonr / RO);
	const T C1 = cos(0.5 * lonr / RO);
	const T D = cos(0.5 * lat / RO) + C1 * sqrt(2.0 * cos(lat / RO));

	//Throw exception
	if (fabs(D) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eisen coordinate function, ", "1.0 / D, D = 0:", D);
	
	const T TT = sin(0.5 * lat / RO) / D;
	const T C = sqrt(2.0 / (1.0 + TT * TT));
	const T E = cos(0.5 * lat / RO) + sqrt(0.5 * cos(lat / RO)) * (C1 + S1);
	const T F = cos(0.5 * lat / RO) + sqrt(0.5 * cos(lat / RO)) * (C1 - S1);

	//Throw exception
	if (fabs(F) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eisen coordinate function, ", "1.0 / E, E = 0:", E);

	const T G = E / F;

	if (G < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(G) in X_eisen coordinate function, ", "G < 0: ", G);

	const T V = sqrt(G);

	//Throw exception
	if (fabs(V) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eisen coordinate function, ", "1.0 / V, V = 0:", V);

	const T X = R * (3.0 + sqrt(8))*(-2.0 * log(V) + C * (V - 1.0 / V)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eisen coordinate function, ", "X_eisen > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eisen(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T S1 = sin(0.5 * lonr / RO);
	const T C1 = cos(0.5 * lonr / RO);
	const T D = cos(0.5 * lat / RO) + C1 * sqrt(2.0 * cos(lat / RO));

	//Throw exception
	if (fabs(D) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eisen coordinate function, ", "1.0 / D, D = 0:", D);

	const T TT = sin(0.5 * lat / RO) / D;
	const T C = sqrt(2.0 / (1.0 + TT * TT));
	const T E = cos(0.5 * lat / RO) + sqrt(0.5 * cos(lat / RO)) * (C1 + S1);
	const T F = cos(0.5 * lat / RO) + sqrt(0.5 * cos(lat / RO)) * (C1 - S1);

	//Throw exception
	if (fabs(F) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eisen coordinate function, ", "1.0 / E, E = 0:", E);

	const T G = E / F;

	//Throw exception
	if (G < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(G) in Y_eisen coordinate function, ", "G < 0: ", G);

	const T V = sqrt(G);

	//Throw exception
	if (fabs(V) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eisen coordinate function, ", "1.0 / V, V = 0:", V);

	//Throw exception
	if ((fmod(fabs(TT), 90) < MIN_FLOAT) && (fabs(TT) > MIN_FLOAT))
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate atan(TT) in Y_eisen coordinate function, ", "TT = PI: ", TT);

	const T Y = R * (3.0 + sqrt(8)) * (-2.0 * atan(TT) + C * TT * (V + 1.0 / V)) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eisen coordinate function, ", "Y_eisen > MAX_FLOAT: ", Y);

	return Y;
}


//Equidistant Cylindrical
template <typename T>
T Projections::X_eqc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr * cos(lat1 / RO) / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eqc coordinate function, ", "X_eqc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eqc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eqc coordinate function, ", "Y_eqc > MAX_FLOAT: ", Y);

	return Y;
}


//Equidistant Conic (true parallel lat1)
template <typename T>
T Projections::X_eqdc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eqdc coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = (R / A + R * (lat1 - lat) / RO) * sin(sin(lat1 / RO) * lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eqdc coordinate function, ", "X_eqdc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eqdc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eqdc coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = R / A - (R / A + R * (lat1 - lat) / RO) * cos(sin(lat1 / RO) * lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eqdc coordinate function, ", "Y_eqdc > MAX_FLOAT: ", Y);

	return Y;
}


//Equidistant Conic (true parallels lat1, lat2)
template <typename T>
T Projections::X_eqdc2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	//Throw exception
	if (fabs(lat1 - 90) < MAX_ANGULAR_DIFF)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eqdc2 coordinate function, ", "1.0 / (90 - lat1), lat1 = 90.", lat1);

	const T X = (R * (90 - lat) / RO) * sin(cos(lat1 / RO) / (90 - lat1) * lonr) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eqdc2 coordinate function, ", "X_eqdc2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eqdc2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	//Throw exception
	if (fabs(lat1 - 90) < MAX_ANGULAR_DIFF)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eqdc2 coordinate function, ", "1.0 / (90 - lat1), lat1 = 90.", lat1);

	const T Y = -(R * (90 - lat) / RO) * cos(cos(lat1 / RO) / (90 - lat1) * lonr) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eqdc2 coordinate function, ", "Y_eqdc2 > MAX_FLOAT: ", Y);

	return Y;
}


//Equidistant Conic (true parallel lat1, pole = point)
template <typename T>
T Projections::X_eqdc3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat1 / RO) - cos(lat2 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_eqdc3 coordinate function, ", "1.0 / A, A = 0:", A);
	
	const T B = A / (lat2 - lat1) * RO;

	const T X = R * ((lat2 / RO * cos(lat1 / RO) - lat1 / RO * cos(lat2 / RO)) / A - lat / RO) * sin( B * lonr/RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_eqdc3 coordinate function, ", "X_eqdc3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_eqdc3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat1 / RO) - cos(lat2 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_eqdc3 coordinate function, ", "1.0 / A, A = 0:", A);
	
	const T B = A / (lat2 - lat1) * RO;

	const T Y = -R * ((lat2 / RO * cos(lat1 / RO) - lat1 / RO * cos(lat2 / RO)) / A - lat / RO) * cos(B * lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_eqdc3 coordinate function, ", "Y_eqdc3 > MAX_FLOAT: ", Y);

	return Y;
}


//Fahey
template <typename T>
T Projections::X_fahey(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat / 2.0 / RO);
	const T B = 1.0 - pow(A, 2);

	//Throw exception
	if (B < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(B) in X_fahey coordinate function, ", "B < 0: ", B);

	const T X = R * lonr / RO * cos(lat1 / RO) * sqrt(B) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_fahey coordinate function, ", "X_fahey > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_fahey(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * (1.0 + cos(lat1 / RO)) * tan(lat / 2.0 / RO) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_fahey coordinate function, ", "Y_fahey > MAX_FLOAT: ", Y);

	return Y;
}


//Foucaut Sine-Tangent
template <typename T>
T Projections::X_fouc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 2.0 / sqrt(M_PI) * R * lonr / RO * cos(lat / RO) * pow((cos(lat / 2.0 / RO)), 2) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_fouc coordinate function, ", "X_fouc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_fouc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = sqrt(M_PI) * R * tan(lat / 2.0 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_fouc coordinate function, ", "Y_fouc > MAX_FLOAT: ", Y);

	return Y;
}


//Foucaut Sinusoidal
template <typename T>
T Projections::X_fouc_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat / RO);
	const T B = sin(lat1 / RO);
	const T C = B + (1.0 - B) * A;

	//Throw exception
	if (fabs(C) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_fouc_s coordinate function, ", "1.0 / C, C = 0:", A);

	const T X = R * lonr / RO * A / C + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_fouc_s coordinate function, ", "X_fouc_s > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_fouc_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T A = sin(lat1 / RO);
	const T Y = R * (A * lat / RO + (1.0 - A) * sin(lat / RO)) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_foucs coordinate function, ", "Y_foucs > MAX_FLOAT: ", Y);

	return Y;
}


//Fournier I Globular
template <typename T>
T Projections::X_fourn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	T X = 0;
	const T lonr = CartTransformation::redLon0(lon, lon0);

	if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if (fabs(fabs(lonr) - 90) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO * cos(lat / RO) + dx;
	}
	
	else
	{
		//Call Y Fournier
		const T Y = Y_fourn(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
		const T N = 1.0 - (2.0 * Y / (M_PI  * R)) * (2.0 * Y / (M_PI  * R));

		//Throw exception
		if (N < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(N) in X_fourn coordinate function, ", "N < 0: ", N);


		X = R * lonr / RO  * sqrt(N) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_fourn coordinate function, ", "X_fourn > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_fourn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	T Y = 0;
	const T lonr = CartTransformation::redLon0(lon, lon0);

	if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF))
	{
		Y = R * lat / RO + dy;
	}

	else if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if (fabs(fabs(lonr) - 90) < MAX_ANGULAR_DIFF)
	{
		Y = R * M_PI / 2.0 * sin(lat / RO) + dy;
	}

	else
	{
		const T C = M_PI * M_PI / 4;
		const T P = M_PI * fabs(sin(lat / RO));
		const T Q = P - 2.0 * fabs(lat) / RO;

		// Throw exception
		if (fabs(Q) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_fourn coordinate function, ", "1.0 / Q, Q = 0.", Q);

		const T F = ((C - (lat / RO) * (lat / RO)) / Q);
		const T G = 2.0 * lonr / (M_PI * RO);
		const T L = (G * G - 1);

		// Throw exception
		if (fabs(L) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_fourn coordinate function, ", "1.0 / L, L = 0.", L);

		const T M = F * F - L * (C - F * P - (lonr / RO) * (lonr / RO));

		//Throw exception
		if (M < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(M) in Y_fourn coordinate function, ", "M < 0: ", M);

		// Throw exception
		if (fabs(lat) == MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_fourn coordinate function, ", "1.0 / lat, lat = 0.", lat);

		Y = R / L  * sign(lat) * (sqrt(M) - F) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_fourn coordinate function, ", "Y_fourn > MAX_FLOAT: ", Y);

	return Y;
}


//Fournier II Elliptical
template <typename T>
T Projections::X_fourn2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R / sqrt(M_PI) * lonr / RO * cos(lat / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_fourn2 coordinate function, ", "X_fourn2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_fourn2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * sqrt(M_PI) / 2.0 * sin(lat / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_fourn2 coordinate function, ", "Y_fourn2 > MAX_FLOAT: ", Y);

	return Y;
}


//Gall Stereographic
template <typename T>
T Projections::X_gall(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * cos(lat1 / RO) * lonr / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_gall coordinate function, ", "X_gall > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_gall(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	//Throw exception
	if (fabs(lat/2) == MAX_LAT)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat/2) in X_fourn coordinate function, ", "lat/2.0 = +-90: ", lat);

	const T Y = R * (1.0 + cos(lat1 / RO)) * tan(lat / 2.0 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_gall coordinate function, ", "Y_gall > MAX_FLOAT: ", Y);

	return Y;
}


//Ginsburg VIII (TsNIIGAiK)
template <typename T>
T Projections::X_gins8(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr / RO * (1.0 - 0.162388 * (lat / RO) * (lat / RO)) * (0.87 - 0.000952426 * pow((lonr / RO), 4)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_gins8 coordinate function, ", "X_gins8 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_gins8(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO * (1.0 + 1.0 / 12 * pow((abs(lat) / RO), 3)) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_gins8 coordinate function, ", "Y_gins8 > MAX_FLOAT: ", Y);

	return Y;
}


//Gnomonic
template <typename T>
T Projections::X_gnom(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	//Throw exception
	//if (lat <= 0)
	//	throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in X_gnom coordinate function, ", "lat <= 0: ", lat);

	const T X = R * tan((90 - lat) / RO) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_gnom coordinate function, ", "X_gnom > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_gnom(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	//Throw exception
	//if (lat <= 0)
	//	throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_gnom coordinate function, ", "lat <= 0: ", lat);

	const T Y = -R * tan((90 - lat) / RO) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_gnom coordinate function, ", "Y_gnom > MAX_FLOAT: ", Y);

	return Y;
}


//Goode Homolosine
template <typename T>
T Projections::X_goode(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lat_int = 40 + 44.0 / 60 + 11.8 / 3600;
	T X = 0;

	if (fabs(lat) <= lat_int)
	{
		X = X_sinu(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}
		
	else
	{
		X = X_moll(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_goode coordinate function, ", "X_goode > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_goode(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lat_int = 40 + 44.0 / 60;
	T Y = 0;

	if (fabs(lat) <= lat_int)
	{
		Y = Y_sinu(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	else
	{
		const T theta0 = lat;
		const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);

		Y = R * (sqrt(2) * sin(theta / RO) - 0.0528 * sign(lat)) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_goode coordinate function, ", "Y_goode > MAX_FLOAT: ", Y);

	return Y;
}


//Hammer
template <typename T>
T Projections::X_hammer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 2.0 * R * cos(lat / RO) * sin(lonr / 2.0 / RO) * sqrt(2) / (sqrt(1.0 + cos(lat / RO) * cos(lonr / 2.0 / RO))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_hammer coordinate function, ", "X_hammer > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_hammer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = 2.0 * R * sin(lat / RO) / (sqrt(2) * sqrt(1.0 + cos(lat / RO) * cos(lonr / 2.0 / RO))) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_hammer coordinate function, ", "Y_hammer > MAX_FLOAT: ", Y);

	return Y;
}


//Hatano Asymmetrical Equal Area
template <typename T>
T Projections::X_hataea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = 0.5 * lat;
	const T theta = NewtonRaphson(FTheta_hataea<T>, FThetaDer_hataea<T>, lat, theta0);
	
	const T X = 0.85 * R * lonr / RO * cos(theta / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_hataea coordinate function, ", "X_hataea > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_hataea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = 0.5 * lat;

	const T theta = NewtonRaphson(FTheta_hataea<T>, FThetaDer_hataea<T>, lat, theta0);

	const T Y = ( lat < 0 ? 1.93052 * R * sin(theta / RO) + dy : 1.56548 * R * sin(theta / RO) + dy );

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_hataea coordinate function, ", "Y_hataea > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_hataea(const T lat, const T theta)
{
	return ( lat < 0 ? 2.0 * theta / RO + sin(2.0 * theta / RO) - 2.43763 * sin(lat / RO) : 2.0 * theta / RO + sin(2.0 * theta / RO) - 2.67595 * sin(lat / RO));
}

template <typename T>
T Projections::FThetaDer_hataea(const T lat, const T theta)
{
	return (2.0 + 2.0 * cos(2.0 * theta / RO)) / RO;
}


//La Hire
template <typename T>
T Projections::X_hire(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * (2.0 + sqrt(2) / 2) * R * cos(lat / RO) / (R * (1.0 + sqrt(2) / 2) + R * sin(lat / RO)) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_hire coordinate function, ", "X_hire > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_hire(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = -R * (2.0 + sqrt(2) / 2) * R * cos(lat / RO) / (R * (1.0 + sqrt(2) / 2) + R * sin(lat / RO)) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_hire coordinate function, ", "Y_hire > MAX_FLOAT: ", Y);

	return Y;
}


//Interrupted Goode Homolosine
template <typename T>
T Projections::X_igh(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	int zone = 0;

	const double lon_10 = 10;
	const double lon_20 = 20;
	const double lon_30 = 30;
	const double lon_40 = 40;
	const double lon_50 = 50;
	const double lon_60 = 60;
	const double lon_80 = 80;
	const double lon_90 = 90;
	const double lon_100 = 100;
	const double lon_140 = 140;
	const double lon_160 = 160;
	const double lon_180 = 180;

	const T lat_int = 40 + 44.0 / 60 + 11.8 / 3600;

	const T lonr = CartTransformation::redLon0(lon, lon0);

	/*
	Zones:

	-180            -40                       180
	+--------------+-------------------------+    Zones 1,2,9,10,11 & 12:
	|1.0             |2.0                        |      Mollweide projection
	|              |                         |
	+--------------+-------------------------+    Zones 3,4,5,6,7 & 8:
	|3.0             |4.0                        |      Sinusoidal projection
	|              |                         |
	0 +-------+------+-+-----------+-----------+
	|5      |6       |7          |8          |
	|       |        |           |           |
	+-------+--------+-----------+-----------+
	|9      |10      |11         |12         |
	|       |        |           |           |
	+-------+--------+-----------+-----------+
	-180    -100      -20         80          180
	*/

	//Detect zone 1|2
	if (lat >= lat_int)
	{
		zone = (lonr <= -lon_40 ? 1.0 : 2);
	}
	//Detect zone 3|4
	else if (lat >= 0)
	{
		zone = (lonr <= -lon_40 ? 3.0 : 4);
	}
	//Detect zone 5|6|7|8
	else if (lat >= -lat_int)
	{
		if (lonr <= -lon_100) zone = 5;		// 5
		else if (lonr <= -lon_20) zone = 6;	// 6
		else if (lonr <= lon_80) zone = 7;	// 7
		else zone = 8;				// 8
	}

	//Detect zone 9|10|11|12
	else
	{
		if (lonr <= -lon_100) zone = 9;		// 9
		else if (lonr <= -lon_20) zone = 10;	// 10
		else if (lonr <= lon_80) zone = 11;	// 11
		else zone = 12;				// 12
	}

	T X = 0;

	//Sinusoidal projection zones: shift X coordinates
	if (zone == 3)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, -lon_100, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_100, 0.0, 0.0, 0.0, c);
	}
	
	else if (zone == 4)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, lon_30, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, lon_30, 0.0, 0.0, 0.0, c);
	}
	
	else if (zone == 5)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, -lon_160, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_160, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 6)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, -lon_60, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_60, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 7)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, lon_20, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, lon_20, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 8)
	{
		X = X_sinu(R, lat1, lat2, lat, lonr, lon_140, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0 ,lon_140, 0.0, 0.0, 0.0, c);
	}
	
	
	//Mollweide projection zones: shift Y coordinates
	if (zone == 1)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, -lon_100, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_100, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 2)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, lon_30, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, lon_30, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 9)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, -lon_160, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_160, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 10)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, -lon_60, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, -lon_60, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 11)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, lon_20, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, lon_20, 0.0, 0.0, 0.0, c);
	}

	else if (zone == 12)
	{
		X = X_moll(R, lat1, lat2, lat, lonr, lon_140, dx, dy, c) + X_sinu(R, lat1, lat2, 0.0, lon_140, 0.0, 0.0, 0.0, c);
	}
	
	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_igh coordinate function, ", "X_igh > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_igh(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	int zone = 0;

	const double lon_10 = 10;
	const double lon_20 = 20;
	const double lon_30 = 30;
	const double lon_40 = 40;
	const double lon_50 = 50;
	const double lon_60 = 60;
	const double lon_80 = 80;
	const double lon_90 = 90;
	const double lon_100 = 100;
	const double lon_140 = 140;
	const double lon_160 = 160;
	const double lon_180 = 180;

	const T lat_int = 40 + 44.0 / 60 + 11.8 / 3600;

	const T lonr = CartTransformation::redLon0(lon, lon0);

	/*
	Zones:

	-180            -40                       180
	+--------------+-------------------------+    Zones 1,2,9,10,11 & 12:
	|1.0             |2.0                        |      Mollweide projection
	|              |                         |
	+--------------+-------------------------+    Zones 3,4,5,6,7 & 8:
	|3.0             |4.0                        |      Sinusoidal projection
	|              |                         |
	0 +-------+------+-+-----------+-----------+
	|5      |6       |7          |8          |
	|       |        |           |           |
	+-------+--------+-----------+-----------+
	|9      |10      |11         |12         |
	|       |        |           |           |
	+-------+--------+-----------+-----------+
	-180    -100      -20         80          180
	*/

	//Detect zone 1|2
	if (lat >= lat_int)
	{
		zone = (lonr <= -lon_40 ? 1.0 : 2);
	}
	//Detect zone 3|4
	else if (lat >= 0)
	{
		zone = (lonr <= -lon_40 ? 3.0 : 4);
	}
	//Detect zone 5|6|7|8
	else if (lat >= -lat_int)
	{
		if (lonr <= -lon_100) zone = 5;		// 5
		else if (lonr <= -lon_20) zone = 6;	// 6
		else if (lonr <= lon_80) zone = 7;	// 7
		else zone = 8;				// 8
	}

	//Detect zone 9|10|11|12
	else
	{
		if (lonr <= -lon_100) zone = 9;		// 9
		else if (lonr <= -lon_20) zone = 10;	// 10
		else if (lonr <= lon_80) zone = 11;	// 11
		else zone = 12;				// 12
	}

	T Y = 0;

	//Sinusoidal projection zones
	if ((zone >= 3) || zone <= 7)
	{
		Y = Y_sinu(R, lat1, lat2, lat, lonr, 0.0, dx, dy, c);
	}

	//Mollweide projection zones
	else
	{
		const T theta0 = lat;
		const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);

		Y = R * (sqrt(2) * sin(theta / RO) - 0.0528 * sign(lat)) + dy;
	}

	return Y;
}


//James Perspective
template <typename T>
T Projections::X_jam(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.5 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_jam coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = R * 2.5 * R * cos(lat / RO) / A * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_jam coordinate function, ", "X_jam > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_jam(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.5 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_jam coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = -R * 2.5 * R * cos(lat / RO) / A * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_jam coordinate function, ", "Y_jam > MAX_FLOAT: ", Y);

	return Y;
}


//Kavrayskiy V
template <typename T>
T Projections::X_kav5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 1.35439 / 1.50488 * R * lonr / RO * cos(lat / RO) / cos(lat / 1.35439 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_kav5 coordinate function, ", "X_kav5 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_kav5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 1.50488 * R * sin(lat / 1.35439 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_kav5 coordinate function, ", "Y_kav5 > MAX_FLOAT: ", Y);

	return Y;
}


//Kavrayskiy VII
template <typename T>
T Projections::X_kav7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T A = 1.0 - 3.0 * pow((lat / (M_PI * RO)), 2);

	//Throw exception
	if (fabs(A) == MIN_FLOAT)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(A) in X_kav7 coordinate function, ", "A = 0: ", A);

	const T X = sqrt(3) / 2.0 * R * sqrt(A) * lonr / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_kav7 coordinate function, ", "X_kav7 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_kav7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_kav7 coordinate function, ", "Y_kav7 > MAX_FLOAT: ", Y);

	return Y;
}


//Lambert Azimuthal Equal Area
template <typename T>
T Projections::X_laea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 2.0 * R * sin((90 - lat) / 2.0 / RO) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_laea coordinate function, ", "X_laea > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_laea(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = -2.0 * R * sin((90 - lat) / 2.0 / RO) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_laea coordinate function, ", "Y_laea > MAX_FLOAT: ", Y);

	return Y;
}


//Lagrange Conformal
template <typename T>
T Projections::X_lagrng(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs(fabs(lat) - MAX_LAT) < MAX_ANGULAR_DIFF) 
	{
		X = dx;
	}

	else
	{
		const T W = 2.0;
		const T A1 = pow((1.0 + sin(lat1 / RO)) / (1.0 - sin(lat1 / RO)), 1.0 / (2.0 * W));
		const T A = pow((1.0 + sin(lat / RO)) / (1.0 - sin(lat / RO)), 1.0 / (2.0 * W));

		//Throw exception
		if (fabs(A1) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_lagrng coordinate function, ", "1.0 / A1, A1 = 0:", A1);

		const T V = A / A1;
		const T C = 0.5 * (V + 1.0 / V) + cos(lonr / W / RO);

		//Throw exception
		if (fabs(C) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_lagrng coordinate function, ", "1.0 / C, C = 0:", C);

		X = 2.0 * R * sin(lonr / W / RO) / C + dx;
	}
		
	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_lagrng coordinate function, ", "X_lagrng > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_lagrng(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	if (fabs(fabs(lat) - MAX_LAT) < MAX_ANGULAR_DIFF)
	{
		Y = 2.0 * R * sign(lat) + dy;
	}

	else
	{
		const T W = 2.0;
		const T A1 = pow((1.0 + sin(lat1 / RO)) / (1.0 - sin(lat1 / RO)), 1.0 / (2.0 * W));
		const T A = pow((1.0 + sin(lat / RO)) / (1.0 - sin(lat / RO)), 1.0 / (2.0 * W));
		
		//Throw exception
		if (fabs(A1) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_lagrng coordinate function, ", "1.0 / A1, A1 = 0:", A1);

		const T V = A / A1;
		const T C = 0.5 * (V + 1.0 / V) + cos(lonr / W / RO);

		//Throw exception
		if (fabs(C) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_lagrng coordinate function, ", "1.0 / C, C = 0:", C);

		//Throw exception
		if (fabs(V) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_lagrng coordinate function, ", "1.0 / V, V = 0:", V);

		Y = R * (V - 1.0 / V) / C + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_lagrng coordinate function, ", "Y_lagrng > MAX_FLOAT: ", Y);

	return Y;
}


//Larrivee
template <typename T>
T Projections::X_larr(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat / RO);

	if (A < 0)
		throw MathInvalidArgumentException <T> ("MathInvalidArgumentException: can not evaluate sqrt(A) in X_larr coordinate function, ", "A < 0: ", A);

	const T X = 0.5 * R * lonr / RO * (1.0 + sqrt(A)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_larr coordinate function, ", "X_larr > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_larr(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = cos(lat / 2.0 / RO) * cos(lonr / 6 / RO);

	if (fabs(A) < MIN_FLOAT)
		throw MathZeroDevisionException<T>("MathDivisonByZeroException: can not evaluate Y_larr coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = R * lat / RO / A + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_larr coordinate function, ", "Y_larr > MAX_FLOAT: ", Y);

	return Y;
}


//Laskowski
template <typename T>
T Projections::X_lask(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * (0.975534 * lonr / RO - 0.119161 * (lonr / RO) * (lat / RO) * (lat / RO) - 0.0143059 * pow((lonr / RO), 3) * (lat / RO) * (lat / RO) -
		0.0547009 * (lonr / RO) * pow((lat / RO), 4)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_lask coordinate function, ", "X_lask > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_lask(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = R * (1.00384 * (lat / RO) + 0.0802894 * (lonr / RO) * (lonr / RO) * (lat / RO) + 0.0998909 * pow((lat / RO), 3) + 
		0.000199025 * pow((lonr / RO), 4) * (lat / RO) - 0.0285500 * (lonr / RO) * (lonr / RO) * pow((lat / RO), 3) - 0.0491032 * pow((lat / RO), 5)) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_lask coordinate function, ", "Y_lask > MAX_FLOAT: ", Y);

	return Y;
}


//Lambert Conformal Conic
template <typename T>
T Projections::X_lcc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = sin(lat1 / RO);
	const T X = R / tan(lat1 / RO) * pow((tan((lat1 / 2.0 + 45) / RO) / tan((lat / 2.0 + 45) / RO)), (A)) * sin(A * lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_lcc coordinate function, ", "X_lcc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_lcc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_lcc coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = sin(lat1 / RO);
	const T Y = R / A - R / A * pow((tan((lat1 / 2.0 + 45) / RO) / tan((lat / 2.0 + 45) / RO)), B) * cos(B * lonr / RO) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_lcc coordinate function, ", "Y_lcc > MAX_FLOAT: ", Y);

	return Y;
}


//Lambert Equal Area Conic (standard parallel lat1)
template <typename T>
T Projections::X_leac(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_leac coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = sin(lat1 / RO);

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_leac coordinate function, ", "1.0 / B, B = 0:", B);

	const T X = sqrt(R / A * R / A + 2.0 * R * R / B * (B - sin(lat / RO))) * sin(B * lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_leac coordinate function, ", "X_leac > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_leac(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = tan(lat1 / RO);

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_leac coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = sin(lat1 / RO);

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_leac coordinate function, ", "1.0 / B, B = 0:", B);

	const T Y = R / A - sqrt(R / A * R / A + 2.0 * R * R / B * (B - sin(lat / RO))) * cos(B * lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_leac coordinate function, ", "Y_leac > MAX_FLOAT: ", Y);

	return Y;
}


//Lambert Equal Area Conic (standard parallels lat1, lat2)
template <typename T>
T Projections::X_leac2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = sin(lat1 / RO);

	//Throw exception
	if (fabs(A + 1) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_leac2 coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = 2.0 / ((1.0 + A) / 2) * (1.0 - sin(lat / RO));

	//Throw exception
	if (B < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(B) in X_leac2 coordinate function, ", "A < 0: ", B);

	const T X = R * sqrt(B) * sin((1.0 + A) / 2.0 * lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_leac2 coordinate function, ", "X_leac2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_leac2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = sin(lat1 / RO);

	//Throw exception
	if (fabs(A + 1) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_leac2 coordinate function, ", "1.0 / A, A = 0:", A);

	const T B = (1.0 - A) / (1.0 + A);

	//Throw exception
	if (B < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(B) in Y_leac2 coordinate function, ", "B < 0: ", B);

	const T C = 2.0 / ((1.0 + A) / 2) * (1.0 - sin(lat / RO));

	//Throw exception
	if (C < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(C) in Y_leac2 coordinate function, ", "C < 0: ", C);

	const T Y = 2.0 * R * sqrt(B) - R * sqrt(C) * cos((1.0 + A) / 2.0 * lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_leac2 coordinate function, ", "Y_leac2 > MAX_FLOAT: ", Y);

	return Y;
}


//Loximuthal
template <typename T>
T Projections::X_loxim(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	const T dlat = lat - lat1;
	if (fabs(dlat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO * cos(lat1 / RO) + dx;
	}

	else
	{
		const T A = 0.5 * lat + 45;
		const T B = 0.5 * lat1 + 45;

		if ((fabs(A) < MAX_ANGULAR_DIFF) || (fabs(fabs(A) - 90) < MAX_ANGULAR_DIFF))
		{
			X = dx;
		}

		else
		{
			X = R * lonr / RO * dlat / RO / (log(tan(A / RO)) - log(tan(B / RO))) + dx;
		}
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_loxim coordinate function, ", "X_loxim > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_loxim(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = R * (lat - lat1) / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_loxim coordinate function, ", "Y_loxim > MAX_FLOAT: ", Y);

	return Y;
}


//McBryde-Thomas Sine I.
template <typename T>
T Projections::X_mbt_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 1.36509 / 1.48875 * R * lonr / RO * cos(lat / RO) / cos(lat / 1.36509 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_mbt_s coordinate function, ", "X_mbt_s > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_mbt_s(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 1.48875 * R * sin(lat / 1.36509 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_mbt_s coordinate function, ", "Y_mbt_s > MAX_FLOAT: ", Y);

	return Y;
}


//McBryde-Thomas Flat-Pole Sine III
template <typename T>
T Projections::X_mbt_s3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lat_int = 55 + 51.0 / 60;
	T X = 0;

	if (fabs(lat) <= lat_int)
	{
		X = X_sinu(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	else
	{
		X = X_mbtfps(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_mbt_s3 coordinate function, ", "X_mbt_s3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_mbt_s3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lat_int = 55 + 51.0 / 60;
	T Y = 0;

	if (fabs(lat) <= lat_int)
	{
		Y = Y_sinu(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	else
	{
		const T theta0 = lat;
		const T theta = NewtonRaphson(FTheta_mbtfps<T>, FThetaDer_mbtfps<T>, lat, theta0);

		Y = R * (sqrt(6 / (4.0 + M_PI)) * theta / RO - 0.069065 * sign(lat)) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_mbt_s3 coordinate function, ", "Y_mbt_s3 > MAX_FLOAT: ", Y);

	return Y;
}


//McBryde-Thomas Flat-Pole Quartic IV
template <typename T>
T Projections::X_mbtfpq(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_mbtfpq<T>, FThetaDer_mbtfpq<T>, lat, theta0);
	
	const T X = R * lonr / RO * (1.0 + 2.0 * cos(theta / RO) / cos(0.5 * theta / RO)) / sqrt(3.0 * sqrt(2) + 6) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_mbtfpq coordinate function, ", "X_mbtfpq > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_mbtfpq(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_mbtfpq<T>, FThetaDer_mbtfpq<T>, lat, theta0);
	
	const T Y = 2.0 * sqrt(3) * R * sin(theta / 2.0 / RO) / sqrt(2.0 + sqrt(2)) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_mbtfpq coordinate function, ", "Y_mbtfpq > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_mbtfpq(const T lat, const T theta)
{
	return sin(theta / 2.0 / RO) + sin(theta / RO) - (1.0 + sqrt(2) / 2) * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_mbtfpq(const T lat, const T theta)
{
	return ( 0.5 * cos(theta / 2.0 / RO) + cos(theta / RO)) / RO;
}


//McBryde-Thomas Flat-Pole Sine II
template <typename T>
T Projections::X_mbtfps(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_mbtfps<T>, FThetaDer_mbtfps<T>, lat, theta0);

	const T X = R * sqrt(6 / (4.0 + M_PI)) * (0.5 + cos(theta / RO)) * lonr / RO / 1.5 + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_mbtfps coordinate function, ", "X_mbtfps > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_mbtfps(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_mbtfps<T>, FThetaDer_mbtfps<T>, lat, theta0);

	const T Y = R * sqrt(6 / (4.0 + M_PI)) * theta / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_mbtfs coordinate function, ", "Y_mbtfps > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_mbtfps(const T lat, const T theta)
{
	return theta / 2.0 / RO + sin(theta / RO) - (1.0 + M_PI / 4) * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_mbtfps(const T lat, const T theta)
{
	return ( 0.5 + cos(theta / RO)) / RO;
}


//Mercator
template <typename T>
T Projections::X_merc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr * cos(lat1 / RO) / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_merc coordinate function, ", "X_merc > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_merc(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T A = tan((lat / 2.0 + 45) / RO);

	//Throw exception
	if (fabs(fabs(lat) - MAX_LAT) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_merc coordinate function, ", "lat = +-90: ", lat);

	//Throw exception
	if (A <= 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate log(A) in Y_merc coordinate function, ", "A <= 0: ", A);

	const T Y = R * log(A) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_merc coordinate function, ", "Y_merc > MAX_FLOAT: ", Y);

	return Y;
}


//Miller Cylindrical
template <typename T>
T Projections::X_mill(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * lonr * cos(lat1 / RO) / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_mill coordinate function, ", "X_mill > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_mill(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T A = tan((0.4 * lat + 45.0) / RO);

	//Throw exception
	if (A <= 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate log(A) in Y_mill coordinate function, ", "A <= 0: ", A);

	const T Y = R * log(A) / 0.8 + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_mill coordinate function, ", "Y_mill > MAX_FLOAT: ", Y);

	return Y;
}


//Mollweide
template <typename T>
T Projections::X_moll(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);

	const T X = 2.0 * R * lonr / RO * sqrt(2) * cos(theta / RO) / M_PI + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_moll coordinate function, ", "X_moll > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_moll(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_moll<T>, FThetaDer_moll<T>, lat, theta0);

	const T Y = R * sqrt(2.0) * sin(theta / RO) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_moll coordinate function, ", "Y_moll > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_moll(const T lat, const T theta)
{
	return 2.0 * theta / RO + sin(2.0 * theta / RO) - M_PI * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_moll(const T lat, const T theta)
{
	return ( 2.0 + 2.0 * cos(2.0 * theta / RO)) / RO;
}


//Nell
template <typename T>
T Projections::X_nell(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_nell<T>, FThetaDer_nell<T>, lat, theta0);

	const T X = R * 0.5 * lonr / RO * (1.0 + cos(theta / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_nell coordinate function, ", "X_nell > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_nell(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_nell<T>, FThetaDer_nell<T>, lat, theta0);

	const T Y = R * theta / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_nell coordinate function, ", "Y_nell > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_nell(const T lat, const T theta)
{
	return  theta / RO + sin(theta / RO) - 2.0 * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_nell(const T lat, const T theta)
{
	return (1.0 + cos(theta / RO)) / RO;
}


//Nell-Hammer
template <typename T>
T Projections::X_nell_h(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T X = 0.5 * R * lonr / RO * (1.0 + cos(lat / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_nell_h coordinate function, ", "X_nell_h > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_nell_h(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * R * (lat / RO - tan(lat / 2.0 / RO)) + dy;

	//Throw exception
	if (fabs(fabs(lat) - 2.0 * MAX_LAT) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in X_nell_h coordinate function, ", "lat = +-180: ", lat);

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_nell_h coordinate function, ", "Y_nell_h > MAX_FLOAT: ", Y);

	return Y;
}


//Nicolosi Globular
template <typename T>
T Projections::X_nicol(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if ((fabs(lonr) < MAX_ANGULAR_DIFF) || ( fabs(fabs(lat) - 90) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if ( fabs(fabs(lonr) - 90) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO * cos(lat / RO) + dx;
	}

	else
	{
		const T B = M_PI / (2.0 * lonr / RO) - 2.0 * lonr / RO / M_PI;
		const T C = 2.0 * lat / RO / M_PI;
		const T D = sin(lat / RO) - C;

		//Throw exception
		if (fabs(D) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_nicol coordinate function, ", "1.0 / D, D = 0:", D);

		const T E = (1.0 - C * C) / D;

		//Throw exception
		if (fabs(E) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_nicol coordinate function, ", "1.0 / E, E = 0:", E);

		const T F = 1.0 + B * B / (E * E);
		const T G = cos(lat / RO);

		//Throw exception
		if (fabs(F) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_nicol coordinate function, ", "1.0 / F, F = 0:", F);

		const T M = (B * sin(lat / RO) / E - B / 2) / F;
		const T N = M * M + G * G / F;

		//Throw exception
		if (N < 0)
		{
			//Avoiding numerical imprecision
			//if (fabs(N) < EPS)
			//	N = 0;
			//else
				throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(N) in X_nicol coordinate function, ", "N < 0: ", N);
		}

		X = R * M_PI / 2.0 * (M + sign(lonr) * sqrt(N)) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_nicol coordinate function, ", "X_nicol > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_nicol(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0; 

	if ((fabs(lonr) < MAX_ANGULAR_DIFF) || ( fabs(fabs(lat) -90) < MAX_ANGULAR_DIFF))
	{
		Y = R * lat / RO + dy;
	}

	else if ( fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if ( fabs(fabs(lonr) - 90) < MAX_ANGULAR_DIFF)
	{
		Y = R * M_PI / 2.0 * sin(lat / RO) + dy;
	}

	else
	{
		const T B = M_PI / (2.0 * lonr / RO) - 2.0 * lonr / RO / M_PI;

		//Throw exception
		if (fabs(B) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_nicol coordinate function, ", "1.0 / B, B = 0:", B);

		const T C = 2.0 * lat / RO / M_PI;
		const T D = sin(lat / RO) - C;

		//Throw exception
		if (fabs(D) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_nicol coordinate function, ", "1.0 / D, D = 0:", D);

		const T E = (1.0 - C * C) / D;
		const T F = 1.0 + E * E / (B * B);

		//Throw exception
		if (fabs(F) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_nicol coordinate function, ", "1.0 / F, F = 0:", F);
		
		const T G = sin(lat / RO);
		const T N = (E * E * G / (B * B) + E / 2.0) / F;
		T P = N * N - (E * E * G * G / (B * B) + E * G - 1.0) / F;

		//Throw exception
		if (P < 0)
		{
			//Avoiding numerical imprecision
			if (fabs(P) < EPS)
				P = 0;
			else
				throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(P) in X_nicol coordinate function, ", "P < 0: ", P);
		}

		Y = R * M_PI / 2.0 * (N + sign(-lat) * sqrt(P)) + dy;
	
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_nicol coordinate function, ", "Y_nicol > MAX_FLOAT: ", Y);

	return Y;
}


//Ortelius
template <typename T>
T Projections::X_ortel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs(lonr) < 90)
	{
		//Call Apian projection
		X = X_api(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	}

	else
	{
		X = R * sign(lonr) * (sqrt(M_PI * M_PI / 4.0 - lat * lat / (RO * RO)) + abs(lonr / RO) - M_PI / 2) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_ortel coordinate function, ", "X_ortel > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_ortel(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	T Y = R * lat / RO + dy;
 
	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_ortho coordinate function, ", "Y_ortho > MAX_FLOAT: ", Y);

	return Y;
}


//Orthographic
template <typename T>
T Projections::X_ortho(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * cos(lat / RO) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_ortho coordinate function, ", "X_ortho > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_ortho(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = -R * cos(lat / RO) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_ortho coordinate function, ", "Y_ortho > MAX_FLOAT: ", Y);

	return Y;
}


//Parabolic
template <typename T>
T Projections::X_parab(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = sqrt(3.0 / M_PI) * R * lonr / RO * (2.0 * cos(2.0 * lat / 3.0 / RO) - 1) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_parab coordinate function, ", "X_parab > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_parab(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = sqrt(3.0 * M_PI) * R * sin(lat / 3.0 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_parab coordinate function, ", "Y_parab > MAX_FLOAT: ", Y);

	return Y;
}


//Perspective Azimuthal
template <typename T>
T Projections::X_pers(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = R * cos(lat / RO) * sqrt(5) / (sqrt(5) + 1.0 - sin(lat / RO)) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_pers coordinate function, ", "X_pers > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_pers(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T Y = -R * cos(lat / RO) * sqrt(5) / (sqrt(5) + 1.0 - sin(lat / RO)) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_pers coordinate function, ", "Y_pers > MAX_FLOAT: ", Y);

	return Y;
}


//Far-side Perspective Azimuthal
template <typename T>
T Projections::X_persf(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * ((1.0 + c) + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_persf coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = R * (2.0 + c) * R * cos(lat / RO) / A * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_persf coordinate function, ", "X_persf > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_persf(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * ((1.0 + c) + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_persf coordinate function, ", "1.0 / A, A = 0:", A);
	
	const T Y = -R * (2.0 + c) * R * cos(lat / RO) / A * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_persf coordinate function, ", "Y_persf > MAX_FLOAT: ", Y);

	return Y;
}


//Near-side Perspective Azimuthal
template <typename T>
T Projections::X_persn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T A = R * ((1.0 + c) + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_persn coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = R * c * R * cos(lat / RO) / A * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_persn coordinate function, ", "X_persn > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_persn(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * ((1.0 + c) + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_persn coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = -R * c * R * cos(lat / RO) / A * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_persn coordinate function, ", "Y_persn > MAX_FLOAT: ", Y);

	return Y;
}


//Hassler Polyconic, American
template <typename T>
T Projections::X_poly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if ( fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else
	{
		const T E = lonr / RO * sin(lat / RO);

		X = R * 1.0 / tan(lat / RO) * sin(E) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_poly coordinate function, ", "X_poly > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_poly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = -R * lat1 / RO + dy;
	}

	else
	{
		const T E = lonr / RO * sin(lat / RO);

		Y = R * (lat/RO - lat1/RO + 1.0 / tan(lat / RO) * (1.0 - cos(E))) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_poly coordinate function, ", "Y_poly > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P1
template <typename T>
T Projections::X_putp1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = 1.0 - 3.0 * pow((lat / (RO * M_PI)), 2);

	//Throw exception
	if (A < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(A) in X_putp1 coordinate function, ", "A < 0: ", A);


	const T X = 0.94745 * R * lonr / RO * sqrt(A) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp1 coordinate function, ", "X_putp1 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 0.94745 * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp1 coordinate function, ", "Y_putp1 > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P2
template <typename T>
T Projections::X_putp2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T latr = lat / RO;
	const T theta0 = ( 0.615709 * latr + 0.00909953 * latr * latr * latr + 0.0046292 * latr * latr * latr * latr * latr) * RO;
	
	const T theta = NewtonRaphson(FTheta_putp2<T>, FThetaDer_putp2<T>, lat, theta0);

	const T X = 1.97949 * R * lonr / RO * (cos(theta / RO) - 0.5) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp2 coordinate function, ", "X_putp2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T latr = lat / RO;
	const T theta0 = (0.615709 * latr + 0.00909953 * latr * latr * latr + 0.0046292 * latr * latr * latr * latr * latr) * RO;

	const T theta = NewtonRaphson(FTheta_putp2<T>, FThetaDer_putp2<T>, lat, theta0);

	const T Y = 1.71848 * R * sin(theta / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp2 coordinate function, ", "Y_putp2 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_putp2(const T lat, const T theta)
{
	return theta / RO + sin(theta/RO) * (cos(theta / RO) - 1) - (4.0 * M_PI - 3.0 * sqrt(3)) / 12.0 * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_putp2(const T lat, const T theta)
{
	return (1.0 + cos(theta / RO) * (cos(theta / RO) - 1) - sin(theta / RO) * sin(theta / RO)) / RO;
}


//Putnins P3
template <typename T>
T Projections::X_putp3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = sqrt(2.0 / M_PI) * R * lonr / RO * (1.0 - 4.0 * pow((lat / (RO * M_PI)), 2)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp3 coordinate function, ", "X_putp3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = sqrt(2.0 / M_PI) * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp3 coordinate function, ", "Y_putp3 > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P3P
template <typename T>
T Projections::X_putp3p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = sqrt(2.0 / M_PI) * R * lonr / RO * (1.0 - 2.0 * pow((lat / (RO * M_PI)), 2.0)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp3p coordinate function, ", "X_putp3p > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp3p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = sqrt(2.0 / M_PI) * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp3p coordinate function, ", "Y_putp3p > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P4P
template <typename T>
T Projections::X_putp4p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = sin(lat / RO);
	const T B = asin(5.0 * sqrt(2) / 8.0 * A);
	const T C = cos(B / 3.0);

	//Throw exception
	if (fabs(C) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_putp4p coordinate function, ", "1.0 / C, C = 0:", C);

	const T X = 2.0 * sqrt(0.6 / M_PI) * R * lonr / RO * cos(B) / C + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp4p coordinate function, ", "X_putp4p > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp4p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * R * sqrt(1.2 * M_PI) * sin(asin(5 * sqrt(2) / 8 * sin(lat / RO)) / 3) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp4p coordinate function, ", "Y_putp4p > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P5
template <typename T>
T Projections::X_putp5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T X = 1.01346 * R * lonr / RO * (2.0 - sqrt(1.0 + 12 * pow((lat / (RO * M_PI)), 2))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp5 coordinate function, ", "X_putp5 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 1.01346 * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp5 coordinate function, ", "Y_putp5 > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P5P
template <typename T>
T Projections::X_putp5p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T X = 1.01346 * R * lonr / RO * (1.5 - 0.5 * sqrt(1.0 + 12 * pow((lat / (RO * M_PI)), 2))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp5p coordinate function, ", "X_putp5p > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp5p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 1.01346 * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp5p coordinate function, ", "Y_putp5p > MAX_FLOAT: ", Y);

	return Y;
}


//Putnins P6
template <typename T>
T Projections::X_putp6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_putp6<T>, FThetaDer_putp6<T>, lat, theta0);

	const T X = 1.01346 * R * lonr / RO * (2.0 - sqrt( 1.0 + theta / RO * theta / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp6 coordinate function, ", "X_putp6 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_putp6<T>, FThetaDer_putp6<T>, lat, theta0);

	const T Y = 0.9191 * R * theta / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp6 coordinate function, ", "Y_putp6 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_putp6(const T lat, const T theta)
{
	const T p = 2.14714 * sin(lat / RO);
	const T r = sqrt(1.0 + theta / RO * theta / RO) ;

	return (4.0 - r) * theta / RO - log(theta / RO + r) - p;
}

template <typename T>
T Projections::FThetaDer_putp6(const T lat, const T theta)
{
	return (4.0 - 2.0 * sqrt(1.0 + theta / RO * theta / RO)) / RO;
}


//Putnins P6p
template <typename T>
T Projections::X_putp6p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_putp6<T>, FThetaDer_putp6<T>, lat, theta0);

	const T X = 0.44329 * R * lonr / RO * (3.0 - sqrt(1.0 + theta / RO * theta / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_putp6p coordinate function, ", "X_putp6p > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_putp6p(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	const T theta0 = lat;
	const T theta = NewtonRaphson(FTheta_putp6p<T>, FThetaDer_putp6p<T>, lat, theta0);

	const T Y = 0.80404 * R * theta / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_putp6p coordinate function, ", "Y_putp6p > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_putp6p(const T lat, const T theta)
{
	const T A = sqrt(1.0 + theta / RO * theta / RO);

	return (6.0 - A) * theta / RO - log(theta / RO + A) - 5.61125 * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_putp6p(const T lat, const T theta)
{
	return (6.0 - 2.0 * sqrt(1.0 + theta / RO * theta / RO)) / RO;
}


//Quartic Authalic
template <typename T>
T Projections::X_qua_aut(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	//Throw exception
	const T X = R * lonr / RO * cos(lat / RO) / cos(lat / 2.0 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_qua_aut coordinate function, ", "X_qua_aut > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_qua_aut(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 2.0 * R * sin(lat / 2.0 / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_qua_aut coordinate function, ", "Y_qua_aut > MAX_FLOAT: ", Y);

	return Y;
}


//Rectangular Polyconic
template <typename T>
T Projections::X_rpoly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0, A = 0;

	if (fabs(lat1) < MAX_ANGULAR_DIFF)
	{
		A = 0.5 * lonr / RO;
	}

	else
	{
		A = tan(0.5 * lonr / RO * sin(lat1 / RO)) / sin(lat1 / RO);
	}

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = 2.0 * R * A + dx;
	}

	else
	{
		const T E = 2.0 * atan(A * sin(lat / RO));

		//Throw exception
		if (fabs(fabs(lat) - MAX_LAT) < MAX_ANGULAR_DIFF)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in X_rpoly coordinate function, ", "lat = +-90: ", lat);

		X = R / tan(lat / RO) * sin(E) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_rpoly coordinate function, ", "X_rpoly > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_rpoly(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0, A = 0;

	if (fabs(lat1) < MAX_ANGULAR_DIFF)
	{
		A = 0.5 * lonr / RO;
	}

	else
	{
		A = tan(0.5 * lonr / RO * sin(lat1 / RO)) / sin(lat1 / RO);
	}

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = -R * lat1 / RO + dy;
	}

	else
	{
		const T E = 2.0 * atan(A * sin(lat / RO));

		//Throw exception
		if (fabs(fabs(lat) - MAX_LAT) < MAX_ANGULAR_DIFF)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_rpoly coordinate function, ", "lat = +-90: ", lat);

		Y = R * (lat / RO - lat1 / RO + 1.0 / tan(lat / RO) * (1.0 - cos(E))) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_rpoly coordinate function, ", "Y_rpoly > MAX_FLOAT: ", Y);

	return Y;
}


//Sinusoidal
template <typename T>
T Projections::X_sinu(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	
	const T X = R * lonr * cos(lat / RO) / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_sinu coordinate function, ", "X_sinu > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_sinu(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_sinu coordinate function, ", "Y_psinu > MAX_FLOAT: ", Y);

	return Y;
}


//Solovyev Azimuthal
template <typename T>
T Projections::X_solo(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 4;

	//Throw exception
	if (fabs(fabs(A) - 90) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in X_solo coordinate function, ", "lat = +-90: ", lat);

	const T X = 4.0 * R * tan(A / RO) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_solo coordinate function, ", "X_solo > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_solo(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 4;

	//Throw exception
	if (fabs(fabs(A) - 90) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_solo coordinate function, ", "lat = +-90: ", lat);

	const T Y = -4.0 * R * tan(A / RO) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_solo coordinate function, ", "Y_solo > MAX_FLOAT: ", Y);

	return Y;
}


//Stereographic
template <typename T>
T Projections::X_stere(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 2;

	//Throw exception
	if (fabs(A - 90) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in X_stere coordinate function, ", "lat = +-90: ", lat);

	const T X = 2.0 * R * tan(A / RO) * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_stere coordinate function, ", "X_stere > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_stere(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = (90 - lat) / 2;

	//Throw exception
	if (fabs(A - 90) < MAX_ANGULAR_DIFF)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate tan(lat) in Y_stere coordinate function, ", "lat = +-90: ", lat);

	const T Y = -2.0 * R * tan(A / RO) * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_stere coordinate function, ", "Y_stere > MAX_FLOAT: ", Y);

	return Y;
}


//Twilight General Vertical Perspective
template <typename T>
T Projections::X_twi(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.4 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_twi coordinate function, ", "1.0 / A, A = 0:", A);

	const T X = R * 2.4 * R * cos(lat / RO) / A * sin(lonr / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_twi coordinate function, ", "X_twi > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_twi(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = R * (1.4 + sin(lat / RO));

	//Throw exception
	if (fabs(A) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_twi coordinate function, ", "1.0 / A, A = 0:", A);

	const T Y = -R * 2.4 * R * cos(lat / RO) / A * cos(lonr / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_twi coordinate function, ", "Y_twi > MAX_FLOAT: ", Y);

	return Y;
}


//Urmaev V
template <typename T>
T Projections::X_urm5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = asin(0.8 *sin(lat / RO));

	const T X = 2.0 * pow(3, 0.25) / 3.0 * R * lonr / RO * cos(A) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_urm5 coordinate function, ", "X_urm5 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_urm5(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{

	const T A = asin(0.8 *sin(lat / RO));
	const T Y = R * A * (1.0 + 0.414524 / 3.0 * pow((A / RO), 2)) / (0.8 * 2.0 * pow(3, 0.25) / 3) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_urm5 coordinate function, ", "Y_urm5 > MAX_FLOAT: ", Y);

	return Y;
}


//Van der Grinten I
template <typename T>
T Projections::X_vandg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	const T B = fabs(2.0 * lat /180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(B - 1.0) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else
	{
		const T A = 0.5 * fabs(180 / lonr - lonr / 180);
		const T D = A * A;
		const T E = B + C - 1;

		//Throw exception
		if (fabs(E) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_vandg coordinate function, ", "1.0 / E, E = 0:", E);

		const T G = C / E;
		const T P = G * (2.0 / B - 1);
		const T Q = D + G;
		const T S = P * P + D;
		const T TT = G - P * P;
		const T U = A * A * TT * TT - S * (G * G - P * P);

		//Throw exception
		if (U < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(U) in X_vandg coordinate function, ", "U < 0: ", U);

		//Throw exception
		if (fabs(S) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_vandg coordinate function, ", "1.0 / S, S = 0:", S);

		X = R * M_PI * sign(lonr) * (A * TT + sqrt(U)) / S + dx;
	}
		
	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_vandg coordinate function, ", "X_vandg > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_vandg(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	const T B = fabs(2.0 * lat / 180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF))
	{
		Y = R * sign(lat) * M_PI * B / (1.0 + C) + dy;
	}

	else if (fabs(B - 1.0) < MAX_ANGULAR_DIFF)
	{
		Y = R * sign(lat) * M_PI + dy;
	}

	else
	{
		const T A = 0.5 * fabs(180 / lonr - lonr / 180);
		const T D = A * A;
		const T E = B + C - 1;

		//Throw exception
		if (fabs(E) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_vandg coordinate function, ", "1.0 / E, E = 0:", E);

		const T G = C / E;
		const T P = G * (2.0 / B - 1);
		const T Q = D + G;
		const T S = P * P + D;
		const T V = (A * A + 1) * S - Q * Q;

		//Throw exception
		if (V < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(V) in Y_vandg coordinate function, ", "V < 0: ", V);

		//Throw exception
		if (fabs(S) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_vandg coordinate function, ", "1.0 / S, S = 0:", S);

		Y = R * M_PI * sign(lat) * (P * Q - A * sqrt(V)) / S + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_vandg coordinate function, ", "Y_vandg > MAX_FLOAT: ", Y);

	return Y;
}


//Van der Grinten II
template <typename T>
T Projections::X_vandg2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	const T B = fabs(2.0 * lat / 180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(B - 1.0) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else
	{
		const T A = 0.5 * fabs(180 / lonr - lonr / 180);
		const T D = 1.0 + A * A * B * B;

		//Throw exception
		if (fabs(D) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_vandg2 coordinate function, ", "1.0 / D, D = 0:", D);

		const T X1 = (C * sqrt(1.0 + A * A) - A * C * C) / D;

		X = R * M_PI * sign(lonr) * X1 + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_vandg2 coordinate function, ", "X_vandg2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_vandg2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	const T B = fabs(2.0 * lat / 180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF))
	{
		Y = R * sign(lat) * M_PI * B / (1.0 + C) + dy;
	}

	else if (fabs(B - 1.0) < MAX_ANGULAR_DIFF)
	{
		Y = R * sign(lat) * M_PI + dy;
	}

	else
	{
		const T A = 0.5 * fabs(180 / lonr - lonr / 180);
		const T D = 1.0 + A * A * B * B;

		//Throw exception
		if (fabs(D) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_vandg2 coordinate function, ", "1.0 / D, D = 0:", D);

		const T X1 = (C * sqrt(1.0 + A * A) - A * C * C) / D;
		const T E = 1.0 - X1 * X1 - 2.0 * A * X1;

		//Throw exception
		if (E < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(E) in Y_vandg2 coordinate function, ", "E < 0: ", E);

		Y = R * M_PI * sign(lat) * sqrt(E) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_vandg2 coordinate function, ", "Y_vandg2 > MAX_FLOAT: ", Y);

	return Y;
}


//Van der Grinten III
template <typename T>
T Projections::X_vandg3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	const T B = fabs(2.0 * lat / 180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(B - 1.0) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else
	{
		const T A = 0.5 * fabs(180.0 / lonr - lonr / 180.0);
		const T Y1 = B / (1.0 + C);
		const T D = A * A + 1.0 - Y1 * Y1;

		//Throw exception
		if (D < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(D) in X_vandg3 coordinate function, ", "D < 0: ", D);

		X = R * M_PI * sign(lonr) * (sqrt(D) - A) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_vandg3 coordinate function, ", "X_vandg3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_vandg3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	const T B = fabs(2.0 * lat / 180);
	const T C = sqrt(1.0 - B * B);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF))
	{
		Y = R * sign(lat) * M_PI * B / (1.0 + C) + dy;
	}

	else if (fabs(B - 1.0) < MAX_ANGULAR_DIFF)
	{
		Y = R * sign(lat) * M_PI + dy;
	}

	else
	{
		const T Y1 = B / (1.0 + C);
		
		Y = R * M_PI * sign(lat) * Y1 + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_vandg3 coordinate function, ", "Y_vandg3 > MAX_FLOAT: ", Y);

	return Y;
}


//Van der Grinten IV
template <typename T>
T Projections::X_vandg4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	const T B = fabs(2.0 * lat / 180);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		X = R * lonr / RO + dx;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(B - 1.0) < MAX_ANGULAR_DIFF))
	{
		X = dx;
	}

	else
	{
		const T C = 0.5 * (B * (8.0 - B * (2.0 + B * B)) - 5) / (B * B * (B - 1));
		const T C1 = 90 / lonr + lonr / 90;
		const T C2 = C1 * C1 - 4;

		//Throw exception
		if (C2 < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(C2) in X_vandg4 coordinate function, ", "C2 < 0: ", C2);

		const T D = sign(fabs(lonr) - 90) * sqrt(C2);
		const T F1 = (B + C) * (B + C);
		const T F2 = (B + 3.0 * C) * (B + 3.0 * C);
		const T F = F1 * (B * B + C * C * D * D - 1) + (1.0 - B * B) * (B * B * (F2 + 4.0 * C * C) + 12 * B * C * C * C + 4.0 * C * C * C * C);

		//Throw exception
		if (F < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(F) in X_vandg4 coordinate function, ", "F < 0: ", F);

		const T G = 4.0 * F1 + D * D;

		//Throw exception
		if (fabs(G) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_vandg4 coordinate function, ", "1.0 / G, G = 0:", G);

		const T X1 = (D * (F1 + C * C - 1) + 2.0 * sqrt(F)) / G;
		
		X = R * M_PI * 0.5 * sign(lonr) * X1 + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_vandg4 coordinate function, ", "X_vandg4 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_vandg4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	const T B = fabs(2.0 * lat / 180);

	if (fabs(lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else if ((fabs(lonr) < MAX_ANGULAR_DIFF) || (fabs(B - 1.0) < MAX_ANGULAR_DIFF))
	{
		Y = R * lat / RO + dy;
	}

	else
	{
		const T C = 0.5 * (B * (8.0 - B * (2.0 + B * B)) - 5) / (B * B * (B - 1));
		const T C1 = 90 / lonr + lonr / 90;
		const T C2 = C1 * C1 - 4;

		//Throw exception
		if (C2 < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(C2) in X_vandg4 coordinate function, ", "C2 < 0: ", C2);

		const T D = sign(fabs(lonr) - 90) * sqrt(C2);
		const T F1 = (B + C) * (B + C);
		const T F2 = (B + 3.0 * C) * (B + 3.0 * C);
		const T F = F1 * (B * B + C * C * D * D - 1) + (1.0 - B * B) * (B * B * (F2 + 4.0 * C * C) + 12 * B * C * C * C + 4.0 * C * C * C * C);

		//Throw exception
		if (F < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(F) in X_vandg4 coordinate function, ", "F < 0: ", F);

		const T G = 4.0 * F1 + D * D;

		//Throw exception
		if (fabs(G) < MIN_FLOAT)
			throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_vandg4 coordinate function, ", "1.0 / G, G = 0:", G);

		const T X1 = (D * (F1 + C * C - 1) + 2.0 * sqrt(F)) / G;
		const T H = 1.0 + D * fabs(X1) - X1 * X1;

		//Throw exception
		if (H < 0)
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(H) in Y_vandg4 coordinate function, ", "H < 0: ", H);


		Y = R * M_PI * 0.5 * sign(lat) * sqrt(H) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_vandg4 coordinate function, ", "Y_vandg4 > MAX_FLOAT: ", Y);

	return Y;
}


//Wagner I
template <typename T>
T Projections::X_wag1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T M = 2.0 * sqrt(sqrt(3)) / 3.0;
	const T N = 0.5 * sqrt(3);
	const T A = N * sin(lat / RO);

	if (fabs(A) > 1.0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate asin(A) in X_wag1 coordinate function, ", "abs(A) > 1: ", A);

	const T theta = asin(A) * RO;

	const T X = R * M * lonr / RO * cos(theta / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag1 coordinate function, ", "X_wag1 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T M = 2.0 * sqrt(sqrt(3)) / 3.0;
	const T N = 0.5 * sqrt(3);
	const T A = N * sin(lat / RO);

	if (fabs(A) > 1.0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate asin(A) in Y_wag1 coordinate function, ", "abs(A) > 1: ", A);

	const T theta = asin(A) * RO;

	const T Y = R * 3.0 * theta / RO * M * N /2.0 + dy;		//Error in Evenden G, I: Cartographic Projection Procedures, page 7

	//Throw exception
	if (fabs(Y) > MAX_FLOAT)
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag1 coordinate function, ", "Y_wag1 > MAX_FLOAT: ", Y);

	return Y;
}


//Wagner II
template <typename T>
T Projections::X_wag2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 0.92483 * R * lonr / RO * cos(asin(0.88022 * sin(0.8855 * lat / RO))) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag2 coordinate function, ", "X_wag2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 1.38725 * R * asin(0.88022 * sin(0.8855 * lat / RO)) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag2 coordinate function, ", "Y_wag2 > MAX_FLOAT: ", Y);

	return Y;
}


//Wagner III
template <typename T>
T Projections::X_wag3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = lat1 / RO;
	const T B = cos(2.0 * A / 3);

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_wag3 coordinate function, ", "1.0 / A, A = 0:", A);


	const T X = R * lonr / RO * (cos(A) / B) * cos(2.0 * lat / 3.0 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag3 coordinate function, ", "X_wag3 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag3(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag3 coordinate function, ", "Y_wag3 > MAX_FLOAT: ", Y);

	return Y;
}


//Wagner IV
template <typename T>
T Projections::X_wag4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = 0.5 * lat;
	const T theta = NewtonRaphson(FTheta_wag4<T>, FThetaDer_wag4<T>, lat, theta0);

	const T X = 0.8631 * R * lonr / RO * cos(theta / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag4 coordinate function, ", "X_wag4 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag4(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = 0.5 * lat;
	const T theta = NewtonRaphson(FTheta_wag4<T>, FThetaDer_wag4<T>, lat, theta0);

	const T Y = 1.5654 * R * sin(theta / RO) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag4 coordinate function, ", "Y_wag4 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_wag4(const T lat, const T theta)
{
	return 2.0 * theta / RO + sin(2.0 * theta / RO) - (4.0 * M_PI + 3.0 * sqrt(3)) / 6 * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_wag4(const T lat, const T theta)
{
	return ( 2.0 + 2.0 * cos(2.0 * theta / RO)) / RO;
}


//Wagner VI
template <typename T>
T Projections::X_wag6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = 1.0 - 3.0 * pow((lat / (M_PI * RO)), 2);

	//Throw exception
	if (A < 0)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(A) in X_wag6 coordinate function, ", "A < 0: ", A);

	const T X = 1.89490 * R * (-0.5 + sqrt(A)) * lonr / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag6 coordinate function, ", "X_wag6 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag6(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = 0.94745 * R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag6 coordinate function, ", "Y_wag6 > MAX_FLOAT: ", Y);

	return Y;
}


//Wagner VII
template <typename T>
T Projections::X_wag7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = pow((0.90631 * sin(lat / RO)), 2);

	//Throw exception
	if (A > 1)
		throw MathInvalidArgumentException <T>("MathInvalidArgumentException: can not evaluate sqrt(A) in X_wag7 coordinate function, ", "A < 0: ", A);

	const T B = sqrt(1.0 - A);
	const T C = 1.0 + B * cos(lonr / 3.0 / RO);

	//Throw exception
	if (fabs(C) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_wag7 coordinate function, ", "1.0 / C, C = 0:", C);

	const T D = 2.0 / C;

	const T X = 2.66723 * R * B * sqrt(D) * sin(lonr / 3.0 / RO) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wag7 coordinate function, ", "X_wag7 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wag7(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = 0.90631 * sin(lat / RO);
	const T B = 1.0 + sqrt(1.0 - pow((A), 2)) * cos(lonr / 3.0 / RO);

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate Y_wag7 coordinate function, ", "1.0 / B, B = 0:", B);

	const T Y = 1.24104 * R * A * sqrt(2.0 / B) + dy;
	
	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wag7 coordinate function, ", "Y_wag7 > MAX_FLOAT: ", Y);

	return Y;
}


//Werner-Staab
template <typename T>
T Projections::X_wer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T X = 0;

	if (fabs (90 - lat) < MAX_ANGULAR_DIFF)
	{
		X = dx;
	}

	else
	{
		X = R * (90 - lat) / RO * sin(lonr * cos(lat / RO) / (90 - lat)) + dx;
	}

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wer coordinate function, ", "X_wer > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wer(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);
	T Y = 0;

	if (fabs(90 - lat) < MAX_ANGULAR_DIFF)
	{
		Y = dy;
	}

	else
	{
		Y = R * (lat - 90) / RO * cos(lonr * cos(lat / RO) / (90 - lat)) + dy;
	}

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wer coordinate function, ", "Y_wer > MAX_FLOAT: ", Y);

	return Y;
}


//Werenskiold I
template <typename T>
T Projections::X_weren(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T A = asin(5 * sqrt(2) / 8 * sin(lat / RO));
	const T B = cos(A / 3);

	//Throw exception
	if (fabs(B) < MIN_FLOAT)
		throw  MathZeroDevisionException <T>("MathZeroDevisionException: can not evaluate X_weren coordinate function, ", "1.0 / B, B = 0:", B);

	const T X = R * lonr / RO * cos(A) / B + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_weren coordinate function, ", "X_weren > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_weren(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T A = asin(5 * sqrt(2) / 8 * sin(lat / RO));
	const T Y = 2.0 * M_PI * sqrt(2) * R * sin(A / 3) + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_weren coordinate function, ", "Y_weren > MAX_FLOAT: ", Y);

	return Y;
}


//Winkel I
template <typename T>
T Projections::X_wink1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T X = 0.5 * R * lonr * (cos(lat1 / RO) + cos(lat / RO)) / RO + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wink1 coordinate function, ", "X_wink1 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wink1(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T Y = R * lat / RO + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wink1 coordinate function, ", "Y_wink1 > MAX_FLOAT: ", Y);

	return Y;
}


//Winkel II
template <typename T>
T Projections::X_wink2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T lonr = CartTransformation::redLon0(lon, lon0);

	const T theta0 = 0.9 * lat;
	const T theta = NewtonRaphson(FTheta_wink2<T>, FThetaDer_wink2<T>, lat, theta0);

	const T X = 0.5 * R * lonr / RO * (cos(theta / RO) + cos(lat1 / RO)) + dx;

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wink2 coordinate function, ", "X_wink2 > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wink2(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	const T theta0 = 0.9 * lat;
	const T theta = NewtonRaphson(FTheta_wink2<T>, FThetaDer_wink2<T>, lat, theta0);
	
	const T Y = M_PI * R * (sin(theta / RO) + 2.0 * lat / M_PI / RO) / 4.0 + dy;

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wink2 coordinate function, ", "Y_wink2 > MAX_FLOAT: ", Y);

	return Y;
}

template <typename T>
T Projections::FTheta_wink2(const T lat, const T theta)
{
	return 2.0 * theta / RO + sin(2.0 * theta / RO) - M_PI * sin(lat / RO);
}

template <typename T>
T Projections::FThetaDer_wink2(const T lat, const T theta)
{
	return (2.0 + 2.0 * cos(2.0 * theta / RO)) / RO;
}


//Winkel Tripel
template <typename T>
T Projections::X_wintri(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	//Aitoff projection
	const T X1 = X_aitoff(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	
	//Equidistant conic
	const T X2 = X_eqc(R, lat1, lat2, lat, lon, lon0, dx, dy, c);
	
	//Average of both projections
	const T X = 0.5 * (X1 + X2);

	//Throw exception
	if (fabs(X) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate X_wintri coordinate function, ", "X_wintri > MAX_FLOAT: ", X);

	return X;
}

template <typename T>
T Projections::Y_wintri(const T R, const T lat1, const T lat2, const T lat, const T lon, const T lon0, const T dx, const T dy, const T c)
{
	//Aitoff projection
	const T Y1 = Y_aitoff(R, lat1, lat2, lat, lon, lon0, dx, dy, c);

	//Equidistant conic
	const T Y2 = Y_eqc(R, lat1, lat2, lat, lon, lon0, dx, dy, c);

	//Average of both projections
	const T Y = 0.5 * (Y1 + Y2);

	//Throw exception
	if (fabs(Y) > MAX_FLOAT )
		throw MathOverflowException <T>("MathOverflowException: can not evaluate Y_wintri coordinate function, ", "Y_wintri > MAX_FLOAT: ", Y);

	return Y;
}


template <typename T, typename FTheta, typename FThetaDer>
T Projections::NewtonRaphson(FTheta ftheta, FThetaDer fthetader, const T lat, const T theta0, const unsigned int MAX_ITERATIONS, const T MAX_DIFF)
{
	//Solve theta = f(theta, lat) by the Newton-Raphson method
	//Used in several pseudocylindrical projections
	unsigned int iterations = 0;
	T theta = theta0;

	//Apply Newton-Raphson method
	do {
		//Compute F(theta) and F'(theta)
		const T ft = ftheta(lat, theta);
		const T ft_der = fthetader(lat, theta);

		//Newton-Raphson step
		T theta_n = theta;
		if (fabs(ft_der) > MIN_FLOAT) {
			theta_n = theta - ft / ft_der;
		}

		//Test the terminal condition
		if (fabs(theta_n - theta) < MAX_DIFF) {
			break;
		}

		//Assign new theta
		theta = theta_n;

	} while (++iterations <= MAX_ITERATIONS);

	return theta;
}

#endif