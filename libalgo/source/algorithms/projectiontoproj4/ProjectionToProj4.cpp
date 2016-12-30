// Description: Connvert Projection to Proj4 output string

// Copyright (c) 2010 - 2015
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


#include "ProjectionToProj4.h"


//List of projections supported by detectproj (some of them have different ID or unsupported by Proj.4)
const char * detectprojNames[] =
{
	"aea",
	"aeqd",
	"aitoff",
	"apian",
	"apiel",		/*Apian elliptic, not supported by Proj. 4 */	

	"armad",		/*Armadillo (orthoapsidal), not supported by Proj. 4*/
	"bacon",
	"bonne",
	"breus",		/*Breusign*, not supported by Proj. 4*/
	"cea",

	"clar",			/*clar, Clark far side perspective, not supported by Proj. 4*/
	"collg",
	"denoy",
	"eck1",
	"eck2",

	"eck3",
	"eck4",
	"eck5",
	"eck6",
	"eqc",

	"eqdc",
	"eqdc2",		/*equidistant conic, pole = point*/
	"eqdc3",		/*equidistant conic, two standard parallels lat1, lat2*/
	"fahey",
	"fouc",

	"fouc_s",
	"fourn",		/*Fournier, not supported by Proj. 4*/
	"fourn2",		/*Fournier 2, not supported by Proj. 4*/
	"fsper",		/*far side perspective, not supported by Proj. 4*/
	"gall",

	"gins8",
	"gnom",
	"hammer",
	"hataea",
	"hire",			/*La Hire far side perspective, not supported by Proj. 4*/

	"jam",			/*James far side perspective, not supported by Proj. 4*/
	"kav5",
	"kav7",
	"laea",
	"larr",

	"lask",

	"lcc",
	"krovak",
	"leac",
	"leac2",		/*conic equal area, pole = point*/
	"mbt_s",

	"mbtfpq",
	"mbtfps",
	"merc",
	"moll",
	"nel_h",

	"nsper",
	"nicol",
	"ortel",
	"orho",
	"parab",

	"pers",			/*perspective, not supported by Proj. 4*/
	"poly",
	"putp1",
	"putp3",
	"putp3p",

	"putp4",
	"putp4p",
	"putp5",
	"putp5p",
	"qua_aut",

	"sinu",
	"sinu2",		/*modified sinusoidal projection, one standard parallel lat0*/
	"solo",			/*Solovjev, not supported by Proj. 4*/
	"stere",
	"twi",			/*Twilight far side perspective, not supported by Proj. 4*/

	"urm5",
	"wag2",
	"wag3",
	"wag4",
	"wag6",

	"wag7",
	"wer",			/*Werner defined as Bonne, lat0=90*/
	"weren",
	"wink1",
	"wink2",

	"wintri"
};


//List of corresponding projections in Proj4: analogous, different ID, or unsupported
const char * proj4Names[] =
{
	"aea",
	"aeqd",
	"aitoff",
	"apian",
	"not_supported",	/*apiel, Apian elliptic*/

	"not_supported",	/*armad, Armadillo (orthoapsidal)*/
	"bacon",
	"bonne",
	"not_supported",	/*breus, Breusign*/
	"cea",

	"not_supported",	/*clar, Clark far side perspective*/
	"collg",
	"denoy",
	"eck1",
	"eck2",

	"eck3",
	"eck4",
	"eck5",
	"eck6",
	"eqc",

	"eqdc",
	"eqdc",			/*eqdc, equidistant conic, pole = point*/
	"eqdc",			/*eqdc, equidistant conic, two standard parallels lat1, lat2*/
	"fahey",
	"fouc",

	"fouc_s",
	"not_supported",	/*fourn, Fournier */
	"not_supported",	/*fourn2, Fournier 2*/
	"not_supported", 	/*fsper, Far side perspective*/
	"gall",

	"gins8",
	"gnom",
	"hammer",
	"hataea",
	"not_supported",	/*hire, La Hire far side perspective*/

	"not_supported",	/*jam, James far side perspective*/
	"kav5",
	"kav7",
	"laea",
	"larr",

	"lask",

	"lcc",
	"krovak",
	"leac",
	"leac",			/*leac, conic equal area, pole = point*/
	"mbt_s",

	"mbtfpq",
	"mbtfps",
	"merc",
	"moll",
	"nel_h",

	"nsper",
	"nicol",
	"ortel",
	"orho",
	"parab",

	"not_supported",	/*pers, Perspective*/
	"poly",
	"putp1",
	"putp3",
	"putp3p",

	"putp4",
	"putp4p",
	"putp5",
	"putp5p",
	"qua_aut",

	"sinu",
	"sinu",			/*modified sinusoidal projection, one standard parallel lat0*/
	"not_supported",	/*solo, Solovjev*/
	"stere",
	"not_supported",	/*twi, Twilight far side perspective*/

	"urm5",
	"wag2",
	"wag3",
	"wag4",
	"wag6",

	"wag7",
	"bonne",		/*Werner defined as Bonne, lat0=90*/
	"weren",
	"wink1",
	"wink2",

	"wintri"
};



void ProjectionToProj4::init( TProjNamesMap &proj_names_list )
{
	//Initialize map of corresponding projections <detectproj, proj4>
	const unsigned int n1 = (sizeof(detectprojNames) / sizeof(*detectprojNames));
	const unsigned int n2 = (sizeof(proj4Names) / sizeof(*proj4Names));

	//Bad detectproj->Proj.4 conversion rules definition
	if (n1 != n2)
		throw BadDataException("BadDataException: different amount of projections for the conversion between detectproj and Proj.4. ", "Conversion to Proj.4 cancelled...");

	//Assign corresponding projections
	proj_names_list[detectprojNames[p_aea]] = proj4Names[p_aea];
	proj_names_list[detectprojNames[p_aeqd]] = proj4Names[p_aeqd];
	proj_names_list[detectprojNames[p_aitoff]] = proj4Names[p_aitoff];
	proj_names_list[detectprojNames[p_apian]] = proj4Names[p_apian];
	proj_names_list[detectprojNames[p_apiel]] = proj4Names[p_apiel];

	proj_names_list[detectprojNames[p_armad]] = proj4Names[p_armad];
	proj_names_list[detectprojNames[p_bacon]] = proj4Names[p_bacon];
	proj_names_list[detectprojNames[p_bonne]] = proj4Names[p_bonne];
	proj_names_list[detectprojNames[p_breus]] = proj4Names[p_breus];
	proj_names_list[detectprojNames[p_cea]] = proj4Names[p_cea];

	proj_names_list[detectprojNames[p_clar]] = proj4Names[p_clar];
	proj_names_list[detectprojNames[p_collg]] = proj4Names[p_collg];
	proj_names_list[detectprojNames[p_denoy]] = proj4Names[p_denoy];
	proj_names_list[detectprojNames[p_eck1]] = proj4Names[p_eck1];
	proj_names_list[detectprojNames[p_eck2]] = proj4Names[p_eck2];

	proj_names_list[detectprojNames[p_eck3]] = proj4Names[p_eck3];
	proj_names_list[detectprojNames[p_eck4]] = proj4Names[p_eck4];
	proj_names_list[detectprojNames[p_eck5]] = proj4Names[p_eck5];
	proj_names_list[detectprojNames[p_eck6]] = proj4Names[p_eck6];
	proj_names_list[detectprojNames[p_eqc]] = proj4Names[p_eqc];

	proj_names_list[detectprojNames[p_eqdc]] = proj4Names[p_eqdc];
	proj_names_list[detectprojNames[p_eqdc2]] = proj4Names[p_eqdc2];
	proj_names_list[detectprojNames[p_eqdc3]] = proj4Names[p_eqdc3];
	proj_names_list[detectprojNames[p_fahey]] = proj4Names[p_fahey];
	proj_names_list[detectprojNames[p_fouc]] = proj4Names[p_fouc];

	proj_names_list[detectprojNames[p_fouc_s]] = proj4Names[p_fouc_s];
	proj_names_list[detectprojNames[p_fourn]] = proj4Names[p_fourn];
	proj_names_list[detectprojNames[p_fourn2]] = proj4Names[p_fourn2];
	proj_names_list[detectprojNames[p_fsper]] = proj4Names[p_fsper];
	proj_names_list[detectprojNames[p_gall]] = proj4Names[p_gall];

	proj_names_list[detectprojNames[p_gins8]] = proj4Names[p_gins8];
	proj_names_list[detectprojNames[p_gnom]] = proj4Names[p_gnom];
	proj_names_list[detectprojNames[p_hammer]] = proj4Names[p_hammer];
	proj_names_list[detectprojNames[p_hataea]] = proj4Names[p_hataea];
	proj_names_list[detectprojNames[p_hire]] = proj4Names[p_hire];

	proj_names_list[detectprojNames[p_jam]] = proj4Names[p_jam];
	proj_names_list[detectprojNames[p_kav5]] = proj4Names[p_kav5];
	proj_names_list[detectprojNames[p_kav7]] = proj4Names[p_kav7];
	proj_names_list[detectprojNames[p_laea]] = proj4Names[p_laea];
	proj_names_list[detectprojNames[p_larr]] = proj4Names[p_larr];

	proj_names_list[detectprojNames[p_lask]] = proj4Names[p_lask];

	proj_names_list[detectprojNames[p_lcc]] = proj4Names[p_lcc];
	proj_names_list[detectprojNames[p_krovak]] = proj4Names[p_krovak];
	proj_names_list[detectprojNames[p_leac]] = proj4Names[p_leac];
	proj_names_list[detectprojNames[p_leac2]] = proj4Names[p_leac2];
	proj_names_list[detectprojNames[p_mbt_s]] = proj4Names[p_mbt_s];

	proj_names_list[detectprojNames[p_mbtfpq]] = proj4Names[p_mbtfpq];
	proj_names_list[detectprojNames[p_mbtfps]] = proj4Names[p_mbtfps];
	proj_names_list[detectprojNames[p_merc]] = proj4Names[p_merc];
	proj_names_list[detectprojNames[p_moll]] = proj4Names[p_moll];
	proj_names_list[detectprojNames[p_nel_h]] = proj4Names[p_nel_h];

	proj_names_list[detectprojNames[p_nsper]] = proj4Names[p_nsper];
	proj_names_list[detectprojNames[p_nicol]] = proj4Names[p_nicol];
	proj_names_list[detectprojNames[p_ortel]] = proj4Names[p_ortel];
	proj_names_list[detectprojNames[p_orho]] = proj4Names[p_orho];
	proj_names_list[detectprojNames[p_parab]] = proj4Names[p_parab];

	proj_names_list[detectprojNames[p_pers]] = proj4Names[p_pers];
	proj_names_list[detectprojNames[p_poly]] = proj4Names[p_poly];
	proj_names_list[detectprojNames[p_putp1]] = proj4Names[p_putp1];
	proj_names_list[detectprojNames[p_putp3]] = proj4Names[p_putp3];
	proj_names_list[detectprojNames[p_putp3p]] = proj4Names[p_putp3p];

	proj_names_list[detectprojNames[p_putp4]] = proj4Names[p_putp4];
	proj_names_list[detectprojNames[p_putp4p]] = proj4Names[p_putp4p];
	proj_names_list[detectprojNames[p_putp5]] = proj4Names[p_putp5];
	proj_names_list[detectprojNames[p_putp5p]] = proj4Names[p_putp5p];
	proj_names_list[detectprojNames[p_qua_aut]] = proj4Names[p_qua_aut];

	proj_names_list[detectprojNames[p_sinu]] = proj4Names[p_sinu];
	proj_names_list[detectprojNames[p_sinu2]] = proj4Names[p_sinu2];
	proj_names_list[detectprojNames[p_solo]] = proj4Names[p_solo];
	proj_names_list[detectprojNames[p_stere]] = proj4Names[p_stere];
	proj_names_list[detectprojNames[p_twi]] = proj4Names[p_twi];

	proj_names_list[detectprojNames[p_urm5]] = proj4Names[p_urm5];
	proj_names_list[detectprojNames[p_wag2]] = proj4Names[p_wag2];
	proj_names_list[detectprojNames[p_wag3]] = proj4Names[p_wag3];
	proj_names_list[detectprojNames[p_wag4]] = proj4Names[p_wag4];
	proj_names_list[detectprojNames[p_wqg6]] = proj4Names[p_wqg6];

	proj_names_list[detectprojNames[p_wag7]] = proj4Names[p_wag7];
	proj_names_list[detectprojNames[p_wer]] = proj4Names[p_wer];
	proj_names_list[detectprojNames[p_weren]] = proj4Names[p_weren];
	proj_names_list[detectprojNames[p_wink1]] = proj4Names[p_wink1];
	proj_names_list[detectprojNames[p_wink2]] = proj4Names[p_wink2];

	proj_names_list[detectprojNames[p_wintri]] = proj4Names[p_wintri];
}
