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


#ifndef ProjectionToProj4_H
#define ProjectionToProj4_H

#include <map>
#include <cstring>
#include <string>
#include <sstream> 


//Extern declarations
extern const char * proj4Names[];
extern const char * detectprojNames[];


//Forward declarations
template <typename T>
class Projection;


//Projections identificators, corresponding to detectproj
enum projections
{
	p_aea = 0,
	p_aeqd,
	p_aitoff,
	p_apian,
	p_apiel,

	p_armad,
	p_bacon,
	p_bonne,
	p_breus,
	p_cea,

	p_clar,
	p_collg,
	p_denoy,
	p_eck1,
	p_eck2,

	p_eck3,
	p_eck4,
	p_eck5,
	p_eck6,
	p_eqc,

	p_eqdc,
	p_eqdc2,
	p_eqdc3,
	p_fahey,
	p_fouc,

	p_fouc_s,
	p_fourn,
	p_fourn2,
	p_fsper,
	p_gall,

	p_gins8,
	p_gnom,
	p_hammer,
	p_hataea,
	p_hire,

	p_jam,
	p_kav5,
	p_kav7,
	p_laea,
	p_larr,

	p_lask,

	p_lcc,
	p_krovak,
	p_leac,
	p_leac2,
	p_mbt_s,

	p_mbtfpq,
	p_mbtfps,
	p_merc,
	p_moll,
	p_nel_h,

	p_nsper,
	p_nicol,
	p_ortel,
	p_orho,
	p_parab,

	p_pers,
	p_poly,
	p_putp1,
	p_putp3,
	p_putp3p,

	p_putp4,
	p_putp4p,
	p_putp5,
	p_putp5p,
	p_qua_aut,

	p_sinu,
	p_sinu2,
	p_solo,
	p_stere,
	p_twi,

	p_urm5,
	p_wag2,
	p_wag3,
	p_wag4,
	p_wqg6,

	p_wag7,
	p_wer,
	p_weren,
	p_wink1,
	p_wink2,

	p_wintri
};


//Comparator for projection names represented by strings
struct compProjNamesMap
{
	bool operator() (const std::string & a, const std::string & b) const
	{
		return std::strcmp(a.c_str(), b.c_str()) < 0;
	}
};


//Map of projection names <detectproj, proj.4>
typedef std::map <std::string, std::string, compProjNamesMap> TProjNamesMap;


//Connvert Projection <T> to Proj.4 output string
class ProjectionToProj4
{

	public:
		template <typename T>
		static std::string ProjectionToProj4String (const Projection <T> *proj);

	private:
		static void init( TProjNamesMap &proj_names_list );
};


template <typename T>
std::string to_string(T const & val) 
{
    //Convert number to string (full support to_string() in C++ X11)
    std::stringstream s_string;
    s_string << val;
    return s_string.str();
}


#include "ProjectionToProj4.hpp"

#endif
