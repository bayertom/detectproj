// Description: Connvert Projection to Proj4 output string

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

#ifndef ProjectionToProj42_HPP
#define ProjectionToProj42_HPP

#include <typeinfo>

#include "libalgo/source/structures/projection2/ProjectionCylindrical.h"
#include "libalgo/source/structures/projection2/ProjectionConic.h"
#include "libalgo/source/structures/projection2/ProjectionPseudoConic.h"


template <typename T>
std::string ProjectionToProj4::ProjectionToProj4String(const std::shared_ptr <Projection <T> > proj)
{
	//Convert projection definition to Proj.4 string
	//Return empty string, if the projection is not supported by Proj.4
	std::string proj4_string = ""; 
	TProjNamesMap proj_names_list;

	//A projection was given
	if (proj != NULL)
	{
		//Initialize projections
		init(proj_names_list);

		//Get projection name
		std::string proj_name = proj->getName();

		//Find projection name in the list
		TProjNamesMap::iterator i_proj_names_list = proj_names_list.find(std::string(proj_name));

		//Projection name found in the list
		if (i_proj_names_list != proj_names_list.end())
		{
			//Get its Proj.4 definition
			std::string proj4_name = i_proj_names_list->second;

			//Is the projection supported by Proj.4
			const std::string not_supported = "not_supported";
			if (proj4_name.compare(not_supported) != 0)
			{
				//Add commnand to the string
				proj4_string += "proj";

				//Get projection properties
				const T R = proj->getR();
				const Point3DGeographic <T> pole = proj->getCartPole();
				const T lat1 = proj->getLat1();
				const T lat2 = proj->getLat2();
				const T lon0 = proj->getLon0();
				const T c = proj->getC();
				const T dx = proj->getDx();
				const T dy = proj->getDy();

				//************************* Process R
				proj4_string += " +R=" + to_string(R);

				//************************* Process projection
				proj4_string += " +proj=" + proj4_name;

				//************************* Process projection aspect
				const T latp = pole.getLat();
				const T lonp = pole.getLon();

				if ((fabs(latp - MAX_LAT) > ANGLE_ROUND_ERROR) && (fabs(lonp - MAX_LON) > ANGLE_ROUND_ERROR))
				{
					proj4_string += " +o_lat_p=" + to_string(latp);
					proj4_string += " +o_lon_b=" + to_string(lonp);
				}

				//************************* Process lat0
				//Conic projections, use lat_1, lat_2
				if (typeid(*(proj)) == typeid(ProjectionConic<T>))
				{
					//Equidistant conic or Albers: 2 different standard parallels lat_1, lat_2 in detectproj
					if((proj_name.compare(detectprojNames[p_eqdc3]) == 0) || (proj_name.compare (detectprojNames[p_aea]) == 0))
					{
						proj4_string += " +lat_1=" + to_string(lat1);
						proj4_string += " +lat_2=" + to_string(lat2);
					}

					// Other conic: 1 standard parallel lat_0 in detectproj
					else
					{
						proj4_string += " +lat_1=" + to_string(lat1);
						proj4_string += " +lat_2=" + to_string(lat1);
					}
				}

				//Cylindrical projections, use lat_ts
				else if (typeid(*(proj)) == typeid(ProjectionCylindrical<T>))
				{
					proj4_string += " +lat_ts=" + to_string(lat1);
				}

				//Pseudoconic projection, use lat_1
				else if (typeid(*(proj)) == typeid(ProjectionPseudoConic<T>))
				{
					proj4_string += " +lat_1=" + to_string(lat1);
				}

				//Werner staab projection = Bonne, lat_0 = 90
				else if (proj_name.compare( detectprojNames[p_wer]) == 0)
				{
					proj4_string += " +lat_0=90";
				}

				//Other projections use current lat_0
				else
				{
					proj4_string += " +lat_0=" + to_string(lat1);
				}

				//************************* Process lon0
				if (lon0 != 0)
				{
					proj4_string += " +lon_0=" + to_string(lon0);
				}

				//************************* Process c (perspective projections)
				if (proj_name.compare (detectprojNames[p_nsper]) == 0)
				{
					const T h = R * c;

					proj4_string += " +h=" + to_string(h);
				}

				//************************* Shifs dx, dy corresponding to false northing and false easting
				proj4_string += " +x_0=" + to_string(dx); 
				proj4_string += " +y_0=" + to_string(dy);
			}
		}
	}

	return proj4_string;
}


#endif