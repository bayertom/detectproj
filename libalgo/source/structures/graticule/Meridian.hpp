// Description: A meridian of the longitude given by the list of points

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


#ifndef Meridian_HPP
#define Meridian_HPP

#include "libalgo/source/algorithms/round/Round.h"

#include "libalgo/source/algorithms/carttransformation2/CartTransformation.h"

template <typename T>
void Meridian<T>::createMeridian()
{
	//Create meridian

 	//Set start value of the latitude as a multiplier of dlat 
	const T lat_start = Round::roundToMultipleFloor(lat_interval.min_value, dlat) + dlat;
	const T lat_end = Round::roundToMultipleCeil(lat_interval.max_value, dlat) - dlat;

	//Add first point (lower bound of the interval)
	lats.push_back(std::max(std::min(lat_interval.min_value + lat_min_shift, MAX_LAT), MIN_LAT));

	//Add intermediate points
	for (T lat_point = lat_start; lat_point <= lat_end; lat_point += dlat)
	{
		lats.push_back(lat_point);
	}

	//Add last point (upper bound of the interval)
	lats.push_back(std::max(std::min(lat_interval.max_value - lat_max_shift, MAX_LAT), MIN_LAT));
}


template <typename T>
template <typename Point>
void Meridian<T>::project(const std::shared_ptr <Projection <T> > proj, const T alpha, TVector <Point> &mer)
{
	//Project meridian
	for (const auto lat : lats)
	{
		T X = 0.0, Y = 0.0;

		//Project point
		try
		{
			CartTransformation::latLontoXY(lat, lon, proj, alpha, X, Y);
		}

		//Throw exception
		catch (MathException <T> &error)
		{
			//Throw new math error
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: error in coordinate function.", "Can not compute meridian points, lat = ", lat);
		}

		//Add projected point to the list
		Point p_proj(X, Y);
		mer.push_back(p_proj);
	}
}


template <typename T>
void Meridian <T>::print(std::ostream &output)
{
	for (const auto lat : lats)
	{
		output << "[" << lat << ", " << lon << "]  ";
	}

	output << '\n' << '\n';
}

#endif
