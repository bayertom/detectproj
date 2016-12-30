// Description: A parallel of the latitude given by the list of points

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


#ifndef Parallel_HPP
#define Parallel_HPP

#include "libalgo/source/algorithms/round/Round.h"

#include "libalgo/source/algorithms/carttransformation2/CartTransformation.h"


template <typename T>
void Parallel<T>::createParallel()
{
	//Create parallel

	//Set start value of the latitude as a multiplier of dlon 
	const T lon_start = Round::roundToMultipleFloor(lon_interval.min, dlon) + dlon;
	const T lon_end = Round::roundToMultipleCeil(lon_interval.max, dlon) - dlon;

	//Add first point (lower bound of the interval)
	lons.push_back(std::max(std::min(lon_interval.min + lon_min_shift, MAX_LON), MIN_LON));

	//Add intermediate points
	for (T lon_point = lon_start; lon_point <= lon_end; lon_point += dlon)
	{
		lons.push_back(lon_point);
	}

	//Add last point (upper bound of the interval)
	lons.push_back(std::max(std::min(lon_interval.max - lon_max_shift, MAX_LON), MIN_LON));
}


template <typename T>
template <typename Point>
void Parallel<T>::project(const std::shared_ptr <Projection <T> > proj, const T alpha, TVector <Point> &par)
{
	//Project parallel
	for (const auto lon : lons)
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
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: error in coordinate function.", "Can not compute parallel points, lon = ", lon);
		}

		//Add projected point to the list
		Point p_proj(X, Y);
		par.push_back(p_proj);
	}
}


template <typename T>
void Parallel<T>::print(std::ostream &output)
{
	for (const auto lon : lons)
	{
		output << "[" << lat << ", " << lon << "]  ";
	}

	output << '\n' << '\n';
}

#endif