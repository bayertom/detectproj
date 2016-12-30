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


#ifndef Parallel_H
#define Parallel_H

#include "libalgo/source/types/TInterval.h"
#include "libalgo/source/types/TVector.h"


template <typename T>
class Parallel
{
	private:
		 T lat;						//Latitude of the parallel
		 TInterval <T> lon_interval;			//Longitude interval of the parallel
		 T dlon;					//Longitude step
		 T lon_min_shift;				//Shift of the initial parallel point: [lat, lon_min + lon_min_shift], default = 0
		 T lon_max_shift;				//Shift of the last parallel point: [lat, lon_max - lon_max_shift], default = 0
		TVector <T> lons;				//Longitude of the parallel points

	public:
		
		Parallel(const T lat_) : lat(lat_), lon_interval{ MIN_LON, MAX_LON }, dlon(10), lon_min_shift(0), lon_max_shift(0) { createParallel(); };
		Parallel(const T lat_, const TInterval <T> lon_interval_, const T dlon_, const T lon_min_shift_, const T lon_max_shift_) :
			lat(lat_), lon_interval(lon_interval_), dlon(dlon_), lon_min_shift(lon_min_shift_), lon_max_shift(lon_max_shift_) {createParallel();};
		
		T getLat() const { return lat; }
		TInterval <T> & getLonInterval() const { return lon_interval; }
		T getDLon() const { return dlon; }
		T getLonMinShift() const { return lon_min_shift; }
		T getLonMaxShift() const { return lon_max_shift; }
		TVector <T> & getLons() const { return lons; }

		template <typename Point>
		void project(const std::shared_ptr <Projection <T> > proj, const T alpha, TVector <Point> &par);

		void print(std::ostream &output = std::cout);

	private:
		void createParallel();
};

#include "Parallel.hpp"

#endif