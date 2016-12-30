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


#ifndef Meridian_H
#define Meridian_H

#include <memory>

#include "libalgo/source/types/TInterval.h"
#include "libalgo/source/types/TVector.h"

#include "libalgo/source/structures/projection2/Projection.h"


template <typename T>
class Meridian
{
	private:
		T lon;						//Longitude of the meridian
		TInterval <T> lat_interval;			//Latitude interval of the meridian
		T dlat;						//Latitude step;
		T lat_min_shift;				//Shift of the initial meridian point: [lat_min + lat_min_shift, lon], default = 0
		T lat_max_shift;				//Shift of the last meridian point: [lat_max - lat_max_shift, lon], default = 0
		TVector <T> lats;				//Latitude of the meridian points

	public:
		Meridian(): lon(0), lat_interval{ MIN_LAT, MAX_LAT }, dlat(10.0), lat_min_shift(0.0), lat_max_shift(0) { createMeridian(); };
		Meridian(const T lon_) : lon(lon_), lat_interval{ MIN_LAT, MAX_LAT }, dlat(10.0), lat_min_shift(0.0), lat_max_shift(0) { createMeridian(); };
		Meridian(const T lon_, const TInterval<T> lat_interval_, const T dlat_, const T lat_min_shift_, const T lat_max_shift_):
			lon(lon_), lat_interval(lat_interval_), dlat(dlat_), lat_min_shift(lat_min_shift_), lat_max_shift(lat_max_shift_) {createMeridian();};
		
	public:
		T getLon() const { return lon; }
		TInterval <T> & getLatInterval() const { return lat_interval; }
		T getDLat() const { return dlat; }
		T getLatMinShift() const { return lat_min_shift; }
		T getLatMaxShift() const { return lat_max_shift; }
		TVector <T> & getLats() const { return lats; }
		
		template <typename Point>
		void project(const std::shared_ptr <Projection <T> > proj, const T alpha, TVector <Point> &mer);

		void print(std::ostream &output = std::cout);

	private:
		void createMeridian();
};

#include "Meridian.hpp"

#endif