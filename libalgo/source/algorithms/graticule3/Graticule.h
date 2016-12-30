// Description: Compute craticule points and apply a projection

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
//3

#ifndef Graticule_H
#define Graticule_H

#include <list>
#include <memory>

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TVector2D.h"
#include "libalgo/source/types/TList.h"
#include "libalgo/source/types/TInterval.h"
#include "libalgo/source/types/TGraticuleError.h"

#include "libalgo/source/structures/graticule/Meridian.h"
#include "libalgo/source/structures/graticule/Parallel.h"


//Compute graticule with specific lat, lon increments and steps
class Graticule
{
	
        public:

                template <typename T>
                static void createGraticule ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> lat_interval, const TInterval<T> lon_interval, const T lat_step, const T lon_step, const T d_lat, const T d_lon, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point3DCartesian <T> > & meridians_proj, TVector <Parallel <T> > &parallels, TVector2D <Point3DCartesian <T> > &parallels_proj );

        private:
		
                template <typename T, typename Point>
                static void createMeridians ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lon_step, const T d_lat, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point> & meridians_proj, T &lat_error, T &lon_error );

                template <typename T, typename Point>
                static void createParallels ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lat_step, const T d_lon, const T alpha, TVector <Parallel <T> > &parallels, TVector2D <Point> & parallels_proj, T &lat_error, T &lon_error );

		template <typename T, typename Point>
		static void createMeridianFragment(const std::shared_ptr <Projection <T> > proj, const T lon, const TInterval <T> lat_interval, const T d_lat, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point> & meridians_proj);

		template <typename T, typename Point>
		static void createParallelFragment(const std::shared_ptr <Projection <T> > proj, const T lat, const TInterval <T> lon_interval, const T d_lon, const T alpha, TVector <Parallel <T> > &parallels, TVector2D <Point> & parallels_proj);

                template <typename T>
		static TVector<TInterval<T> > splitLatInterval(const TInterval <T> &lat_interval, const T lon, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection lon_dir, const T lon0);

		template <typename T>
		static TVector<TInterval<T> > splitLonInterval(const TInterval <T> &lon_interval, const T lat, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection lon_dir, const T lon0);
};

#include "Graticule.hpp"

#endif
