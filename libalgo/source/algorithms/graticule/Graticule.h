// Description: Compute craticule and export to DXF

// Copyright (c) 2010 - 2013
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


#ifndef Graticule_H
#define Graticule_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


//Forward declarations
template <typename T>
class Meridian;

template <typename T>
class Parallel;

template <typename T>
struct TMeridiansListF;

template <typename T>
struct TParallelsListF;

template <typename T>
class Projection;


//Graticule type: computed for normal aspect given by (lat, lon) or for transverse/oblique aspect given by (lat_trans, lon_trans)
typedef enum
{
        NormalGraticule = 0,
        TransformedGraticule,
} TGraticuleType;


//Graticule error: Was error caused by lat or lon?
typedef enum
{
        LatitudeError = 0,
        LongitudeError,
} TGraticuleError;


//List of intervals for meridians and parallels
template <typename T>
struct TIntervals
{
        typedef std::list < TMinMax <T> > Type;
};


//Auxiliary structure prevents the compiler error in VS2010+
template<typename T>
struct nil
{
        typedef T type;
};


//Compute graticule with specific lat, lon increments and steps
class Graticule
{
        public:

                template <typename T>
                static void computeGraticule ( const Projection <T> *proj, const TMinMax<T> lat_interval, const TMinMax<T> lon_interval, const T lat_step, const T lon_step, const T d_lat, const T d_lon, const T alpha, const TGraticuleType grat_type,  typename TMeridiansListF <T> ::Type & meridians, typename TParallelsListF <T> ::Type & parallels, Container <Node3DCartesian <T> *> *points, unsigned int & index );

        private:
                template <typename T>
                static void computeMeridians ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const TMinMax<T> &lon_interval, const T lon_step, const T d_lat, const T alpha, const TGraticuleType grat_type, typename TMeridiansListF <T> ::Type & meridians, Container <Node3DCartesian <T> *> *points, unsigned int & index, T &lat_error, T &lon_error );

                template <typename T>
                static void computeParallels ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const TMinMax<T> &lon_interval, const T lat_step, const T d_lon, const T alpha, const TGraticuleType grat_type, typename TParallelsListF <T> ::Type & parallels, Container <Node3DCartesian <T> *> *points, unsigned int & index, T &lat_error, T &lon_error );

                template <typename T>
                static void computeMeridian ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const T d_lat, const T alpha, const TGraticuleType grat_type, Meridian <T> &meridian, Container <Node3DCartesian <T> *> *points, unsigned int & index );

                template <typename T>
                static void computeParallel ( const Projection <T> *proj, const TMinMax<T> &lon_interval, const T d_lon, const T alpha, const TGraticuleType grat_type, Parallel <T> &parallel, Container <Node3DCartesian <T> *> *points, unsigned int & index );

                template <typename T>
                static void createLatIntervals ( const TMinMax <T> &lat_interval, const T latp, typename TIntervals <T>::Type &lat_intervals );

                template <typename T>
                static void createLonIntervals ( const TMinMax <T> &lon_interval, const T lonp, typename TIntervals <T>::Type &lon_intervals );

                template <typename Intervals>
                static void splitIntervals ( Intervals &intervals, typename Intervals ::iterator &i_intervals, const typename nil <typename Intervals::value_type>::type::value_type &lat_lon_error );
};

#include "Graticule.hpp"

#endif
