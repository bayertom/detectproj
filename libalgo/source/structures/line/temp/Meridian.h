// Description: Cartographic meridian defined by start / end point

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


#ifndef Meridian_H
#define Meridian_H

#include <stdio.h>
#include <set>
#include <ostream>

#include "libalgo/source/structures/list/Container.h"
#include "libalgo/source/comparators/sortMeridiansByLon.h"
#include "libalgo/source/structures/list/IndexLists.h"

//Forward declaration
template <typename T>
class Node3DCartesianProjected;

template <typename T>
struct TFittingLine;

//Meridian given by lon and indices of points
template <typename T>
class Meridian
{
        private:
                T lon;					//Longitude of the meridian
                TIndexList points_indices;		//Indices of points in list

        public:
                typedef T Type;				//Store type representing meridian lon

        public:
                Meridian () : lon ( 0 ), points_indices ( 0 ) {}
                Meridian ( const T lon_ ) : lon ( lon_ ), points_indices ( 0 ) {}
                Meridian ( const T lon_, const TIndexList & meridian_points_indices_ ) :
                        lon ( lon_ ), points_indices ( meridian_points_indices_ ) {}

                template <typename Point>
                Meridian ( const TFittingLine <Point> &line );

        public:
                Meridian ( const Meridian <T> *m );
                bool operator < ( const Meridian <T> &m ) const {return this->lon < m.lon;}

        public:
                //Other methods
                T getLon() const {return lon;}
                T getCoord() const {return lon;}
                TIndexList const & getPointsIndices () const {return points_indices;}
                TIndexList &  getPointsIndices () {return points_indices;}
                unsigned int getPointsSize() const {return points_indices.size();}

                void setLon ( const T lon_ ) {lon = lon_;}
                void setPointsIndices ( const TIndexList & meridian_points_indices_ ) {points_indices = meridian_points_indices_;}

        public:

                void print ( std::ostream * output = &std::cout ) const;

                template <typename Point>
                void print ( const Container <Point> *points, std::ostream * output = &std::cout ) const;
};

#include "Meridian.hpp"


//List of meridians sorted by the lon
template <typename T>
struct TMeridiansList
{
        typedef std::set <Meridian <T>, sortMeridiansByLon <Meridian <T> > > Type;					//List of meridians
};


//List of meridian fragmensts sorted by the lon
template <typename T>
struct TMeridiansListF
{
        typedef std::multiset <Meridian <T>, sortMeridiansByLon <Meridian <T> > > Type;					//List of meridians
};


#endif

