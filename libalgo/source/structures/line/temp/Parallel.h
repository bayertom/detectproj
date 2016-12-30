
// Description: Cartographic parallel given by start / end pointers

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


#ifndef Parallel_H
#define Parallel_H

#include <stdio.h>
#include <set>

#include <ostream>

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/comparators/sortParallelsByLat.h"


//Forward declaration
template <typename T>
class Node3DCartesianProjected;

template <typename T>
struct TFittingLine;

//Parallel given by lat and indices of points
template <typename T>
class Parallel
{
        protected:
                T lat;					//Latitude of the parallel
                TIndexList points_indices;		//Indices of points in list

        public:
                typedef T Type;				//Store type representing parallel lat

        public:
                Parallel () : lat ( 0 ), points_indices ( 0 ) {}
                Parallel ( const T lat_ ) : lat ( lat_ ), points_indices ( 0 ) {}
                Parallel ( const T lat_, const TIndexList & points_indices_ ) :
                        lat ( lat_ ), points_indices ( points_indices_ ) {}

                template <typename Point>
                Parallel ( const TFittingLine <Point> &line );

        public:
                Parallel ( const Parallel <T> &p );
                bool operator < ( const Parallel <T> &p ) const {return this->lat < p->lat;}


        public:
                T getLat() const {return lat;}
                T getCoord() const {return lat;}
                TIndexList const &  getPointsIndices () const {return points_indices;}
                TIndexList &  getPointsIndices () {return points_indices;}
                unsigned int getPointsSize() const {return points_indices.size();}

                void setLat ( const T lat_ ) {lat = lat_;}
                void setPointsIndices ( const TIndexList & points_indices_ ) {points_indices = points_indices_;}

        public:

                void print ( std::ostream * output = &std::cout ) const;

                template <typename Point>
                void print ( const Container <Point> *points, std::ostream * output = &std::cout ) const;
};

#include "Parallel.hpp"

//List of parallels sorted by the lat
template <typename T>
struct TParallelsList
{
        typedef std::set <Parallel <T>, sortParallelsByLat <Parallel <T> > > Type;					//List of parallels
};


//List of parallel fragments sorted by the lat
template <typename T>
struct TParallelsListF
{
        typedef std::multiset <Parallel <T>, sortParallelsByLat <Parallel <T> > > Type;					//List of parallels
};

#endif

