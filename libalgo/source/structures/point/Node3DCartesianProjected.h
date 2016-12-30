// Description: 3D Cartesian node storing cartographic parameters and indices to adjacent meridian/parallel points

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


#ifndef Node3DCartesianProjected_H
#define Node3DCartesianProjected_H

#include <stdio.h>
#include <iostream>
#include <ostream>
#include <set>

#include "libalgo/source/algorithms/cartdistortion/CartDistortion.h"

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class Point3DCartesian;

template <typename T>
class Point3DCartesianProjected;

template <typename T>
class HalfEdge;

template <typename T>
class Face;


//Projected point with indices to adjacent points in geographic network (to meridians and parallel points)
template <typename T>
class Node3DCartesianProjected : public Node3DCartesian <T>, public Point3DCartesianProjected <T>
{
        protected:

                int parallel_point_index_prev;		//Index of the previous parallel point
                int parallel_point_index_next;		//Index of the next parallel point
                int meridian_point_index_next;		//Index of the previous meridian point
                int meridian_point_index_prev;		//Index of the next meridian point

        public:
                Node3DCartesianProjected() : Point3DCartesian <T> (), Node3DCartesian <T>(), Point3DCartesianProjected <T>(),
                        parallel_point_index_prev ( -1 ), parallel_point_index_next ( -1 ), meridian_point_index_next ( -1 ), meridian_point_index_prev ( -1 ) {}

                Node3DCartesianProjected ( const T x_, const T y_, const T z_ = 0 ) : Point3DCartesian <T> ( x_, y_, z_ ), Node3DCartesian <T> ( x_, y_, z_ ), Point3DCartesianProjected <T> ( x_, y_, z_ ),
                        parallel_point_index_prev ( -1 ), parallel_point_index_next ( -1 ), meridian_point_index_next ( -1 ), meridian_point_index_prev ( -1 ) {}

                Node3DCartesianProjected ( const T x_, const T y_, const T z_, const T h_, const T k_, const T s_, const T theta_, const TTissotIndicatrix <T> tiss_, const T w_ ) : Point3DCartesian <T> ( x_, y_, z_ ), Node3DCartesian <T> ( x_, y_, z_ ),
                        Point3DCartesianProjected <T> ( x_, y_, z_, h_, k_, s_, theta_, tiss_, w_ ), parallel_point_index_prev ( -1 ), parallel_point_index_next ( -1 ),
                        meridian_point_index_next ( -1 ), meridian_point_index_prev ( -1 ) {}

                Node3DCartesianProjected ( const T x_, const T y_, const T z_, const HalfEdge <T> *half_edge, const Face <T> *face ) : Point3DCartesian <T> ( x_, y_, z_ ), Node3DCartesian <T> ( x_, y_, z_, half_edge, face ),
                        Point3DCartesianProjected <T> ( x_, y_ , z_ ), parallel_point_index_prev ( -1 ), parallel_point_index_next ( -1 ), meridian_point_index_next ( -1 ), meridian_point_index_prev ( -1 ) {}

                Node3DCartesianProjected ( const T x_, const T y_, const T z_, const int parallel_point_index_prev_, const int parallel_point_index_next_, const int meridian_point_index_next_, const int meridian_point_index_prev_ ) :
                        Point3DCartesian <T> ( x_, y_, z_ ), Node3DCartesian <T> ( x_, y_, z_ ), Point3DCartesianProjected <T> ( x_, y_, z_ ), parallel_point_index_prev ( parallel_point_index_prev_ ), parallel_point_index_next ( parallel_point_index_next_ ),
                        meridian_point_index_next ( meridian_point_index_next_ ), meridian_point_index_prev ( meridian_point_index_prev_ ) {}

                Node3DCartesianProjected ( const std::string & point_label_, const T x_, const T y_, const T z_, const T h_, const T k_, const T s_, const T theta_, const TTissotIndicatrix <T> &tiss_, const T w_, HalfEdge <T> *half_edge, Face <T> *face,
                                           const int parallel_point_index_prev_, const int parallel_point_index_next_, const int meridian_point_index_prev_, const int meridian_point_index_next_ ) : Point3DCartesian <T> ( point_label_, x_, y_, z_ ), Node3DCartesian <T> ( point_label_, x_, y_, z_, half_edge, face ),
                        Point3DCartesianProjected <T> ( point_label_, x_, y_, z_, h_, k_, s_, theta_, tiss_, w_ ), parallel_point_index_prev ( parallel_point_index_prev_ ), parallel_point_index_next ( parallel_point_index_next_ ), meridian_point_index_next ( meridian_point_index_next_ ),
                        meridian_point_index_prev ( meridian_point_index_prev_ ) {}

                Node3DCartesianProjected ( const Node3DCartesianProjected <T> *n );
                virtual ~Node3DCartesianProjected();

        public:
                //Operators
                bool operator == ( const Node3DCartesianProjected <T> &n ) const;
                bool operator != ( const Node3DCartesianProjected <T> &n ) const {return ! ( *this == n );};

        public:
                //Other functions
                int getParallelPointIndexPrev() const {return parallel_point_index_prev;}
                int getParallelPointIndexNext() const {return parallel_point_index_next;}
                int getMeridianPointIndexNext() const {return meridian_point_index_next;}
                int getMeridianPointIndexPrev() const {return meridian_point_index_prev;}

                void setParallelPointIndexPrev ( const int parallel_point_index_prev_ ) {parallel_point_index_prev = parallel_point_index_prev_;}
                void setParallelPointIndexNext ( const int parallel_point_index_next_ ) {parallel_point_index_next = parallel_point_index_next_;}
                void setMeridianPointIndexNext ( const int meridian_point_index_next_ ) {meridian_point_index_next = meridian_point_index_next_;}
                void setMeridianPointIndexPrev ( const int meridian_point_index_prev_ ) {meridian_point_index_prev = meridian_point_index_prev_;}

        public:
                virtual void print ( std::ostream * output = &std::cout ) const;

                virtual Node3DCartesianProjected <T> *clone()
                {
                        return new Node3DCartesianProjected <T> ( this->point_label, this->x, this->y, this->z,
                                        this->h, this->k, this->s, this->theta, this->tiss, this->w, this->half_edge, this->face, this->parallel_point_index_prev,
                                        this->parallel_point_index_prev, this->meridian_point_index_prev, this->meridian_point_index_next );
                }
};

#include "Node3DCartesianProjected.hpp"

#endif
