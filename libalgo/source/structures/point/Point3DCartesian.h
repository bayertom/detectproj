// Description: 3D cartesian point

// Copyright (c) 2010 - 2016
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


#ifndef Point3DCartesian_H
#define Point3DCartesian_H

#include <iostream>
#include <ostream>
#include <iomanip>
#include <string>


template <typename T>
class Point3DCartesian
{
        private:
                static unsigned int points_cart_id_counter;	//Static variable: counter of created points

        protected:
                unsigned int point_id;				//Internal point ID (start from 0)
                std::string point_label;			//Point label ( or point ID ) loaded from file ( may contain alphanumeric characters )
                T x;						//Point x coordinate
                T y;						//Point y coordinate
                T z;						//Point z coordinate

        public:
                //Store type of the point coordiantes
                typedef T Type;

        public:
                //Constructors
                Point3DCartesian () : point_id ( points_cart_id_counter ++ ), point_label ( "" ), x ( 0 ),  y ( 0 ), z ( 0 ) {}
                Point3DCartesian ( const T x_, const T y_, const T z_ = 0 ) : point_id ( points_cart_id_counter ++ ), point_label ( "" ), x ( x_ ), y ( y_ ), z ( z_ ) {}
		Point3DCartesian ( const std::string &point_label_, const T x_, const T y_, const T z_ = 0) : point_id(points_cart_id_counter++), point_label (point_label_), x(x_), y(y_), z(z_) {}
		Point3DCartesian(const Point3DCartesian <T> *p) : point_label(p->point_label), point_id(p->point_id), x(p->x), y(p->y), z(p->z) {}
		Point3DCartesian(const Point3DCartesian <T> &p) : point_label(p.point_label), point_id(p.point_id), x(p.x), y(p.y), z(p.z) {}
		virtual ~Point3DCartesian() {};

        public:
                //Operators
                Point3DCartesian & operator = ( const Point3DCartesian <T> &p );
                bool operator == ( const Point3DCartesian <T> &p ) const;
                bool operator != ( const Point3DCartesian <T> &p ) const {return ! ( *this == p );}

        public:
                //Other functions
                unsigned int getPointID() const {return point_id; }
                std::string & getPointLabel() const {return point_label; }
                std::string & getPointLabel() {return point_label; }

                T getX() const {return x;}
                T getY() const {return y;}
                T getZ() const {return z;}
                T getCoordinate ( bool swap ) const {return ( swap ? y : x );}

                void updateID () {point_id = points_cart_id_counter ++;}
		void setPointLabel(const std::string & point_label_) { point_label = point_label_; }
                void setX ( const T x_ ) {x = x_;}
                void setY ( const T y_ ) {y = y_;}
                void setZ ( const T z_ ) {z = z_;}

        public:
                //Virtual functions
                virtual void print ( std::ostream * output = &std::cout ) const;
                virtual Point3DCartesian <T> *clone() const {return new Point3DCartesian <T> ( this->point_label, this->x, this->y, this->z );}
};


template <typename T>
void operator << ( std::ostream & output, const Point3DCartesian <T> &p )
{
        //Print point: overloaded operator, common for all derived class
        p.print ( &output ) ;
}


#include "Point3DCartesian.hpp"

#endif
