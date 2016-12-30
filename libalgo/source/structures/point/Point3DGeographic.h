// Description: 3D geographic point

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


#ifndef Point3DGeographic_H
#define Point3DGeographic_H

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <ostream>
#include <iomanip>


//3D geographic point (lat, lon, elevation)
template <typename T>
class Point3DGeographic
{
        private:
                static unsigned int points_geo_id_counter;	//Static variable: counter of created points

        protected:
                unsigned int point_id;				//Internal point ID (start from 0)
                std::string point_label;			//Point label ( or point ID ) loaded from file ( may contain alphanumeric characters )
                T lat;						//Point latitude
                T lon;						//Point longitude
                T H;						//Point height

        public:
                //Store type of the point coordinates
                typedef T Type;
	
        public:
                //Constructors
		Point3DGeographic() : point_id(points_geo_id_counter++), point_label(""), lat(0), lon(0), H(0) {}
                Point3DGeographic ( const T lat_, const T lon_, const T H_ = 0 ) : point_id ( points_geo_id_counter ++ ), point_label ( "" ), lat ( lat_ ), lon ( lon_ ), H ( H_ ) { }
		Point3DGeographic(const std::string &point_label_, const T lat_, const T lon_, const T H_ = 0) : point_label (point_label_), point_id(points_geo_id_counter++), lat(lat_), lon(lon_), H(H_) {}
		Point3DGeographic(const Point3DGeographic <T> *p) : point_label(p.point_label), point_id(p->point_id), lat(p->lat), lon(p->lon), H(p->H) {}
		Point3DGeographic(const Point3DGeographic <T> &p) : point_label(p.point_label), point_id(p.point_id), lat(p.lat), lon(p.lon), H(p.H) {}
		virtual ~Point3DGeographic() {}

        public:
                //Operators
                Point3DGeographic & operator = ( const Point3DGeographic <T> &p );
                bool operator == ( const Point3DGeographic <T> &p ) const;
                bool operator != ( const Point3DGeographic <T> &p ) const {return ! ( *this == p );}

        public:
                //Other methods
                unsigned short getPointID() const {return point_id;}
                std::string & getPointLabel() const {return point_label; }
                T getLat() const {return lat;}
                T getLon() const {return lon;}
                T getH() const {return H;}

                void updateID () {point_id = points_geo_id_counter ++;}
                void setPointID ( const unsigned int point_id_ )  {point_id = point_id_;}
		void setPointLabel(const std::string &point_label_) { point_label = point_label_; }
                void setLat ( const T lat_ ) {lat = lat_;}
                void setLon ( const T lon_ ) {lon = lon_;}
                void setH ( const T H_ ) {H = H_;}

        public:
                virtual void print ( std::ostream * output = &std::cout ) const;
                virtual Point3DGeographic <T> *clone() const {return new Point3DGeographic <T> ( this->point_label, this->lat, this->lon, this->H ) ;}
};

#include "Point3DGeographic.hpp"

#endif
