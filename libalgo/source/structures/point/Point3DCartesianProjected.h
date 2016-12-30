// Description: 3D cartesian cartographic point (cartographic distortions are stored)

// Copyright (c) 2010 - 2016
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can reribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is ributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


#ifndef Point3DCartesianProjected_H
#define Point3DCartesianProjected_H

#include <iostream>
#include <ostream>

#include "Point3DCartesian.h"

#include "libalgo/source/algorithms/cartdistortion/CartDistortion.h"


template <typename T>
class Point3DCartesianProjected : virtual public Point3DCartesian <T>
{
        protected:

                T h;				//Meridian scale
                T k;				//Parallel scale
                T s;				//Areal scale
                T theta;			//Angle between meridian and parallel
                TTissotIndicatrix <T> tiss;	//Tissot Indicatrix parameters
                T w;				//Max angular ortion

        public:
                Point3DCartesianProjected() : Point3DCartesian <T> (),
                        h ( 0 ), k ( 0 ), s ( 0 ), theta ( 0 ), tiss ( ), w ( 0 ) {}
                Point3DCartesianProjected ( const T x_, const T y_, const T z_ = 0 ) : Point3DCartesian <T> ( x_, y_, z_ ),
                        h ( 0 ), k ( 0 ), s ( 0 ), theta ( 0 ), tiss ( ), w ( 0 ) {}
                Point3DCartesianProjected ( const T x_, const T y_, const T z_, const T h_, const T k_, const T s_, const T theta_, const TTissotIndicatrix <T> &tiss_, const T w_ ) : Point3DCartesian <T> ( x_, y_, z_ ),
                        h ( h_ ), k ( k_ ), s ( s_ ), theta ( theta_ ), tiss ( tiss_ ), w ( w_ ) {}
                Point3DCartesianProjected ( const std::string &point_label, const T x_, const T y_, const T z_, const T h_, const T k_, const T s_, const T theta_, const TTissotIndicatrix <T> &tiss_, const T w_ ) : Point3DCartesian <T> ( point_label, x_, y_, z_ ),
                        h ( h_ ), k ( k_ ), s ( s_ ), theta ( theta_ ), tiss ( tiss_ ), w ( w_ ) {}

                virtual ~Point3DCartesianProjected() {}

        public:
                //Operators
                bool operator == ( const Point3DCartesianProjected <T> &p ) const;
                bool operator != ( const Point3DCartesianProjected <T> &p ) const {return ! ( *this == p );}

        public:

                //Other functions

                T getH() const {return h;}
                T getK() const {return k;}
                T getS() const {return s;}
                T getTheta() const {return theta;}
                TTissotIndicatrix <T> const & getTiss() const {return tiss;}
                T const getW() const {return w;}

        public:

                void setH ( const T h_ ) {h = h_;}
                void setK ( const T k_ ) {k = k_;}
                void setS ( const T s_ ) {s = s_;}
                void setTheta ( const T theta_ ) {theta = theta_;}
                void setTiss ( const TTissotIndicatrix <T> & tiss_ ) {tiss = tiss_;}
                void setW ( const T w_ ) {w = w_;}

        public:
                virtual void print ( std::ostream * output = &std::cout ) const;

                virtual Point3DCartesianProjected <T> *clone() const
                {
                        return new Point3DCartesianProjected <T> ( this->point_label, this->x, this->y, this->z,
                                        this->h, this->k, this->s, this->theta, this->tiss, this->w ) ;
                }
};

#include "Point3DCartesianProjected.hpp"

#endif
