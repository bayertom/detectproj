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


#ifndef Graticule_HPP
#define Graticule_HPP

#include <math.h>

#include "libalgo/source/exceptions/MathException.h"

template <typename T>
void Graticule::computeGraticule ( const Projection <T> *proj, const TMinMax<T> lat_interval, const TMinMax<T> lon_interval, const T lat_step, const T lon_step, const T d_lat, const T d_lon, const T alpha, const TGraticuleType grat_type, typename TMeridiansListF<T> ::Type & meridians, typename TParallelsListF <T> ::Type & parallels, Container <Node3DCartesian <T> *> *points, unsigned int & index )
{
        //Compute graticule given by meridians and parallels
        TMinMax <T> lat_interval_part, lon_interval_part;

        //Cartographic pole
        const T latp = proj->getCartPole().getLat(), lonp = proj->getCartPole().getLon();
	const T lon0 = proj->getLon0();

        //Get size of the interval
        const T lat_interval_width = lat_interval.max_val - lat_interval.min_val;
        const T lon_interval_width = lon_interval.max_val - lon_interval.min_val;

        //Generate meridians and parallels
        if ( lat_step <= lat_interval_width || lon_step <= lon_interval_width )
        {
                //Create lat and lon intervals to avoid basic singular points
                typename TIntervals <T>::Type lat_intervals, lon_intervals;
                createLatIntervals ( lat_interval, latp, lat_intervals );
		
		//Normal aspect
		if ( (fabs(latp - MAX_LAT) < MIN_FLOAT) && (fabs(lonp) < MIN_FLOAT) )
			createLonIntervals ( lon_interval, lon0, lon_intervals );

		//Transverse or oblique aspects
		else
			createLonIntervals(lon_interval, lonp, lon_intervals);

                //Process all meridians and parallels: automatic detection of additional singular points
                unsigned int split_amount = 0;

                for ( typename TIntervals <T>::Type::iterator i_lat_intervals = lat_intervals.begin(); i_lat_intervals != lat_intervals.end(); i_lat_intervals ++ )
                {
                        for ( typename TIntervals <T>::Type::iterator i_lon_intervals = lon_intervals.begin() ; i_lon_intervals != lon_intervals.end(); )
                        {
                                //Create temporary meridians and parallels
                                typename TMeridiansListF <T> ::Type meridians_temp;
                                typename TParallelsListF <T> ::Type parallels_temp;

                                //Compute meridians and parallels
                                T lat_error = 0.0, lon_error = 0.0;

                                try
                                {
                                        //Compute meridians and parallels
                                        computeMeridians ( proj, *i_lat_intervals, *i_lon_intervals, lon_step, d_lat, alpha, grat_type, meridians_temp, points, index, lat_error, lon_error );
                                        computeParallels ( proj, *i_lat_intervals, *i_lon_intervals, lat_step, d_lon, alpha, grat_type, parallels_temp, points, index, lat_error, lon_error );

                                        //Copy temporary meridians and parallels to the output data structure
                                        std::copy ( meridians_temp.begin(), meridians_temp.end(), std::inserter ( meridians, meridians.begin() ) );
                                        std::copy ( parallels_temp.begin(), parallels_temp.end(), std::inserter ( parallels, parallels.begin() ) );

                                        //Increment lon intervals only if comuptation was successful
                                        i_lon_intervals ++;
                                }

                                //Exception
                                catch ( MathException <T> &error )
                                {
                                        //Too many splits, projection is suspected, stop computations
                                        if ( split_amount > 100 )
                                        {
                                                meridians.clear(); parallels.clear();
                                                return;
                                        }

					//Empty lat interval, delete
					if (fabs((*i_lat_intervals).max_val - (*i_lat_intervals).min_val) < 10 * GRATICULE_ANGLE_SHIFT)
					{
						i_lat_intervals = lat_intervals.erase(i_lat_intervals);
						++i_lat_intervals;
					}

					//Lat value is lower bound: shift lower bound
					else if ((fabs((*i_lat_intervals).min_val - lat_error) <  2 * GRATICULE_ANGLE_SHIFT) && (fabs((*i_lat_intervals).max_val - (*i_lat_intervals).min_val) > 10 * GRATICULE_ANGLE_SHIFT))
					{
						(*i_lat_intervals).min_val += 2 *GRATICULE_ANGLE_SHIFT;
					}

					//Lat value is upper bound: shift upper bound
					else if ((fabs((*i_lat_intervals).max_val - lat_error) <  2 * GRATICULE_ANGLE_SHIFT) && (fabs((*i_lat_intervals).max_val - (*i_lat_intervals).min_val) > 10 * GRATICULE_ANGLE_SHIFT))
					{
						(*i_lat_intervals).max_val -= 2 * GRATICULE_ANGLE_SHIFT;
					}

					//Lat value inside interval: split intervals
					else if (((*i_lat_intervals).min_val < lat_error) && ((*i_lat_intervals).max_val > lat_error))
					{
						splitIntervals(lat_intervals, i_lat_intervals, lat_error);

						//Increment split amount
						split_amount++;
					}

					//Empty lon interval, delete
					if (fabs((*i_lon_intervals).max_val - (*i_lon_intervals).min_val) < 10 * GRATICULE_ANGLE_SHIFT)
					{
						i_lon_intervals = lon_intervals.erase(i_lon_intervals);
						++i_lon_intervals;
					}

					//Lon value is lower bound : shift lower bound
					else if ((fabs((*i_lon_intervals).min_val - lon_error) <  2 *GRATICULE_ANGLE_SHIFT) && (fabs((*i_lon_intervals).max_val - (*i_lon_intervals).min_val) > 10 * GRATICULE_ANGLE_SHIFT))
					{
						(*i_lon_intervals).min_val += 2 * GRATICULE_ANGLE_SHIFT;
					}
					
					//Lon value is upper bound: shift upper bound
					else if ((fabs((*i_lon_intervals).max_val - lon_error) <  2 * GRATICULE_ANGLE_SHIFT) && (fabs((*i_lon_intervals).max_val - (*i_lon_intervals).min_val) > 10 * GRATICULE_ANGLE_SHIFT))
					{
						(*i_lon_intervals).max_val -= 2 * GRATICULE_ANGLE_SHIFT;
					}

                                        //Lon value inside interval: split intervals
                                        else if ( ( ( *i_lon_intervals ).min_val < lon_error ) && ( ( *i_lon_intervals ).max_val > lon_error ) )
                                        {
                                                splitIntervals ( lon_intervals, i_lon_intervals, lon_error );

						//Increment split amount
						split_amount++;
                                        }

					//std::cout << (*i_lat_intervals).min_val << "   " << (*i_lat_intervals).max_val << '\n';
					//std::cout << (*i_lon_intervals).min_val << "   " << (*i_lon_intervals).max_val << '\n';
                                }

                                //Other than math error: error in equation or in parser
                                catch ( Exception &error )
                                {
                                        //Clear meridians and parallels
                                        meridians.clear(); parallels.clear();
                                        return;
                                }
                        }
                }
        }
}


template <typename T>
void Graticule:: createLatIntervals ( const TMinMax <T> &lat_interval, const T latp, typename TIntervals <T>::Type &lat_intervals )
{
        //Split given lat interval into several sub-intervals
        TMinMax <T> lat_interval_part;

        //First lat interval ( -90, latp )
        lat_interval_part.min_val = std::max ( MIN_LAT, lat_interval.min_val ); lat_interval_part.max_val = std::min ( latp, lat_interval.max_val );

        //Test interval and add to the list
        if ( lat_interval_part.max_val > lat_interval_part.min_val )
        {
                lat_intervals.push_back ( lat_interval_part );
        }

        //Second lat interval ( latp, 90 )
        lat_interval_part.min_val = std::max ( latp, lat_interval.min_val ); lat_interval_part.max_val = std::min ( lat_interval.max_val, MAX_LAT );

        //Test interval and add to the list
        if ( lat_interval_part.max_val > lat_interval_part.min_val )
        {
                lat_intervals.push_back ( lat_interval_part );
        }
}


template <typename T>
void Graticule::createLonIntervals ( const TMinMax <T> &lon_interval, const T lonp, typename TIntervals <T>::Type &lon_intervals )
{
        //Split given lon interval into several sub-intervals
        TMinMax <T> lon_interval_part;

        //Split given interval into sub intervals: lonp >=0
        if ( lonp >= 0 )
        {
                //First interval <-180, lonp - 180)
                lon_interval_part.min_val = std::max ( MIN_LON, lon_interval.min_val ); lon_interval_part.max_val =  std::min ( lonp - 180.0, lon_interval.max_val );

                //Test first interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }

                //Second interval (lonp - 180, lonp>
                lon_interval_part.min_val = std::max ( lonp - 180.0, lon_interval.min_val ); lon_interval_part.max_val = std::min ( lonp, lon_interval.max_val );

                //Test second interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }

                //Third interval (lonp, MAX_LON>
                lon_interval_part.min_val = std::max ( lonp, lon_interval.min_val ); lon_interval_part.max_val = std::min ( MAX_LON, lon_interval.max_val );

                //Test third interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }
        }

        //Split given interval into sub intervals: lonp < 0
        else
        {
                //First interval <-180, lonp >
                lon_interval_part.min_val = std::max ( MIN_LON, lon_interval.min_val ); lon_interval_part.max_val = std::min ( lonp, lon_interval.max_val );

                //Test first interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }

                //Second interval <lonp, lonp + 180>
                lon_interval_part.min_val = std::max ( lonp, lon_interval.min_val ); lon_interval_part.max_val = std::min ( lonp + 180.0, lon_interval.max_val );

                //Test second interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }

                //Third interval (lonp + 180, 180>
                lon_interval_part.min_val = std::max ( lonp + 180.0, lon_interval.min_val ); lon_interval_part.max_val = std::min ( MAX_LON, lon_interval.max_val );

                //Test third interval and add to the list
                if ( lon_interval_part.max_val > lon_interval_part.min_val )
                {
                        lon_intervals.push_back ( lon_interval_part );
                }
        }
}


template <typename Intervals>
void Graticule::splitIntervals ( Intervals &intervals, typename Intervals::iterator &i_intervals, const typename nil <typename Intervals::value_type>::type::value_type &lat_lon_error )
{
        //Split one interval in 2
        TMinMax <typename Intervals::value_type::value_type> interval_temp;
        const typename Intervals::value_type::value_type lat_lon_max = ( *i_intervals ).max_val;

        //Split old interval
        ( *i_intervals ).max_val = lat_lon_error;

        //Create new interval
        interval_temp.min_val = lat_lon_error;
        interval_temp.max_val = lat_lon_max;

        //Add interval to the list
        intervals.insert ( i_intervals, interval_temp );

        //Move iterator to the previous item
        i_intervals--;
}


template <typename T>
void Graticule::computeMeridians ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const TMinMax<T> &lon_interval, const T lon_step, const T d_lat, const T alpha, const TGraticuleType grat_type, typename TMeridiansListF<T> ::Type & meridians, Container <Node3DCartesian <T> *> *points, unsigned int & index, T &lat_error, T &lon_error )
{
        //Compute all meridians inside <lon_min, lon_max> interval
        T lon = lon_interval.min_val;

        //Correct longitude: set start value
        if ( fmod ( lon_interval.min_val, lon_step ) != 0.0 )
        {
                if ( lon < 0.0 )
                        lon = lon_interval.min_val -  fmod ( lon_interval.min_val, lon_step );
                else
                        lon = lon_interval.min_val -  fmod ( lon_interval.min_val, lon_step ) + lon_step ;
        }

        //Compute all meridians
        for ( unsigned int i = 0; ( lon  < lon_interval.max_val + lon_step ) && ( lat_interval.min_val < lat_interval.max_val ) && ( lon_interval.min_val <= lon_interval.max_val ); lon += lon_step, i++ )
        {
                //Create new meridian
                Meridian <T> m;

                //Moved meridian lon = lon_min
                if ( fabs ( lon - lon_interval.min_val ) < ANGLE_ROUND_ERROR )
                {
                        //Corect longitude
                        m.setLon ( lon_interval.min_val + GRATICULE_ANGLE_SHIFT );
                }

                //Create new meridian to fill the gap, point closer to lon_min than a lon_step
                else if ( ( ( fabs ( lon - lon_interval.min_val ) > ANGLE_ROUND_ERROR ) && ( fabs ( lon - lon_interval.min_val ) < lon_step ) ) && ( i == 0 ) )
                {
                        //Corect longitude
                        m.setLon ( lon_interval.min_val  + GRATICULE_ANGLE_SHIFT );

                        //Decrement lon_step, compute again with same lon
                        lon -= lon_step;
                }

                //Moved meridian lon = lon_max (last meridian), otherwise do not move
                else if ( fabs ( lon - lon_interval.max_val ) < ANGLE_ROUND_ERROR )
                {
                        //Corect longitude
                        m.setLon ( lon_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Crete new meridian to fill the gap, point closer to lon_max than a lon_step
                else if ( ( lon > lon_interval.max_val ) && ( fabs ( lon - lon_interval.max_val ) < lon_step ) )
                {
                        //Create new meridian and move it before the interval break
                        m.setLon ( lon_interval.max_val  - GRATICULE_ANGLE_SHIFT );
                }

                //Common meridian: not moved
                else
                {
                        m.setLon ( lon );
                }

                try
                {
                        //Compute meridian
                        computeMeridian ( proj, lat_interval, d_lat, alpha, grat_type, m, points, index );

                        //Does it contain more than 1 point?
                        if ( m.getPointsSize() > 1 )
                                //Add meridian to the list
                                meridians.insert ( m );
                }

                //Throw math exception: get error values
                catch ( MathException <T> &error )
                {
                        //Get posible error values
                        lat_error = error.getArg();
                        lon_error = m.getLon();

                        //Throw exception
                        throw;
                }

                //Throw other exception: error in parser or incorrect equation
                catch ( Exception &error )
                {
                        throw;
                }
        }
}


template <typename T>
void Graticule::computeParallels ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const TMinMax<T> &lon_interval, const T lat_step, const T d_lon, const T alpha, const TGraticuleType grat_type, typename TParallelsListF <T> ::Type & parallels, Container <Node3DCartesian <T> *> *points, unsigned int & index, T &lat_error, T &lon_error )
{
        //Compute all parallels
        T lat = lat_interval.min_val;

        //Correct latitude: set start value
        if ( fmod ( lat_interval.min_val, lat_step ) != 0.0 )
        {
                if ( lat < 0.0 )
                        lat = lat_interval.min_val -  fmod ( lat_interval.min_val, lat_step );
                else
                        lat = lat_interval.min_val -  fmod ( lat_interval.min_val, lat_step ) + lat_step;
        }

        //Compute all parallels given by lat interval
        for ( unsigned int i = 0; ( lat  < lat_interval.max_val + lat_step ) && ( lon_interval.min_val < lon_interval.max_val ) && ( lat_interval.min_val <= lat_interval.max_val );  lat += lat_step, i++ )
        {
                //Create new parallel
                Parallel <T> p ( lat );

                //Moved parallel lat = lat_min:
                if ( fabs ( lat - lat_interval.min_val ) < ANGLE_ROUND_ERROR )
                {
                        //Corect longitude
                        p.setLat ( lat_interval.min_val + GRATICULE_ANGLE_SHIFT );
                }

                //Crete new parallel to fill the gap, point closer to lat_min than a lat_step
                else if ( ( ( fabs ( lat - lat_interval.min_val ) > ANGLE_ROUND_ERROR ) && ( fabs ( lat - lat_interval.min_val ) < lat_step ) ) && ( i == 0 ) )
                {
                        //Correct latitude
                        p.setLat ( lat_interval.min_val + GRATICULE_ANGLE_SHIFT );

                        //Decrement lat, compute again with the same value
                        lat -= lat_step;
                }

                //Move parallel lat = lat_max
                else if ( fabs ( lat - lat_interval.max_val ) < ANGLE_ROUND_ERROR )
                {
                        //Corect latitude
                        p.setLat ( lat_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Crete new parallel to fill the gap, point closer to lat_min than a lat_step
                else if ( ( lat > lat_interval.max_val ) && ( fabs ( lat - lat_interval.max_val ) < lat_step ) )
                {
                        //Corect latitude
                        p.setLat ( lat_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                try
                {
                        //Compute parallel
                        computeParallel ( proj, lon_interval, d_lon, alpha, grat_type, p, points, index );

                        //Does it contain more than 1 point?
                        if ( p.getPointsSize() > 1 )

                                //Add parallel to the list
                                parallels.insert ( p );
                }

                //Throw math exception: get error values
                catch ( MathException <T> &error )
                {
                        //Get posible error values
                        lat_error = p.getLat();
                        lon_error = error.getArg();

                        //Throw exception
                        throw;
                }

                //Throw other exception: error in parser or incorrect equation
                catch ( Exception &error )
                {
                        throw;
                }
        }
}


template <typename T>
void Graticule::computeMeridian ( const Projection <T> *proj, const TMinMax<T> &lat_interval, const T d_lat, const T alpha, const TGraticuleType grat_type, Meridian <T> &meridian, Container <Node3DCartesian <T> *> *points, unsigned int & index )
{
        //Compute one meridian
        TIndexList mer_point_ind;

        //Initialize latitude
        T lat = lat_interval.min_val;

        //Correct latitude: set start value
        if ( fmod ( lat_interval.min_val, d_lat ) != 0.0 )
        {
                if ( lat < 0.0 )
                        lat = lat_interval.min_val -  fmod ( lat_interval.min_val, d_lat );
                else
                        lat = lat_interval.min_val -  fmod ( lat_interval.min_val, d_lat ) + d_lat ;
        }

        //Create meridian points
        for ( unsigned int i = 0; ( lat  < lat_interval.max_val + d_lat ) && ( lat_interval.min_val < lat_interval.max_val ); lat += d_lat, i++ )
        {
                //Create new geographic point
                Point3DGeographic <T> point_geo_temp ( lat, meridian.getLon () );

                //Point closer to lon_min than ANGLE_ROUND_ERROR: fabs(lon-lon_min) < ANGLE_ROUND_ERROR,  move to the  east to avoid the singularity
                if ( fabs ( lat - lat_interval.min_val ) < ANGLE_ROUND_ERROR )
                {
                        point_geo_temp.setLat ( lat_interval.min_val + GRATICULE_ANGLE_SHIFT );
                }

                //Point close to lon_max than ANGLE_ROUND_ERROR: fabs(lon-lon_max) < ANGLE_ROUND_ERROR, move to the west to avoid the singularity
                else if ( fabs ( lat_interval.max_val - lat ) < ANGLE_ROUND_ERROR )
                {
                        point_geo_temp.setLat ( lat_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Point closer to lat_min than a dlat: fabs(lat-lat_min) < d_lat, create a new point to fill the gap
                else if ( ( fabs ( lat - lat_interval.min_val ) > ANGLE_ROUND_ERROR ) && ( fabs ( lat - lat_interval.min_val ) < d_lat ) && ( i == 0 ) )
                {
                        point_geo_temp.setLat ( lat_interval.min_val + GRATICULE_ANGLE_SHIFT );

                        //Decrement lat, compute again with the same lat
                        lat -= d_lat;
                }

                //Point closer to lat_max than d_lat step and lat > lat_max: create new point to fill the gap
                else if ( lat >  lat_interval.max_val )
                {
                        point_geo_temp.setLat ( lat_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Reduce longitude
                T lon_red = ( proj->getLon0() != 0 ? CartTransformation::redLon0 ( meridian.getLon (), proj->getLon0() ) :  meridian.getLon () );

		//Graticule in transverse/oblique aspect: convert the point
		if (grat_type == TransformedGraticule)
		{
			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Convert geographic point to oblique position
			const T lat_trans = CartTransformation::latToLatTrans(point_geo_temp.getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon());
			const T lon_trans = CartTransformation::lonToLonTrans(point_geo_temp.getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon(), trans_lon_dir);

			//Change coordinates of the point
			point_geo_temp.setLat(lat_trans);
			point_geo_temp.setLon(lon_trans);
		}

		//Compute coordinates x, y
		T x = 0, y = 0;

		try
		{
			//Compute equations
			T x_temp = CartTransformation::latLonToX(&point_geo_temp, proj, false);
			T y_temp = CartTransformation::latLonToY(&point_geo_temp, proj, false);

			//Get shifts: need to be additionally subtracted
			const T dx = proj->getDx();
			const T dy = proj->getDy();

			//Non-rotated projection
			if (alpha == 0.0)
			{
				x = x_temp;
				y = y_temp;
			}

			//Rotated projection
			else
			{
				x = (x_temp - dx)* cos(alpha * M_PI / 180) - (y_temp - dy) * sin(alpha * M_PI / 180) + dx;
				y = (x_temp - dx)* sin(alpha * M_PI / 180) + (y_temp - dy)* cos(alpha * M_PI / 180) + dy;
			}
		}

		//Get error argument
		catch (MathException <T> &error)
		{
			//Throw new math error
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: error in coordinate functions (lat/lon).", "Can not compute meridian points.", point_geo_temp.getLat());
		}

                //Create new cartographic point
                Node3DCartesianProjected <T> *point_proj_temp  = new Node3DCartesianProjected <T> ( x, y );

                //Add projected point to the list of points
                points->push_back ( point_proj_temp );

                //Add index of the point to the meridian
                mer_point_ind.push_back ( index ++ );
        }

        //Set all meridian points
        meridian.setPointsIndices ( mer_point_ind );
}


template <typename T>
void Graticule::computeParallel ( const Projection <T> *proj, const TMinMax<T> &lon_interval, const T d_lon, const T alpha, const TGraticuleType grat_type, Parallel <T> &parallel, Container <Node3DCartesian <T> *> *points, unsigned int & index )
{
        //Compute one parallel
        TIndexList par_point_ind;

        //Set initial value of longitude
        T lon = lon_interval.min_val;

        //Correct longitude: set start value
        if ( fmod ( lon_interval.min_val, d_lon ) != 0.0 )
        {
                if ( lon < 0.0 )
                        lon = lon_interval.min_val -  fmod ( lon_interval.min_val, d_lon );
                else
                        lon = lon_interval.min_val -  fmod ( lon_interval.min_val, d_lon ) + d_lon ;
        }

        //Create parallel points
        for ( unsigned int i = 0; ( lon  < lon_interval.max_val + d_lon ) && ( lon_interval.min_val < lon_interval.max_val ); lon += d_lon, i++ )
        {
                //Create new geographic point
                Point3DGeographic <T> point_geo_temp ( parallel.getLat(), lon );

                //Point closer to lon_min than ANGLE_ROUND_ERROR: fabs(lon-lon_min) < ANGLE_ROUND_ERROR,  move to the  east to avoid the singularity
                if ( fabs ( lon - lon_interval.min_val ) < ANGLE_ROUND_ERROR )
                {
                        point_geo_temp.setLon ( lon_interval.min_val + GRATICULE_ANGLE_SHIFT );
                }

                //Point close to lon_max than ANGLE_ROUND_ERROR: fabs(lon-lon_max) < ANGLE_ROUND_ERROR, move to the west to avoid the singularity
                else if ( fabs ( lon_interval.max_val - lon ) < ANGLE_ROUND_ERROR )
                {
                        point_geo_temp.setLon ( lon_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Point closer to lon_min than d_lon step: create new point
                else if ( ( ( fabs ( lon - lon_interval.min_val ) > ANGLE_ROUND_ERROR ) ) && ( fabs ( lon - lon_interval.min_val ) < d_lon ) && ( i == 0 ) )
                {
                        point_geo_temp.setLon ( lon_interval.min_val + GRATICULE_ANGLE_SHIFT );

                        //Decrement lon, compute again with the same lon
                        lon -= d_lon;
                }

                //Point closer to lon_max than d_lon step and lon > lon_max: create a new point to fill the gap
                else if ( lon > lon_interval.max_val )
                {
                        point_geo_temp.setLon ( lon_interval.max_val - GRATICULE_ANGLE_SHIFT );
                }

                //Reduce longitude
		const T lon_point = point_geo_temp.getLon();
		T lon_red = (proj->getLon0() != 0 ? CartTransformation::redLon0(lon_point, proj->getLon0()) : lon_point);
		
		//Graticule in transverse/oblique aspect: convert the point
		if (grat_type == TransformedGraticule)
		{
			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Convert geographic point to oblique position
			const T lat_trans = CartTransformation::latToLatTrans(parallel.getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon());
			const T lon_trans = CartTransformation::lonToLonTrans(parallel.getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon(), trans_lon_dir);

			//Change coordinates of the temporary point
			point_geo_temp.setLat(lat_trans);
			point_geo_temp.setLon(lon_trans);
		}

		//Compute coordinates x, y
		T x = 0, y = 0;

		try
		{
			//Compute equations
			T x_temp = CartTransformation::latLonToX(&point_geo_temp, proj, false);
			T y_temp = CartTransformation::latLonToY(&point_geo_temp, proj, false);

			//Get shifts: need to be additionally subtracted
			const T dx = proj->getDx();
			const T dy = proj->getDy();

			//Non-rotated projection
			if (alpha == 0.0)
			{
				x = x_temp;
				y = y_temp;
			}

			//Rotated projection
			else
			{
				x = (x_temp - dx)* cos(alpha * M_PI / 180) - (y_temp - dy) * sin(alpha * M_PI / 180) + dx;
				y = (x_temp - dx)* sin(alpha * M_PI / 180) + (y_temp - dy) * cos(alpha * M_PI / 180) + dy;
			}
		}

		//Get error argument
		catch (MathException <T> &error)
		{
			//Throw new math error
			throw MathInvalidArgumentException <T>("MathInvalidArgumentException: error in coordinate functions (lat/lon).", "Can not compute parallel points.", lon_point);
		}


                //Create new cartographic point
                Node3DCartesianProjected <T> *point_proj_temp  = new Node3DCartesianProjected <T> ( x, y );

                //Add projected point to the list of points
                points->push_back ( point_proj_temp );

                //Add index of the point to the parallel
                par_point_ind.push_back ( index ++ );
        }

        //Set all parallel points
        parallel.setPointsIndices ( par_point_ind );
}

#endif
