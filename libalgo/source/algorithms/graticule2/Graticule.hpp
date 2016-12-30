// Description: Create projection graticule given by the lat/lon intervals

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
// 2

#ifndef Graticule_HPP
#define Graticule_HPP

#include <math.h>
#include <algorithm>
#include <iterator>

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/algorithms/round/Round.h"
#include "libalgo/source/exceptions/MathException.h"


template <typename T>
void Graticule::createGraticule ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> lat_extent, const TInterval<T> lon_extent, const T lat_step, const T lon_step, const T d_lat, const T d_lon, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point3DCartesian <T> > & meridians_proj, TVector <Parallel <T> > &parallels, TVector2D <Point3DCartesian <T> > &parallels_proj )
{
        //Compute graticule given by meridians and parallels
        TInterval <T> lat_interval_part, lon_interval_part;

        //Parameters of the projection
        const T latp = proj->getCartPole().getLat();
	const T lonp = proj->getCartPole().getLon();
	const T lon0 = proj->getLon0();

        //Get size of the interval
        const T lat_interval_width = lat_extent.max - lat_extent.min;
        const T lon_interval_width = lon_extent.max - lon_extent.min;

        //Generate meridians and parallels
        if ( lat_step <= lat_interval_width || lon_step <= lon_interval_width )
        {
                //Create lat and lon intervals to avoid basic singular points
                TList <TInterval<T>> lat_intervals, lon_intervals;
                createLatIntervals ( lat_extent, latp, lat_intervals );

		const T lon_break = CartTransformation::redLon0(lonp, -lon0);
		createLonIntervals(lon_extent, lon_break, lon_intervals);

                //Process all meridians and parallels: automatic detection of additional singular points
                unsigned int split_amount = 0;

		//Process all lat intervals
                for ( typename TList <TInterval<T>>::iterator i_lat_intervals = lat_intervals.begin(); i_lat_intervals != lat_intervals.end(); i_lat_intervals ++ )
                {
			//std::cout << "lat: <" << i_lat_intervals->min << ", " << i_lat_intervals->max << "> \n";
                        //Process all lon intervals
			for ( typename TList <TInterval<T>>::iterator i_lon_intervals = lon_intervals.begin() ; i_lon_intervals != lon_intervals.end(); )
                        {
				//std::cout << "   lon: <" << i_lon_intervals->min << ", " << i_lon_intervals->max << "> \n";
                                //Create temporary meridians and parallels
				TVector <Meridian <T> > meridians_temp;
				TVector <Parallel <T> > parallels_temp;
				TVector2D <Point3DCartesian<T> > meridians_temp_proj,  parallels_temp_proj;

                                //Compute meridians and parallels
                                T lat_error = 0.0, lon_error = 0.0;

				try
				{
					//Compute meridians and parallels
					createMeridians ( proj, *i_lat_intervals, *i_lon_intervals, lon_step, d_lat, alpha, meridians_temp, meridians_temp_proj, lat_error, lon_error );
					createParallels(proj, *i_lat_intervals, *i_lon_intervals, lat_step, d_lon, alpha, parallels_temp, parallels_temp_proj, lat_error, lon_error);

					//Copy temporary meridians and parallels to the output data structure
					std::copy(meridians_temp.begin(), meridians_temp.end(), std::inserter(meridians, meridians.begin()));
					std::copy(parallels_temp.begin(), parallels_temp.end(), std::inserter(parallels, parallels.begin()));

					//Copy temporary projected meridians and parallels to the output data structure
					std::copy(meridians_temp_proj.begin(), meridians_temp_proj.end(), std::inserter(meridians_proj, meridians_proj.begin()));
					std::copy(parallels_temp_proj.begin(), parallels_temp_proj.end(), std::inserter(parallels_proj, parallels_proj.begin()));

                                        //Increment lon intervals only if computation was successful
                                        i_lon_intervals ++;
                                }
				
                                //Exception
                                catch ( MathException <T> &error )
                                {
                                        //Too many splits, projection is suspected, stop computations
                                        if ( split_amount > 100 )
                                        {
                                                meridians.clear(); parallels.clear();
						meridians_proj.clear(); parallels_proj.clear();
                                                return;
                                        }

					//Empty lat interval, delete
					if (fabs((*i_lat_intervals).max - (*i_lat_intervals).min) < 10 * GRATICULE_LAT_LON_SHIFT)
					{
						i_lat_intervals = lat_intervals.erase(i_lat_intervals);
						++i_lat_intervals;
					}

					//Lat value is lower bound: shift lower bound
					else if ((fabs((*i_lat_intervals).min - lat_error) <  2 * GRATICULE_LAT_LON_SHIFT) && (fabs((*i_lat_intervals).max - (*i_lat_intervals).min) > 10 * GRATICULE_LAT_LON_SHIFT))
					{
						(*i_lat_intervals).min += 2 * GRATICULE_LAT_LON_SHIFT;
					}

					//Lat value is upper bound: shift upper bound
					else if ((fabs((*i_lat_intervals).max - lat_error) <  2 * GRATICULE_LAT_LON_SHIFT) && (fabs((*i_lat_intervals).max - (*i_lat_intervals).min) > 10 * GRATICULE_LAT_LON_SHIFT))
					{
						(*i_lat_intervals).max -= 2 * GRATICULE_LAT_LON_SHIFT;
					}

					//Lat value inside interval: split intervals
					else if ((((*i_lat_intervals).min) < lat_error) && ((*i_lat_intervals).max > lat_error))
					{
						splitIntervals(lat_intervals, i_lat_intervals, lat_error);

						//Increment split amount
						split_amount++;
					}

					//Empty lon interval, delete
					if (fabs((*i_lon_intervals).max - (*i_lon_intervals).min) < 10 * GRATICULE_LAT_LON_SHIFT)
					{
						i_lon_intervals = lon_intervals.erase(i_lon_intervals);
						++i_lon_intervals;
					}

					//Lon value is lower bound : shift lower bound
					else if ((fabs((*i_lon_intervals).min - lon_error) <  2 *GRATICULE_LAT_LON_SHIFT) && (fabs((*i_lon_intervals).max - (*i_lon_intervals).min) > 10 * GRATICULE_LAT_LON_SHIFT))
					{
						(*i_lon_intervals).min += 2 * GRATICULE_LAT_LON_SHIFT;
					}
					
					//Lon value is upper bound: shift upper bound
					else if ((fabs((*i_lon_intervals).max - lon_error) <  2 * GRATICULE_LAT_LON_SHIFT) && (fabs((*i_lon_intervals).max - (*i_lon_intervals).min) > 10 * GRATICULE_LAT_LON_SHIFT))
					{
						(*i_lon_intervals).max -= 2 * GRATICULE_LAT_LON_SHIFT;
					}

                                        //Lon value inside interval: split intervals
                                        else if (((( *i_lon_intervals ).min ) < lon_error ) && (( *i_lon_intervals ).max > lon_error ))
                                        {
                                                splitIntervals ( lon_intervals, i_lon_intervals, lon_error );

						//Increment split amount
						split_amount++;
                                        }

					//std::cout << (*i_lat_intervals).min << "   " << (*i_lat_intervals).max << '\n';
					//std::cout << (*i_lon_intervals).min << "   " << (*i_lon_intervals).max << '\n';
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
void Graticule:: createLatIntervals ( const TInterval <T> &lat_extent, const T latp, TList <TInterval<T>> &lat_intervals )
{
        //Split the geographic extent into several sub-intervals
        TInterval <T> lat_interval_part;

        //First lat interval ( -90, latp )
        lat_interval_part.min = std::max ( MIN_LAT, lat_extent.min ); 
	lat_interval_part.max = std::min ( latp, lat_extent.max );

        //Test interval and add to the list
        if ( lat_interval_part.max > lat_interval_part.min )
                lat_intervals.push_back ( lat_interval_part );

        //Second lat interval ( latp, 90 )
        lat_interval_part.min = std::max ( latp, lat_extent.min );
	lat_interval_part.max = std::min ( lat_extent.max, MAX_LAT );

        //Test interval and add to the list
        if ( lat_interval_part.max > lat_interval_part.min )
                lat_intervals.push_back ( lat_interval_part );
}


template <typename T>
void Graticule::createLonIntervals ( const TInterval <T> &lon_extent, const T lonp, TList <TInterval<T>> &lon_intervals )
{
        //Split the geographic extent into several sub-intervals
        TInterval <T> lon_interval_part;

        //Split given interval into sub intervals: lonp >=0
        if ( lonp >= 0 )
        {
                //First interval <-180, lonp - 180)
                lon_interval_part.min = std::max ( MIN_LON, lon_extent.min ); 
		lon_interval_part.max =  std::min ( lonp - 180.0, lon_extent.max );

                //Test first interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );

                //Second interval (lonp - 180, lonp>
                lon_interval_part.min = std::max ( lonp - 180.0, lon_extent.min ); 
		lon_interval_part.max = std::min ( lonp, lon_extent.max );

                //Test second interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );

                //Third interval (lonp, MAX_LON>
                lon_interval_part.min = std::max ( lonp, lon_extent.min ); 
		lon_interval_part.max = std::min ( MAX_LON, lon_extent.max );

                //Test third interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );
        }

        //Split given interval into sub intervals: lonp < 0
        else
        {
                //First interval <-180, lonp >
                lon_interval_part.min = std::max ( MIN_LON, lon_extent.min ); 
		lon_interval_part.max = std::min ( lonp, lon_extent.max );

                //Test first interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );

                //Second interval <lonp, lonp + 180>
                lon_interval_part.min = std::max ( lonp, lon_extent.min ); 
		lon_interval_part.max = std::min ( lonp + 180.0, lon_extent.max );

                //Test second interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );

                //Third interval (lonp + 180, 180>
                lon_interval_part.min = std::max ( lonp + 180.0, lon_extent.min ); 
		lon_interval_part.max = std::min ( MAX_LON, lon_extent.max );

                //Test third interval and add to the list
                if ( lon_interval_part.max > lon_interval_part.min )
                        lon_intervals.push_back ( lon_interval_part );
        }
}


template <typename T>
void Graticule::splitIntervals ( TList <TInterval<T>> &intervals, typename TList <TInterval<T>>::iterator i_intervals, const T error )
{
        //Split one interval in 2
        const T max_val = ( *i_intervals ).max;

        //Resize old interval <min, max> to <min, error)
        ( *i_intervals ).max = error;

        //Create new interval (error, max>
	TInterval <T> interval_temp{ error, max_val };

        //Add the new interval to the list
        intervals.insert ( i_intervals, interval_temp );

        //Move iterator to the previous item
        i_intervals--;
}


template <typename T, typename Point>
void Graticule::createMeridians (const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lon_step, const T d_lat, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point> & meridians_proj, T &lat_error, T &lon_error )
{
	//Set start value of the longitude as a multiplier of lon_step 
	const T lon_start = Round::roundToMultipleFloor(lon_interval.min, lon_step) + lon_step;
	const T lon_end = Round::roundToMultipleCeil(lon_interval.max, lon_step) - lon_step;
	T lon = std::max(std::min(lon_interval.min + GRATICULE_LAT_LON_SHIFT, MAX_LON), MIN_LON);

	try
	{
		//Create first meridian (lower bound of the interval, not lower than MIN_LON)
		Meridian <T> m1(lon, lat_interval, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> m1_proj;
		m1.project(proj, alpha, m1_proj);

		//Add meridian and projected meridian to the list
		meridians.push_back(m1);
		meridians_proj.push_back(m1_proj);

		//Create intermediate meridians
		for (lon = lon_start; lon <= lon_end; lon += lon_step)
		{
			Meridian <T> m2(lon, lat_interval, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> m2_proj;
			m2.project(proj, alpha, m2_proj);

			//Add meridian and projected meridian to the list
			meridians.push_back(m2);
			meridians_proj.push_back(m2_proj);
		}

		//Create last meridian (upper bound of the interval, not higher than MAX_LON)
		lon = std::max(std::min(lon_interval.max - GRATICULE_LAT_LON_SHIFT, MAX_LON), MIN_LON);

		Meridian <T> m3(lon, lat_interval, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> m3_proj;
		m3.project(proj, alpha, m3_proj);

		//Add meridian and projected meridian to the list
		meridians.push_back(m3);
		meridians_proj.push_back(m3_proj);
	}

	//Throw math exception: get error values
	catch (MathException <T> &error)
	{
		//Get possible error values
		lat_error = error.getArg();
		lon_error = lon;

		//Throw exception
		throw;
	}

	//Throw other exception: error in parser or incorrect equation
	catch (Exception &error)
	{
		throw;
	}
}


template <typename T, typename Point>
void Graticule::createParallels ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lat_step, const T d_lon, const T alpha, TVector <Parallel <T> > &parallels, TVector2D <Point> & parallels_proj, T &lat_error, T &lon_error )
{
	//Set start value of the longitude as a multiplier of lon_step 
	const T lat_start = Round::roundToMultipleFloor(lat_interval.min, lat_step) + lat_step;
	const T lat_end = Round::roundToMultipleCeil(lat_interval.max, lat_step) - lat_step;
	T lat = std::max(std::min(lat_interval.min + GRATICULE_LAT_LON_SHIFT, MAX_LAT), MIN_LAT);

	try
	{
		//Create first parallel (lower bound of the interval, nor lower than MIN_LAT)
		Parallel <T> p1(lat, lon_interval, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> p1_proj;
		p1.project(proj, alpha, p1_proj);

		//Add parallel and projected parallel to the list
		parallels.push_back(p1);
		parallels_proj.push_back(p1_proj);

		//Create intermediate parallels
		for (lat = lat_start; lat <= lat_end; lat += lat_step)
		{
			Parallel <T> p2(lat, lon_interval, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> p2_proj;
			p2.project(proj, alpha, p2_proj);

			//Add parallel and projected parallel to the list
			parallels.push_back(p2);
			parallels_proj.push_back(p2_proj);
		}

		//Create last parallel (upper bound of the interval, not greater than MAX_LAT)
		lat = std::max(std::min(lat_interval.max - GRATICULE_LAT_LON_SHIFT, MAX_LAT), MIN_LAT);
		Parallel <T> p3(lat, lon_interval, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> p3_proj;
		p3.project(proj, alpha, p3_proj);

		//Add parallel and projected parallel to the list
		parallels.push_back(p3);
		parallels_proj.push_back(p3_proj);
	}

	//Throw math exception: get error values
	catch (MathException <T> &error)
	{
		//Get possible error values
		lat_error = lat;
		lon_error = error.getArg();

		//Throw exception
		throw;
	}

	//Throw other exception: error in parser or incorrect equation
	catch (Exception &error)
	{
		throw;
	}
}


#endif
