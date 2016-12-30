// Description: Create projection graticule given by the lat/lon intervals
// Support lon0 shift for the oblique aspect

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

#ifndef Graticule_HPP
#define Graticule_HPP

#include <math.h>
#include <algorithm>
#include <iterator>

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/algorithms/round/Round.h"
#include "libalgo/source/exceptions/MathException.h"


template <typename T>
void Graticule::computeGraticule ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> lat_extent, const TInterval<T> lon_extent, const T lat_step, const T lon_step, const T d_lat, const T d_lon, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point3DCartesian <T> > & meridians_proj, TVector <Parallel <T> > &parallels, TVector2D <Point3DCartesian <T> > &parallels_proj )
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
                //Create lat and lon intervals
                TList <TInterval<T>> lat_intervals, lon_intervals;
		lat_intervals.push_back(lat_extent);
		lon_intervals.push_back(lon_extent);

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
					//computeMeridians ( proj, *i_lat_intervals, *i_lon_intervals, lon_step, d_lat, alpha, meridians_temp, meridians_temp_proj, lat_error, lon_error );
					computeParallels(proj, *i_lat_intervals, *i_lon_intervals, lat_step, d_lon, alpha, parallels_temp, parallels_temp_proj, lat_error, lon_error);

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
void Graticule::computeMeridians (const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lon_step, const T d_lat, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point> & meridians_proj, T &lat_error, T &lon_error )
{
	//Set start value of the longitude as a multiplier of lon_step 
	const T lon_start = Round::roundToMultipleFloor(lon_interval.min, lon_step) + lon_step;
	const T lon_end = Round::roundToMultipleCeil(lon_interval.max, lon_step) - lon_step;
	
	T lon = std::max(std::min(lon_interval.min + GRATICULE_LAT_LON_SHIFT, MAX_LON), MIN_LON);

	try
	{
		//Create first meridian (lower bound of the interval, not lower than MIN_LON)
		//Split the interval, if a meridian is intersected by the prime meridian of the transformed system [lat_trans, lon_trans]
		computeMeridianFragment(proj, lon, lat_interval, d_lat, alpha, meridians, meridians_proj);
		
		/*
		for (TInterval <T> lat_interval_split : splitLatInterval(lat_interval, lon, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
		{
			Meridian <T> m1(lon, lat_interval_split, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> m1_proj;
			m1.project(proj, alpha, m1_proj);
			//m1.print();
			//std::cout << '\n' << m1.getLon() << " " << '\n';
			//std::copy(m1_proj.begin(), m1_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));

			//Add meridian and projected meridian to the list
			meridians.push_back(m1);
			meridians_proj.push_back(m1_proj);
		}
		*/

		//Create intermediate meridians
		//Split the interval, if a meridian is intersected by the prime meridian of the transformed system [lat_trans, lon_trans]
		for (lon = lon_start; lon <= lon_end; lon += lon_step)
		{
			computeMeridianFragment(proj, lon, lat_interval, d_lat, alpha, meridians, meridians_proj);
			/*
			for (TInterval <T> lat_interval_split : splitLatInterval(lat_interval, lon, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
			{
				Meridian <T> m2(lon, lat_interval_split, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
				TVector <Point> m2_proj;
				m2.project(proj, alpha, m2_proj);
				//m2.print();
				//std::cout << '\n' << m2.getLon() << " " << '\n';
				//std::copy(m2_proj.begin(), m2_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));

				//Add meridian and projected meridian to the list
				meridians.push_back(m2);
				meridians_proj.push_back(m2_proj);
			}
			*/
		}

		//Create last meridian (upper bound of the interval, not higher than MAX_LON)
		//Split the interval, if a meridian is intersected by the prime meridian of the transformed system [lat_trans, lon_trans]
		lon = std::max(std::min(lon_interval.max - GRATICULE_LAT_LON_SHIFT, MAX_LON), MIN_LON);

		computeMeridianFragment(proj, lon, lat_interval, d_lat, alpha, meridians, meridians_proj);
		/*
		for (TInterval <T> lat_interval_split : splitLatInterval(lat_interval, lon, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
		{
			Meridian <T> m3(lon, lat_interval_split, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> m3_proj;
			m3.project(proj, alpha, m3_proj);
			//m3.print();
			//std::cout << '\n' << m3.getLon() << " " << '\n';
			//std::copy(m3_proj.begin(), m3_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));

			//Add meridian and projected meridian to the list
			meridians.push_back(m3);
			meridians_proj.push_back(m3_proj);
		}
		*/

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
void Graticule::computeParallels ( const std::shared_ptr <Projection <T> > proj, const TInterval<T> &lat_interval, const TInterval<T> &lon_interval, const T lat_step, const T d_lon, 
	const T alpha, TVector <Parallel <T> > &parallels, TVector2D <Point> & parallels_proj, T &lat_error, T &lon_error )
{
	//Set start value of the longitude as a multiplier of lon_step 
	const T lat_start = Round::roundToMultipleFloor(lat_interval.min, lat_step) + lat_step;
	const T lat_end = Round::roundToMultipleCeil(lat_interval.max, lat_step) - lat_step;
	
	T lat = std::max(std::min(lat_interval.min + GRATICULE_LAT_LON_SHIFT, MAX_LAT), MIN_LAT);
	
	try
	{
		//Create first parallel (lower bound of the interval, nor lower than MIN_LAT)
		//Split the interval, if a parallel is intersected by the prime meridiann of the transformed system [lat_trans, lon_trans]
		computeParallelFragment(proj, lat, lon_interval, d_lon, alpha, parallels, parallels_proj);
		/*
		for (TInterval <T> lon_interval_split : splitLonInterval(lon_interval, lat, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
		{
			Parallel <T> p1(lat, lon_interval_split, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> p1_proj;
			p1.project(proj, alpha, p1_proj);

			//Add parallel and projected parallel to the list
			parallels.push_back(p1);
			parallels_proj.push_back(p1_proj);
			//p1.print();
			//std::cout << '\n' << p1.getLat() << " " << '\n';
			//std::copy(p1_proj.begin(), p1_proj.end(), std::ostream_iterator<Point3DCartesian<T> >(std::cout, " "));
		}
		*/

		//Create intermediate parallels
		//Split the interval, if a parallel is intersected by the prime meridian of the transformed system [lat_trans, lon_trans]
		for (lat = lat_start; lat <= lat_end; lat += lat_step)
		{
			computeParallelFragment(proj, lat, lon_interval, d_lon, alpha, parallels, parallels_proj);
			/*
			for (TInterval <T> lon_interval_split : splitLonInterval(lon_interval, lat, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
			{
				Parallel <T> p2(lat, lon_interval_split, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
				TVector <Point> p2_proj;
				p2.project(proj, alpha, p2_proj);
				//p2.print();
				//std::cout << '\n' << p2.getLat() << " " << '\n';
				//std::copy(p2_proj.begin(), p2_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));

				//Add parallel and projected parallel to the list
				parallels.push_back(p2);
				parallels_proj.push_back(p2_proj);
			}
			*/
		}

		//Create last parallel (upper bound of the interval, not greater than MAX_LAT)
		//Split the interval, if a parallel is intersected by the prime meridian of the transformed system [lat_trans, lon_trans]
		lat = std::max(std::min(lat_interval.max - GRATICULE_LAT_LON_SHIFT, MAX_LAT), MIN_LAT);

		computeParallelFragment(proj, lat, lon_interval, d_lon, alpha, parallels, parallels_proj);

		/*
		for (TInterval <T> lon_interval_split : splitLonInterval(lon_interval, lat, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
		{
			Parallel <T> p3(lat, lon_interval_split, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
			TVector <Point> p3_proj;
			p3.project(proj, alpha, p3_proj);
			//p3.print();
			//std::cout << '\n' << p3.getLat() << " " << '\n';
			//std::copy(p3_proj.begin(), p3_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));


			//Add parallel and projected parallel to the list
			parallels.push_back(p3);
			parallels_proj.push_back(p3_proj);
		}
		*/

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


template <typename T, typename Point>
void Graticule::computeMeridianFragment(const std::shared_ptr <Projection <T> > proj, const T lon, const TInterval <T> lat_interval, const T d_lat, const T alpha, TVector <Meridian <T> > &meridians, TVector2D <Point> & meridians_proj)
{
	//Create meridian fragment given by the latitude interval
	//A meridian is interrupted at its intersection with the prime meridian of the transfomed system [lat_trans, lon_trans]
	for (TInterval <T> lat_interval_split : splitLatInterval(lat_interval, lon, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
	{
		Meridian <T> mer(lon, lat_interval_split, d_lat, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> mer_proj;
		mer.project(proj, alpha, mer_proj);
		//mer.print();
		//std::cout << '\n' << m3.getLon() << " " << '\n';
		//std::copy(mer_proj.begin(), mer_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));

		//Add meridian and projected meridian to the list
		meridians.push_back(mer);
		meridians_proj.push_back(mer_proj);
	}
}


template <typename T, typename Point>
void Graticule::computeParallelFragment(const std::shared_ptr <Projection <T> > proj, const T lat, const TInterval <T> lon_interval, const T d_lon, const T alpha, TVector <Parallel <T> > &parallels, TVector2D <Point> & parallels_proj)
{
	//Create parallel fragment given by the longitude interval [lon_minm lon_max]
	//A parallel is interrupted at its intersection with the prime meridian of the transformed system [lat_trans, lon_trans]
	for (TInterval <T> lon_interval_split : splitLonInterval(lon_interval, lat, proj->getCartPole(), proj->getLonDir(), proj->getLon0()))
	{
		Parallel <T> par(lat, lon_interval_split, d_lon, GRATICULE_LAT_LON_SHIFT, GRATICULE_LAT_LON_SHIFT);
		TVector <Point> par_proj;
		par.project(proj, alpha, par_proj);
		//par.print();
		//std::cout << '\n' << par.getLat() << " " << '\n';
		//std::copy(par_proj.begin(), par_proj.end(), std::ostream_iterator<Point3DCartesian<T>>(std::cout, " "));


		//Add parallel and projected parallel to the list
		parallels.push_back(par);
		parallels_proj.push_back(par_proj);
	}
}




template <typename T>
TVector<TInterval<T> > Graticule::splitLatInterval(const TInterval <T> &lat_interval, const T lon, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection lon_dir, const T lon0)
{
	//Sample meridian of the created graticule in the normal aspect
	TVector<TInterval<T> > lats_split;

	const Point3DGeographic <T> p1(lat_interval.min, lon);
	const Point3DGeographic <T> p2(0.5 * (lat_interval.min + lat_interval.max), lon);
	const Point3DGeographic <T> p3(lat_interval.max, lon);
	
	//Sample central meridian of the transformed system in the oblique aspect
	//Transform it to the normal aspect
	const T lat4 = CartTransformation::latTransToLat(-80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon4 = CartTransformation::lonTransToLon(-80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//P5: (lat_trans, lon_trans) -> (lat, lon)
	const T lat5 = CartTransformation::latTransToLat(0.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon5 = CartTransformation::lonTransToLon(0.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//P6: (lat_trans, lon_trans) -> (lat, lon)
	const T lat6 = CartTransformation::latTransToLat(80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon6 = CartTransformation::lonTransToLon(80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//Create sampled points
	const Point3DGeographic <T> p4(lat4, lon4);
	const Point3DGeographic <T> p5(lat5, lon5);
	const Point3DGeographic <T> p6(lat6, lon6);

	//Intersection of both great circles
	Point3DGeographic<double> i1, i2;
	bool intersection_exists = GreatCircleIntersection::getGreatCirclePlainIntersection(p1, p2, p3, p4, p5, p6, i1, i2, pole, lon_dir);

	//Is there any intersection ?
	if (intersection_exists)
	{
		//Get both intersections
		const T lat_inters1 = i1.getLat();
		const T lat_inters2 = i2.getLat();
		const T lon_inters1 = i1.getLon();
		const T lon_inters2 = i2.getLon();

		//Is an intersection point i1 latitude inside the lat interval? Create 2 intervals.
		if ((lat_inters1 > lat_interval.min) && (lat_inters1 < lat_interval.max) && (fabs(lon - lon_inters1) < ANGLE_ROUND_ERROR))
		{
			lats_split.push_back(TInterval <T> { lat_interval.min, lat_inters1 - GRATICULE_LAT_LON_SHIFT });
			lats_split.push_back(TInterval <T> { lat_inters1 + GRATICULE_LAT_LON_SHIFT, lat_interval.max });
			
			return lats_split;
		}
		
		//Is an intersection point i1 latitude inside the lat interval? Create 2 intervals.
		else if ((lat_inters2 > lat_interval.min) && (lat_inters2 < lat_interval.max) && (fabs(lon - lon_inters2) < ANGLE_ROUND_ERROR))
		{		
			lats_split.push_back(TInterval <T> { lat_interval.min, lat_inters2 - GRATICULE_LAT_LON_SHIFT});
			lats_split.push_back(TInterval <T> { lat_inters2 + GRATICULE_LAT_LON_SHIFT, lat_interval.max });

			return lats_split;
		}
	}

	//No intersection, no split
	lats_split.push_back(TInterval <T> { lat_interval.min, lat_interval.max });
	
	return lats_split;
}


template <typename T>
TVector<TInterval<T> > Graticule::splitLonInterval(const TInterval <T> &lon_interval, const T lat, const Point3DGeographic <T> &pole, const TTransformedLongitudeDirection lon_dir, const T lon0)
{
	//Sample parallel of the created graticule in the normal aspect
	TVector<TInterval<T> > lons_split;

	const Point3DGeographic <T> p1(lat, lon_interval.min + 1.0);
	const Point3DGeographic <T> p2(lat, 0.5 * (lon_interval.min + lon_interval.max));
	const Point3DGeographic <T> p3(lat, lon_interval.max - 1.0);

	//Sample central meridian of the transformed system in the oblique aspect
	//Transform it to the normal aspect
	const T lat4 = CartTransformation::latTransToLat(-80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon4 = CartTransformation::lonTransToLon(-80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//P5: (lat_trans, lon_trans) -> (lat, lon)
	const T lat5 = CartTransformation::latTransToLat(0.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon5 = CartTransformation::lonTransToLon(0.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//P6: (lat_trans, lon_trans) -> (lat, lon)
	const T lat6 = CartTransformation::latTransToLat(80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);
	const T lon6 = CartTransformation::lonTransToLon(80.0, lon0, pole.getLat(), pole.getLon(), lon_dir);

	//Create sampled points
	const Point3DGeographic <T> p4(lat4, lon4);
	const Point3DGeographic <T> p5(lat5, lon5);
	const Point3DGeographic <T> p6(lat6, lon6);

	//Intersection of both great circles
	Point3DGeographic<double> i1, i2;
	bool intersection_exists = GreatCircleIntersection::getGreatCirclePlainIntersection(p1, p2, p3, p4, p5, p6, i1, i2, pole, lon_dir);

	//Is there any intersection ?
	if (intersection_exists)
	{
		//Sort both intersections according to lon
		const T lon_inters1 = std::min(i1.getLon(), i2.getLon());
		const T lon_inters2 = std::max(i1.getLon(), i2.getLon());

		//Both intersections inside the lon interval? Create 3 intervals.
		if ((lon_inters1 > lon_interval.min && lon_inters1 < lon_interval.max) && (lon_inters2 > lon_interval.min && lon_inters2 < lon_interval.max))
		{
			lons_split.push_back(TInterval <T> { lon_interval.min, lon_inters1 - GRATICULE_LAT_LON_SHIFT});
			lons_split.push_back(TInterval <T> { lon_inters1 + GRATICULE_LAT_LON_SHIFT, lon_inters2 - GRATICULE_LAT_LON_SHIFT});
			lons_split.push_back(TInterval <T> { lon_inters2 + GRATICULE_LAT_LON_SHIFT, lon_interval.max });

			return lons_split;
		}

		//Only the first intersection lon1 inside the lon interval. Create 2 intervals.
		else if (lon_inters1 > lon_interval.min && lon_inters1 < lon_interval.max)
		{
			lons_split.push_back(TInterval <T> { lon_interval.min, lon_inters1 - GRATICULE_LAT_LON_SHIFT});
			lons_split.push_back(TInterval <T> { lon_inters1 + GRATICULE_LAT_LON_SHIFT, lon_interval.max });

			return lons_split;
		}

		//Only the second intersection lon2 inside the lon interval. Create 2 intervals.
		else if (lon_inters2 > lon_interval.min && lon_inters2 < lon_interval.max)
		{
			lons_split.push_back(TInterval <T> { lon_interval.min, lon_inters2 - GRATICULE_LAT_LON_SHIFT});
			lons_split.push_back(TInterval <T> { lon_inters2 + GRATICULE_LAT_LON_SHIFT, lon_interval.max });

			return lons_split;
		}
	}

	//No intersection, no split
	lons_split.push_back(TInterval <T> { lon_interval.min, lon_interval.max });

	return lons_split;
}


#endif
