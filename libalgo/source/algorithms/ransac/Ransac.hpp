// Description: RANSAC algorithm implementation

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


#ifndef Ransac_HPP
#define Ransac_HPP

#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <vector>

#include "libalgo/source/structures/point/Point3DCartesian.h"

#include "libalgo/source/structures/line/Meridian.h"
#include "libalgo/source/structures/line/Parallel.h"
#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/algorithms/hash/Hash.h"

#include "libalgo/source/comparators/sortFittingLinesByHash.h"
#include "libalgo/source/comparators/sortPointsIndicesByX.h"


template <typename Point>
bool Ransac::ransacFitLine ( const Container <Point> &input, TFittingLine <Point> &acceptable_solution, const typename Point::Type acceptable_error, const bool find_best, const bool print_message, const bool print_exception, std::ostream * output )
{
        //Perform RANSAC algorithm: fit regression line through the dataset
        //Options:
        //    find_acceptable = find best fitting regression line
        //   !find_acceptable = find first acceptable regression line
        //

        //Set parameters of the algorithm
        const unsigned int n = input.size();					//Total points of the input data set
        const unsigned int max_iterations = 10 * n ;				//Total iterations

        bool status = false;							//Has already any solution been found?
        unsigned short iterations = 0;						//Total iterations under input data set
        typename Point::Type best_found_solution_error = acceptable_error;	//Error of the best found solution

        //Repeat RANSAC algorithm
        while ( iterations ++ < max_iterations )
        {
                //Maybe indices and solution
                std::vector <bool> maybe_covered_indices ( n, 0 );
                TFittingLine <Point> maybe_solution ( input ), maybe_solution_old ( input );

                //Get 2 randomly generated different indices that are not present in vector
                int index2 = 0;
                int index1 = rand() % n ;

                do
                {
                        index2 = rand() % n;
                }
                while ( index1 == index2 );


                //Add those indices into maybe solution
                maybe_solution.points_indices.insert ( index1 );
                maybe_solution.points_indices.insert ( index2 );

                //Check vertices as coveres
                maybe_covered_indices[index1] = 1;
                maybe_covered_indices[index2] = 1;

                //Compute regression line
                LeastSquaresFitting::fitLine ( input, maybe_solution );

                //Process all indices that are in maybe solution
                for ( unsigned int i = 0; i < n; i++ )
                {
                        if ( !maybe_covered_indices[i] )
                        {
                                //Remember old solution
                                maybe_solution_old = maybe_solution;

                                //Add point to maybe solution
                                maybe_solution.points_indices.insert ( i );

                                //Compute regression line
                                LeastSquaresFitting::fitLine ( input, maybe_solution );

                                //Compute max_point error
                                const typename Point::Type max_point_error = acceptable_error / sqrt ( ( double ) maybe_solution.points_indices.size() );

                                //Does not worsen the result too much? Test difference between old and actual solution
                                if ( fabs ( maybe_solution.error - maybe_solution_old.error ) > max_point_error )
                                {
                                        //Assign old solution
                                        maybe_solution = maybe_solution_old;
                                }

                                //Set point as used
                                maybe_covered_indices[i] = 1;
                        }
                }


                //Test maybe solution: does it satisfies the condition of the 4 points
                if ( maybe_solution.points_indices.size() >= RANSAC_MIN_LINE_POINTS )
                {
                        //We found a better solution
                        if ( maybe_solution.error < best_found_solution_error )
                        {
                                //Assign error
                                best_found_solution_error = maybe_solution.error;

                                //Assign found solution to the acceptable solution
                                acceptable_solution = maybe_solution;

                                //Compute hash for acceptable solution
                                acceptable_solution.hash_val = Hash::hashXY ( acceptable_solution.alpha, acceptable_solution.b, 1, 1 );

                                //Find first acceptable solution and stop processing
                                if ( !find_best )
                                {
                                        return true;
                                }

                                //Set status = true, some solution has been found
                                status = true;
                        }
                }
        }

        //Has already any solution been found
        return status;

}


template <typename Point>
void Ransac::ransacFitAllLines ( const Container <Point> &input, typename TRansacResults <Point> ::Type & ransac_results, const typename Point::Type acceptable_error, const bool print_message, const bool print_exception,  std::ostream * output )
{
        //Fit all regression lines through the dataset
        try
        {
                unsigned int max_loops = 10 * input.size() ;

                if ( print_message )
                {
                        *output << "> Starting RANSAC, detection of lines... " << std::endl ;
                }

                //Initialize random number generator
                srand ( ( unsigned ) time ( 0 ) );

                //Fit "all" regression lines
                for ( unsigned int loops = 0; loops < max_loops; loops++ )
                {
                        TFittingLine <Point> acceptable_solution ( input );

                        //Perform RANSAC: any first acceptable solution has been found
                        if ( ransacFitLine ( input, acceptable_solution, acceptable_error, false ) )
                        {
                                //Print . for every  10 th found item
                                if ( ( int ) ( 100 * ( 0.01 * loops ) ) == ( int ) ( 100 * ( ( int ) ( 0.01 * loops ) ) ) )
                                {
                                        std::cout.flush();
                                        std::cout << ".";
                                }

                                //Test if such solution is stored in the map
                                typename TRansacResults <Point>::Type ::iterator i_found = ransac_results.find ( acceptable_solution );

                                //List of results does not contain acceptable solution
                                if ( i_found != ransac_results.end() )
                                {
                                        //Acceptable solution with the same hash has been found
                                        if ( acceptable_solution.hash_val == ( *i_found ).hash_val )
                                        {
                                                //Actual acceptable solution is better than stored solution (less error) or has more points and same error
                                                if ( acceptable_solution.error < ( *i_found ).error  ||
                                                                acceptable_solution.error == ( ( *i_found ).error ) && ( acceptable_solution.points_indices.size() < ( *i_found ).points_indices.size() ) )
                                                {
                                                        //Erase found item
                                                        ransac_results.erase ( i_found );

                                                        //Add new item with better cost
                                                        ransac_results.insert ( acceptable_solution );
                                                }
                                        }
                                }

                                //List of results contains acceptable solution
                                else
                                {
                                        ransac_results.insert ( acceptable_solution );
                                }
                        }
                }

                //Print message
                if ( print_message )
                {
                        *output << ransac_results.size() << " lines have been found... Completed. \n" << std::endl << std::endl;
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                if ( print_exception )
                {
                        error.printException();
                }

                *output << "RANSAC cancelled..." << std::endl;

                throw;
        }
}


template <typename T>
void Ransac::ransacFitMeridiansAndParallels ( const Container <Point3DGeographic <T> *> &pl_geographic,  typename TMeridiansList <T> ::Type & meridians,  typename TParallelsList <T> ::Type & parallels, const T acceptable_error, const T angle_tolerance,
                const bool print_message, const bool print_exception, std::ostream * output )
{
        //Fit all meridians and parallels using RANSAC algorithm
        try
        {
                //Results of the detection
                typename TRansacResults <Point3DCartesian <T> > ::Type ransac_results;

                if ( print_message )
                {
                        std::cout << ">> Detection of meridians and parallels... \n" ;
                }

                //Create tempory Point3DCartesian list
                Container <Point3DCartesian <T> > pl_temp;

                //Load (lat, lon) coordinates, different orientation from (x, y)
                for ( unsigned int i = 0; i < pl_geographic.size(); i++ )
                {
                        pl_temp.push_back ( Point3DCartesian <T> ( pl_geographic [i]->getLon() , pl_geographic [i]->getLat() ) );
                }

                //Perforfm RANSAC under geographic points
                ransacFitAllLines ( pl_temp, ransac_results, acceptable_error, !print_message, print_exception, output );

                //Store only  meridians and parallels from RANSAC results
                for ( typename TRansacResults <Point3DCartesian <T> >::Type ::iterator i_ransac_results = ransac_results.begin(); i_ransac_results != ransac_results.end(); ++i_ransac_results )
                {
                        //Does the detected line represent a meridian? : orientation (lat, lon) = orientation (x, y) - PI/2
                        if ( fabs ( i_ransac_results->alpha - MAX_LAT ) <= angle_tolerance )
                        {
                                //Create meridian
                                Meridian <T> m ( *i_ransac_results );

                                //Add meridian to the list;
                                meridians.insert ( m );
                        }

                        //Does the detected line represent a parallel?: : orientation (lat, lon) = orientation (x, y) - PI/2
                        else if ( fabs ( i_ransac_results->alpha ) <= angle_tolerance )
                        {
                                //Create new parallel
                                Parallel <T> p ( *i_ransac_results );

                                //Add parallel to the list;
                                parallels.insert ( p );
                        }
                }

                //Print message
                if ( print_message )
                {
                        *output << '\n' << meridians.size() << " meridians and "  << parallels.size() << " parallels have been found... Completed." << std::endl;
                }
        }


        //Throw exception
        catch ( Exception & error )
        {
                if ( print_exception )
                {
                        error.printException();
                }

                *output << "Can not detect meridians and parallels..." << std::endl;

                throw;
        }
}


#endif
