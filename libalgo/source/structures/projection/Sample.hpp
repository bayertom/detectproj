// Description: Definition of the cartographic sample used for further analysis

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


#ifndef Sample_HPP
#define Sample_HPP
#include <iostream>
#include <iomanip>

#include "libalgo/source/structures/point/Node3DCartesianProjected.h"
#include "libalgo/source/structures/projection/Projection.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

#include "libalgo/source/io/Format.h"

template <typename T>
Sample <T> ::~ Sample ( )
{
        if ( proj != NULL )
                proj = NULL;
}


template <typename T>
void Sample <T> ::print ( std::ostream * output ) const
{
        //Print properties of the sample
}


template <typename T>
T Sample <T> ::getSampleCost ( const typename TAnalysisParameters <T>::TAnalysisType & analysis_type ) const
{
        //Return cost of the sample: geometric mean of all criteria
        unsigned int n = 0;
        
                /*return exp ( ( ( analysis_type.a_cnd ? ( cross_nearest_neighbour_distance_ratio_position > 0 ? log ( cross_nearest_neighbour_distance_ratio * ( float ) ( ++ n || 1 ) ) : 0.0 ) : 0.0 ) +
                               ( analysis_type.a_homt ? ( homothetic_transformation_ratio_position > 0 ? log ( homothetic_transformation_ratio * ( float ) ( ++ n || 1 ) ) : 0.0 ) : 0.0 ) +
                               ( analysis_type.a_helt ? ( helmert_transformation_ratio_position > 0 ? log ( helmert_transformation_ratio * ( float ) ( ++ n || 1 ) ) : 0.0 ) : 0.0 ) +
                               ( analysis_type.a_gn_tf ? ( gn_turning_function_ratio_position > 0 ? log ( gn_turning_function_ratio * ( float ) ( ++ n || 1 ) ) : 0.0 ) : 0.0 ) +
                               ( analysis_type.a_vd_tf ? ( voronoi_cell_turning_function_ratio_position > 0 ? log ( voronoi_cell_turning_function_ratio * ( float ) ( ++ n || 1 ) ) : 0.0 ) : 0.0 ) )
                             / ( float ) n );*/
        	
		return homothetic_transformation_ratio;
}


template <typename T>
void Sample <T> ::printSampleRatios ( const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, const unsigned int n, std::ostream * output ) const
{
        //Print results of the cartometric analysis (ratios) in line
        const char DASH[7] = "------";

        //Get category
	const char *category = (proj != NULL ? proj->getProjectionFamily() : "Unknown");

        //Print line
        *output << std::setw ( 5 ) << std::fixed << std::right;
        *output << std::setw ( 4 ) << position;

        //Print '*' for rotated sample
        if ( rotated_sample ) *output << std::setw ( 1 ) << "*";
        else *output << std::setw ( 1 ) << " ";

	*output << std::setw(7) << (proj != NULL ? proj->getName(): "xxx");
        *output << std::setw ( 7 ) << category;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << latp;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lonp;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << lat0;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lon0;

        char constant_text [255];
        Format::modScientific ( c, constant_text );
        *output << std::setw ( 9 ) << std::setprecision ( 2 ) << constant_text;

        //Print results, b-best points (remove outliers)
        *output << std::setw ( 5 ) << std::setprecision ( 0 ) << ( int ) ( 100 * k_best_points_indices.size() / ( float ) n ) << "%";

        *output << std::scientific;

        //Print results: Cross nearest distance
        if ( ( analysis_type.a_cnd ) && ( cross_nearest_neighbour_distance_ratio_position > 0 ) )
        {
                //Print number in shorten semilogaritmic format
                char cross_nearest_neighbour_distance_ratio_text [255];
                Format::modScientific ( cross_nearest_neighbour_distance_ratio, cross_nearest_neighbour_distance_ratio_text );
                *output << std::setw ( 9 ) << cross_nearest_neighbour_distance_ratio_text;
        }

        //Test has not been performed
        else *output << std::setw ( 9 ) << DASH ;

        //Print results: Homothetic transformation
        if ( ( analysis_type.a_homt ) && ( homothetic_transformation_ratio_position > 0 ) )
        {
                //Print number in shorten semilogaritmic format
                char homothetic_transformation_ratio_text [255];
                Format::modScientific ( homothetic_transformation_ratio, homothetic_transformation_ratio_text );
                *output	<< std::setw ( 9 ) << homothetic_transformation_ratio_text;
                *output << std::setw ( 5 ) << std::setw ( 4 ) << homothetic_transformation_perc_match << "%";
        }

        //Test has not been performed
        else  *output << std::setw ( 14 ) << DASH ;

        //Print results: Helmert transformation
        if ( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) )
        {
                //Print number in shorten semilogaritmic format
                char helmert_transformation_ratio_text [255];
                Format::modScientific ( helmert_transformation_ratio, helmert_transformation_ratio_text );
                *output << std::setw ( 9 ) << helmert_transformation_ratio_text;
                *output << std::setw ( 5 ) << std::setw ( 4 ) << helmert_transformation_perc_match << "%" ;
        }

        //Test has not been performed
        else  *output << std::setw ( 14 ) << DASH ;

        //Print results: Geographic network turning function
        if ( ( analysis_type.a_gn_tf ) && ( gn_turning_function_ratio_position > 0 ) )
        {
                //Print number in shorten semilogaritmic format
                char gn_turning_function_ratio_text [255];
                Format::modScientific ( gn_turning_function_ratio, gn_turning_function_ratio_text );
                *output << std::setw ( 9 ) << gn_turning_function_ratio_text;
        }

        //Test has not been performed
        else  *output << std::setw ( 9 ) << DASH ;

        //Print results: Voronoi diagram turning function
        if ( ( analysis_type.a_vd_tf ) && ( voronoi_cell_turning_function_ratio_position > 0 ) )
        {
                //Print number in shorten semilogaritmic format
                char voronoi_cell_turning_function_ratio_text [255];
                Format::modScientific ( voronoi_cell_turning_function_ratio, voronoi_cell_turning_function_ratio_text );
                *output << std::setw ( 9 ) << voronoi_cell_turning_function_ratio_text;
        }

        //Test has not been performed
        else *output << std::setw ( 9 ) << DASH ;

		
		//Print results: Voronoi diagram inner distance
		if ( ( analysis_type.a_vd_id ) && (voronoi_cell_inner_distance_ratio_position > 0 ) )
		{
				//Print number in shorten semilogaritmic format
				char voronoi_cell_inner_distance_ratio_text[255];
				Format::modScientific(voronoi_cell_inner_distance_ratio, voronoi_cell_inner_distance_ratio_text);
				*output << std::setw(9) << voronoi_cell_inner_distance_ratio_text;
		}

		//Test has not been performed
		else *output << std::setw(9) << DASH;
		

        //Print empty line
        *output << std::endl;
}


template <typename T>
void Sample <T> ::printSampleRatios2 ( const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, const unsigned int n, std::ostream * output ) const
{
        //Print results of the cartometric analysis (ratios) in line
        const char DASH[7] = "------";

        //Get category
	const char *category = (proj != NULL ? proj->getFamily() : "Unknown");

        //Print line
        *output << std::setw ( 5 ) << std::fixed << std::right;
        *output << std::setw ( 4 ) << position;
	*output << std::setw( 20) << category;
	*output << std::setw ( 8 ) << (proj != NULL ? proj->getName(): "xxx");
	*output << std::setw ( 10 ) << std::setprecision( 3 ) << R;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << latp;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lonp;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << lat0;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lon0;

        char constant_text [255];
        Format::modScientific ( c, constant_text );
        *output << std::setw ( 9 ) << std::setprecision ( 2 ) << constant_text;

	*output << std::setw(10) << std::setprecision(1) << alpha;
	*output << std::setw(13) << int(scale_homt);

	char homothetic_transformation_ratio_text[255];
	Format::modScientific(homothetic_transformation_ratio, homothetic_transformation_ratio_text);
	*output << std::setw(9) << homothetic_transformation_ratio_text;

	*output << std::setw(10) << std::setprecision(0) << iterations;
	*output << std::setw(10) << std::setprecision(0) << residual_eval;


        //Print end of line
        *output << std::endl;
}



template <typename T>
void Sample <T> ::printSamplePositions ( const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, std::ostream * output ) const
{
        //Print results of the cartometric analysis (positions) in line
        const char DASH[5] = "----";

        //Get category
	const char *category = (proj != NULL ? proj->getFamily() : "Unknown");

        *output << std::showpoint << std::fixed << std::right;

        //Print projection properties
        *output << std::setw ( 4 ) << position;

        //Print '*' for rotated sample
        if ( rotated_sample ) *output << std::setw ( 1 ) << "*";
        else *output << std::setw ( 1 ) << " ";

	*output << std::setw(7) << (proj != NULL ? proj->getName() : "xxx");
        *output << std::setw ( 7 ) << category ;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << latp;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lonp;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << lat0;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lon0;

        char constant_text [255];
        Format::modScientific ( c, constant_text );
        *output << std::setw ( 9 ) << std::setprecision ( 2 ) << constant_text;


        *output << std::scientific;

        //Print results
        ( ( analysis_type.a_cnd ) && ( cross_nearest_neighbour_distance_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 0 ) << cross_nearest_neighbour_distance_ratio_position ) : ( *output << std::setw ( 6 ) << DASH );
        ( ( analysis_type.a_homt ) && ( homothetic_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 0 ) << homothetic_transformation_ratio_position ) : ( *output << std::setw ( 6 ) << DASH );
        ( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 0 ) << helmert_transformation_ratio_position ) : ( *output << std::setw ( 6 ) << DASH );
        ( ( analysis_type.a_gn_tf ) && ( gn_turning_function_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 0 ) << gn_turning_function_ratio_position ) : ( *output << std::setw ( 6 ) << DASH );
        ( ( analysis_type.a_vd_tf ) && ( voronoi_cell_turning_function_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 0 ) << voronoi_cell_turning_function_ratio_position ) : ( *output << std::setw ( 6 ) << DASH );
		( ( analysis_type.a_vd_id) && (voronoi_cell_inner_distance_ratio_position > 0)) ? (*output << std::setw(6) << std::setprecision(0) << voronoi_cell_inner_distance_ratio_position) : (*output << std::setw(6) << DASH);

        //Print empty line
        *output << std::endl;
}


template <typename T>
void Sample <T> ::printSampleMatchedPoints ( const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> & nl_reference, const int position,
                const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, std::ostream * output ) const
{
        //Print results of the cartometric analysis: all matched points
        const char CHECK[2] = "x";
        const char DASH_SHORT[6] = "-----";
        const char DASH_LONG[31] = " -----------------------------";

        //Print projection properties
        *output << std::fixed << std::right;
        *output << std::setw ( 4 ) << position;

        if ( rotated_sample ) *output << std::setw ( 1 ) << "*";
        else *output << std::setw ( 1 ) << " ";

	*output << std::setw(6) << (proj != NULL ? proj->getName() : "xxx");
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << latp;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lonp;
        *output << std::setw ( 6 ) << std::setprecision ( 1 ) << lat0;
        *output << std::setw ( 7 ) << std::setprecision ( 1 ) << lon0 << ':';
        *output << std::endl << DASH_LONG << std::endl << std::endl;

        //Print scales and rotation for analyzed map
        *output	<< std::setw ( 20 ) << "Scale HOMT:";
        ( ( analysis_type.a_homt ) && ( homothetic_transformation_ratio_position > 0 ) ) ? ( *output  << std::setw ( 12 ) << std::setprecision ( 1 ) << scale_homt ) : ( *output << DASH_SHORT );
        *output << std::endl;

        *output	<< std::setw ( 20 ) << "Scale HELT:";
        ( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 12 ) << std::setprecision ( 1 ) << scale_helt ) : ( *output << DASH_SHORT );
        *output << std::endl;

        *output	<< std::setw ( 20 ) << "Rotation HELT:";

        if ( !rotated_sample )	( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 2 ) << rotation << " deg" ) : ( *output << DASH_SHORT );
        else
        {
                ( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 2 ) << 0.0 << " deg" ) : ( *output << DASH_SHORT );
                *output << std::endl;

                *output	<< std::setw ( 20 ) << "Old rotation HELT:";
                ( ( analysis_type.a_helt ) && ( helmert_transformation_ratio_position > 0 ) ) ? ( *output << std::setw ( 6 ) << std::setprecision ( 2 ) << rotation << " deg" ) : ( *output << DASH_SHORT );
        }

        *output << std::endl << std::endl;

        //Create header for the list
        *output	<< std::setw ( 3 ) << "#"
                << std::setw ( 6 ) << "NSING"
                << std::setw ( 6 ) << "BKEY"
                << std::setw ( 6 ) << "HOMT"
                << std::setw ( 6 ) << "HELT" << std::endl;

        //Create index variables
        unsigned int index_non_singular_points = 0, index_k_best_points = 0, index_homt = 0, index_helt = 0;
        const unsigned int n_points = nl_test.size(), n_non_singular = non_singular_points_indices.size(), n_k_best = k_best_points_indices.size(),
                           n_homt = homothetic_transformation_matched_points_indices.size(), n_helt = helmert_transformation_matched_points_indices.size();

        /*
        std::cout << n_k_best<< "  " << n_homt <<"   "  << n_helt << '\n';

        for ( int i = 0; i < n_points;i++ ) std::cout << non_singular_points_indices[i] << " ";
        std::cout <<'\n';
        for ( int i = 0; i < n_k_best;i++ ) std::cout << k_best_points_indices[i] << " ";
        std::cout <<'\n';
        for ( int i = 0; i < n_homt;i++ ) std::cout << homothetic_transformation_matched_points_indices[i] << " ";
        std::cout <<'\n';
        for ( int i = 0; i < n_homt;i++ ) std::cout << helmert_transformation_matched_points_indices[i] << " ";
        std::cout <<'\n';
        */
        //Print all items
        for ( unsigned int i = 0; i < n_points; i++ )
        {
                //Print point index and * for rotated sample
                *output << std::setw ( 3 ) << i;

                //Print non singular points
                if ( ( index_non_singular_points < n_non_singular ) && ( i == non_singular_points_indices [index_non_singular_points ] ) )
                {
                        //*output << " insb" << index_non_singular_points;
                        *output << std::setw ( 6 ) << CHECK;
                        index_non_singular_points++;
                }

                else *output << std::setw ( 6 ) << " ";


                //Print k_best points without outliers
                //*output << "f1_ " << ( i==k_best_points_indices [index_k_best_points ]) << " " <<non_singular_points_indices [ k_best_points_indices [index_k_best_points ] ];
                if ( ( index_k_best_points < n_k_best ) && ( i == non_singular_points_indices [ k_best_points_indices [index_k_best_points ] ] ) )
                {
                        //*output <<"f2 ";
                        //*output << " ikb" << index_k_best_points;
                        *output << std::setw ( 6 ) << CHECK;
                        index_k_best_points++;
                }

                else *output << std::setw ( 6 ) << " ";

                //*output << std::setw ( 6 ) << " kbest " << n_k_best;
                //*output <<"f3 " << index_homt <<"  " << homothetic_transformation_matched_points_indices[index_homt] << "  " << k_best_points_indices [ homothetic_transformation_matched_points_indices[index_homt]];
                //Print matched points: homothetic transformation
                if ( ( index_homt < n_homt ) && ( i == non_singular_points_indices [ k_best_points_indices [ homothetic_transformation_matched_points_indices[index_homt] ] ] ) )
                {
                        //*output <<"f4 ";
                        //*output << " ihot" << index_homt <<" ht " <<  homothetic_transformation_matched_points_indices[index_homt];
                        *output << std::setw ( 6 ) << CHECK;
                        index_homt++;
                }

                else *output << std::setw ( 6 ) << " ";

                //Print matched points: helmert transformation
                //*output <<"f5 ";
                if ( ( index_helt < n_helt ) && ( i == non_singular_points_indices [ k_best_points_indices [ helmert_transformation_matched_points_indices[index_helt] ] ] ) )
                {
                        // *output <<"f6 ";
                        // *output << " ihet" << index_helt;
                        *output << std::setw ( 6 ) << CHECK;
                        index_helt++;
                }

                else *output << std::setw ( 6 ) << " ";

                *output << std::endl;
        }

        *output << std::endl << std::endl;
}


#endif
