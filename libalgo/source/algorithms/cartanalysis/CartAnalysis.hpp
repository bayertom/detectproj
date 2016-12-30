// Description: Performs cartometric analysis (i.e. detection of the cartographic projection)

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


#ifndef CartAnalysis_HPP
#define CartAnalysis_HPP

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <typeinfo>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Node3DCartesianProjected.h"
#include "libalgo/source/structures/projection/Sample.h"

#include "libalgo/source/algorithms/cartdistortion/CartDistortion.h"
#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"
#include "libalgo/source/algorithms/geneticalgorithms/DifferentialEvolution.h"
#include "libalgo/source/algorithms/transformation/HomotheticTransformation2D.h"
#include "libalgo/source/algorithms/transformation/HelmertTransformation2D.h"
#include "libalgo/source/algorithms/angle3points/Angle3Points.h"
#include "libalgo/source/algorithms/facearea/FaceArea.h"

#include "libalgo/source/algorithms/turningfunction/TurningFunction.h"
#include "libalgo/source/algorithms/nndistance/NNDistance.h"
#include "libalgo/source/algorithms/innerdistance/InnerDistance.h"
#include "libalgo/source/algorithms/voronoi2D/Voronoi2D.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjJ.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjJ2.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjJ3.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjJ4.h"

#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV2.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV3.h"
#include "libalgo/source/algorithms/nonlinearleastsquares/FAnalyzeProjV4.h"

#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjVS.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV2S.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV3S.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV4S.h"

#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjVS2.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV2S2.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV3S2.h"
#include "libalgo/source/algorithms/simplexmethod/FAnalyzeProjV4S2.h"

#include "libalgo/source/algorithms/geneticalgorithms/FAnalyzeProjVDE.h"
#include "libalgo/source/algorithms/geneticalgorithms/FAnalyzeProjV2DE.h"
#include "libalgo/source/algorithms/geneticalgorithms/FAnalyzeProjV3DE.h"
#include "libalgo/source/algorithms/geneticalgorithms/FAnalyzeProjV4DE.h"

#include "libalgo/source/algorithms/nonlinearleastsquares/NonLinearLeastSquares.h"
#include "libalgo/source/algorithms/simplexmethod/SimplexMethod.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"
#include "libalgo/source/algorithms/randompermutation/RandomPermutation.h"

#include "libalgo/source/comparators/sortPointsByID.h"
#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByLat.h"
#include "libalgo/source/comparators/sortPointsByLon.h"
#include "libalgo/source/comparators/sortProjectionPolePositionsByLat.h"
#include "libalgo/source/comparators/sortProjectionPolePositionsByCompCriterium.h"

#include "libalgo/source/comparators/removeUnequalMeridianParallelPointIndices.h"
#include "libalgo/source/comparators/removeProjectionPolePositions.h"
#include "libalgo/source/comparators/findMeridianParallelPointIndices.h"
#include "libalgo/source/comparators/getSecondElementInPair.h"

#include "libalgo/source/comparators/sortSamplesByCrossNearestNeighbourDistanceRatio.h"
#include "libalgo/source/comparators/sortSamplesByAverageNearestNeighbourDistanceRatio.h"
#include "libalgo/source/comparators/sortSamplesByHomotheticTransformationRatio.h"
#include "libalgo/source/comparators/sortSamplesByHelmertTransformationRatio.h"
#include "libalgo/source/comparators/sortSamplesByAngularDifferencesRatio.h"
#include "libalgo/source/comparators/sortSamplesByGNTurningFunctionRatio.h"
#include "libalgo/source/comparators/sortSamplesByNNNGraphRatio.h"
#include "libalgo/source/comparators/sortSamplesBySphereOfInfluenceGraphRatio.h"
#include "libalgo/source/comparators/sortSamplesByVoronoiCellAreaLengthRatio.h"
#include "libalgo/source/comparators/sortSamplesByVoronoiCellTurningFunctionRatio.h"
#include "libalgo/source/comparators/sortSamplesByVoronoiCellTARRatio.h"
#include "libalgo/source/comparators/sortSamplesByVoronoiCellInnerDistanceRatio.h"
#include "libalgo/source/comparators/sortSamplesByAllRatios.h"

#include "libalgo/source/exceptions/MathException.h"

#include "libalgo/source/io/DXFExport.h"

using namespace MatrixOperations;

template <typename T>
void CartAnalysis::computeAnalysisForAllSamplesGS(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
                                                  typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
                                                  unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
    //Find mimum using the global jitterpling of the function

    //Total computed analysis ( successful + thrown by the heuristic )
    total_created_or_thrown_samples = 0;

    //Total successfully computed analysis for one cartographic projection
    unsigned int total_created_and_analyzed_samples_projection = 0;

    //Create sample for analyzed projection from command line and set flag for this sample
    if (analysis_parameters.analyzed_projections.size() > 0)
    {
        //Analyze all projections specified in command line
        for (typename TItemsList<Projection <T> *> ::Type::iterator i_projections = analysis_parameters.analyzed_projections.begin(); i_projections != analysis_parameters.analyzed_projections.end(); ++i_projections)
        {
            //Get analyzed projection
            Projection <T> *analyzed_proj = *i_projections;

            //List of points using new central meridian redefined in projection file
            Container <Point3DGeographic <T> *> pl_reference_red;

            //Reduce lon using a new central meridian redefined in projection file, if necessary
            if ((*i_projections)->getLon0() != 0.0) redLon(pl_reference, (*i_projections)->getLon0(), pl_reference_red);

            //Set pointer to processed file: reduced or non reduced
            Container <Point3DGeographic <T> *> * p_pl_reference = ((*i_projections)->getLon0() == 0.0 ? &pl_reference : &pl_reference_red);

            //Create temporary containers for non singular points
            Container <Node3DCartesian <T> *> nl_test_non_sing;
            Container <Point3DGeographic <T> *> pl_reference_non_sing;

            typename TDevIndexPairs <T>::Type non_singular_pairs;
            TIndexList non_singular_points;

            //Initialize non singular indices
            for (unsigned int i = 0; i < p_pl_reference->size(); i++) non_singular_points.push_back(i);

            //Set pointer to processed file: with or wihout singular points
            Container <Node3DCartesian <T> *> *p_nl_test = &nl_test;

            //Remove singular points to prevent throwing a sample
            bool singular_points_found = false;
            removeSingularPoints(*p_nl_test, *p_pl_reference, *i_projections, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

            //Some singular points have been found
            if (nl_test.size() != nl_test_non_sing.size())
            {
                //Set pointers to files without singular points
                p_nl_test = &nl_test_non_sing;
                p_pl_reference = &pl_reference_non_sing;

                //Set flag to true, used for a sample using non-singular sets
                singular_points_found = true;

                //Correct meridians and parallels
                correctMeridiansAndParrallels <T>(meridians, parallels, non_singular_pairs);

                //Convert non singular pairs to index list: indices will be printed in output
                non_singular_points.clear();
                std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());
            }

            //Compute analysis
            try
            {
                Sample <T> analyzed_sample;
                (void) computeAnalysisForOneSample(*p_nl_test, *p_pl_reference, meridians, parallels, faces_test, analyzed_proj, analysis_parameters, analyzed_sample, singular_points_found, total_created_and_analyzed_samples_projection, output);

                //Add result to the list
                if (total_created_and_analyzed_samples_projection) sl.push_back(analyzed_sample);
            }

            //Throw exception
            catch (Exception & error)
            {
                if (analysis_parameters.print_exceptions)
                    error.printException();
            }

            //Sample with analyzed projection has been successfully created (not thrown by the heuristic)
            if (total_created_and_analyzed_samples_projection > 0) sl[sl.size() - 1].setAnalyzedProjectionSample(true);
        }

        if (total_created_and_analyzed_samples_projection == 0) throw BadDataException("BadDataException: no analyzed projection has been used because of dissimilarity.", "Analysis has been stopped.");
    }

	//Process all cartographic projections from the list one by one
	for (typename TItemsList <Projection <T> *> ::Type::const_iterator i_projections = pl.begin(); i_projections != pl.end(); ++i_projections)
	{
		//Get limits of the cartographic pole latitude and longitude: some projections are defined only in normal position
		total_created_and_analyzed_samples_projection = 0;

		//Print actual projection name to the log
		*output << (*i_projections)->getName() << ": ";

		//List of points using new central meridian redefined in projection file
		Container <Point3DGeographic <T> *> pl_reference_red;

		//Reduce lon using a new central meridian redefined in projection file, if necessary
		if ((*i_projections)->getLon0() != 0.0) redLon(pl_reference, (*i_projections)->getLon0(), pl_reference_red);

		//Set pointer to processed sets: reduced or non reduced
		Container <Point3DGeographic <T> *> * p_pl_reference = ((*i_projections)->getLon0() == 0.0 ? &pl_reference : &pl_reference_red);

		//Create list of possible pole positions
		typename TItemsList <TProjectionPolePosition<T> >::Type proj_pole_positions_list;

		//Get both latp and lonp intervals
		TMinMax <T> latp_interval_heur = (*i_projections)->getLatPInterval();
		TMinMax <T> lonp_interval_heur = (*i_projections)->getLonPInterval();

		//Find intervals of latp, lonp
		if (analysis_parameters.perform_heuristic)
			findLatPLonPIntervals(*p_pl_reference, *i_projections, latp_interval_heur, lonp_interval_heur);

		//Normal aspect
		if (analysis_parameters.analyze_normal_aspect)
			createOptimalLatPLonPPositions(*p_pl_reference, *i_projections, latp_interval_heur, lonp_interval_heur, analysis_parameters, NormalAspect, proj_pole_positions_list, output);

		//Transverse aspect
		if (analysis_parameters.analyze_transverse_aspect)
			createOptimalLatPLonPPositions(*p_pl_reference, *i_projections, latp_interval_heur, lonp_interval_heur, analysis_parameters, TransverseAspect, proj_pole_positions_list, output);

		//Oblique aspect
		if (analysis_parameters.analyze_oblique_aspect)
			createOptimalLatPLonPPositions(*p_pl_reference, *i_projections, latp_interval_heur, lonp_interval_heur, analysis_parameters, ObliqueAspect, proj_pole_positions_list, output);


		//Test, if some singular points has been found
		bool singular_points_found = false;

		//Set pointers to processed non singular sets
		Container <Point3DGeographic <T> *> *p_pl_reference_non_sing = p_pl_reference;
		Container <Node3DCartesian <T> *> *p_nl_test_non_sing = &nl_test;

		//Create temporary containers for non singular points: containers of points
		Container <Node3DCartesian <T> *> nl_test_non_sing;
		Container <Point3DGeographic <T> *> pl_reference_non_sing;

		//Non singular pairs
		typename TDevIndexPairs <T>::Type non_singular_pairs;
		TIndexList non_singular_points;

		//Pointer to meridians and parallels
		typename TMeridiansList <T> ::Type * p_meridians_non_sing = &meridians;
		typename TParallelsList <T> ::Type * p_parallels_non_sing = &parallels;

		//Create temporary meridians and parallels
		typename TMeridiansList <T> ::Type meridians_non_sing;
		typename TParallelsList <T> ::Type parallels_non_sing;

		//Process all found positions
		for (unsigned int i = 0; i < proj_pole_positions_list.size(); i++)
		{
			//Set projection parameters: cartographic pole
			(*i_projections)->setCartPole(proj_pole_positions_list[i].cart_pole);
			(*i_projections)->setLat0(proj_pole_positions_list[i].lat0);

			//Try to remove singular points to prevent a throwing of sample only if  cartographic pole coordinates latp, lonp change
			if ((i == 0) || (i > 0) && (proj_pole_positions_list[i].cart_pole != proj_pole_positions_list[i - 1].cart_pole))
			{
				//Set pointers to old sets (they do not contain singular data)
				p_pl_reference_non_sing = p_pl_reference;
				p_nl_test_non_sing = &nl_test;

				//Remove singular points: empty containers
				nl_test_non_sing.clear(); pl_reference_non_sing.clear(); non_singular_pairs.clear();
				removeSingularPoints(*p_nl_test_non_sing, *p_pl_reference_non_sing, *i_projections, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

				//Singular points found
				if (nl_test.size() != nl_test_non_sing.size())
				{
					//Set flag to true, used for all samples using non-singular sets
					singular_points_found = true;

					//Create copy of meridians / parallels
					meridians_non_sing = meridians;
					parallels_non_sing = parallels;

					//Correct meridians and parallels
					correctMeridiansAndParrallels <T>(meridians_non_sing, parallels_non_sing, non_singular_pairs);

					//Convert non singular pairs to index list: indices will be printed in output
					non_singular_points.clear();
					std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());

					//Set pointers to non-singular meridians and parallels
					p_meridians_non_sing = &meridians_non_sing;
					p_parallels_non_sing = &parallels_non_sing;

					//Set pointers to newly created non-singular sets
					p_nl_test_non_sing = &nl_test_non_sing;
					p_pl_reference_non_sing = &pl_reference_non_sing;
				}

				//Singular points did not find
				else
				{
					//Set flag to false
					singular_points_found = false;

					//Initialize non singular indices: all indices are valid
					non_singular_points.clear();

					for (unsigned int j = 0; j < p_pl_reference->size(); j++) non_singular_points.push_back(j);
				}
			}

			//Compute analysis
			unsigned int created_samples = 0;

			try
			{
				//Analyze normal aspect
				Sample <T> analyzed_sample;

				computeAnalysisForOneSample(*p_nl_test_non_sing, *p_pl_reference_non_sing, meridians, parallels, faces_test, *i_projections, analysis_parameters, analyzed_sample, singular_points_found, created_samples, output);

				//Add result to the list
				if (total_created_and_analyzed_samples_projection) sl.push_back(analyzed_sample);
			}

			//Throw exception
			catch (Exception & error)
			{
				if (analysis_parameters.print_exceptions)
					error.printException();
			}

			//Increment amount of created and thrown samples
			if (created_samples == 0) total_created_or_thrown_samples++;
			else
			{
				total_created_and_analyzed_samples_projection += created_samples;
				total_created_or_thrown_samples += created_samples;
			}

			//Print "." for every 500-th sample to the log
			if (total_created_or_thrown_samples % 500 == 0)
			{
				std::cout.flush();
				std::cout << ".";
			}
		}

		//Print successfully analyzed samples for one cartographic projection
		*output << " [" << total_created_and_analyzed_samples_projection << " iterations]" << std::endl;
	}
}


template <typename T>
void CartAnalysis::computeAnalysisForAllSamplesSim(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
	unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
	//Find mimum using the Simplex method ( Nelder-Mead algorithm )
	const unsigned int m = nl_test.size();

	//Total successfully computed analysis for one cartographic projection
	unsigned int total_created_and_analyzed_samples_projection = 0;

	//Find extreme values
	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	//Compute the approximate map scale: 1:1000 covers 15 x 15 sec
	const T map_scale1 = 1000 / (15.0 / 3600) * std::max(lat_interval.max_val - lat_interval.min_val,
		(lon_interval.max_val - lon_interval.min_val) * cos(0.5 * (lat_interval.max_val + lat_interval.min_val) * M_PI / 180));

	//Create sample for analyzed projection from command line and set flag for this sample
	if (analysis_parameters.analyzed_projections.size() > 0)
	{
		//Analyze all projections specified in command line
		for (typename TItemsList<Projection <T> *> ::Type::iterator i_projections = analysis_parameters.analyzed_projections.begin(); i_projections != analysis_parameters.analyzed_projections.end(); ++i_projections)
		{
			//Get analyzed projection
			Projection <T> *analyzed_proj = *i_projections;

			//List of points using new central meridian redefined in projection file
			Container <Point3DGeographic <T> *> pl_reference_red;

			//Reduce lon using a new central meridian redefined in projection file, if necessary
			if ((*i_projections)->getLon0() != 0.0) redLon(pl_reference, (*i_projections)->getLon0(), pl_reference_red);

			//Set pointer to processed file: reduced or non reduced
			Container <Point3DGeographic <T> *> * p_pl_reference = ((*i_projections)->getLon0() == 0.0 ? &pl_reference : &pl_reference_red);

			//Create temporary containers for non singular points
			Container <Node3DCartesian <T> *> nl_test_non_sing;
			Container <Point3DGeographic <T> *> pl_reference_non_sing;

			typename TDevIndexPairs <T>::Type non_singular_pairs;
			TIndexList non_singular_points;

			//Initialize non singular indices
			for (unsigned int i = 0; i < p_pl_reference->size(); i++) non_singular_points.push_back(i);

			//Set pointer to processed file: with or wihout singular points
			Container <Node3DCartesian <T> *> *p_nl_test = &nl_test;

			//Remove singular points to prevent throwing a sample
			bool singular_points_found = false;
			removeSingularPoints(*p_nl_test, *p_pl_reference, *i_projections, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

			//Some singular points have been found
			if (nl_test.size() != nl_test_non_sing.size())
			{
				//Set pointers to files without singular points
				p_nl_test = &nl_test_non_sing;	p_pl_reference = &pl_reference_non_sing;

				//Set flag to true, used for a sample using non-singular sets
				singular_points_found = true;

				//Correct meridians and parallels
				correctMeridiansAndParrallels <T>(meridians, parallels, non_singular_pairs);

				//Convert non singular pairs to index list: indices will be printed in output
				non_singular_points.clear();
				std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());
			}

			//Create sample
			Sample <T> analyzed_sample;

			//Compute analysis
			try
			{
				(void)computeAnalysisForOneSample(*p_nl_test, *p_pl_reference, meridians, parallels, faces_test, analyzed_proj, analysis_parameters, analyzed_sample, singular_points_found, total_created_and_analyzed_samples_projection, output);
			}

			//Throw error
			catch (Exception & error)
			{
				if (analysis_parameters.print_exceptions) error.printException();
			}

			//Sample with analyzed projection has been successfully created (not thrown by the heuristic)
			if (total_created_and_analyzed_samples_projection > 0) sl[sl.size() - 1].setAnalyzedProjectionSample(true);
		}

		if (total_created_and_analyzed_samples_projection == 0) throw BadDataException("BadDataException: no analyzed projection has been used because of dissimilarity.", "Analysis has been stopped.");
	}

	//Process all cartographic projections from the list one by one
	for (typename TItemsList <Projection <T> *> ::Type::const_iterator i_projections = pl.begin(); i_projections != pl.end(); ++i_projections)
	{
		//Get limits of the cartographic pole latitude and longitude: some projections are defined only in normal position
		total_created_and_analyzed_samples_projection = 0;

		//Get defined radius
		const T R_def = (*i_projections)->getR(), R_est = R_def;

		//Print actual projection name to the log
		std::cout << (*i_projections)->getName() << ": ";
		*output << (*i_projections)->getName() << ": ";

		//Get both latp and lonp intervals: lonp intervals are set to the moved central meridian (further must be reduced)
		TMinMax <T> latp_interval_heur = (*i_projections)->getLatPIntervalH(lat_interval);
		TMinMax <T> lonp_interval_heur = (*i_projections)->getLonPIntervalH(lon_interval);
		TMinMax <T> lat0_interval = (*i_projections)->getLat0Interval();

		const T latp_min = (*i_projections)->getLatPInterval().min_val;
		const T latp_max = (*i_projections)->getLatPInterval().max_val;
		const T lonp_min = (*i_projections)->getLonPInterval().min_val;
		const T lonp_max = (*i_projections)->getLonPInterval().max_val;
		const T lat0_min = lat0_interval.min_val;
		const T lat0_max = lat0_interval.max_val;

		//Create sample
		Sample <T> best_sample;

		//Create matrix of determined variables
		unsigned short n_par = 6;
		if (analysis_parameters.analysis_method == SimplexRotMethod)
			n_par = 7;
		else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
			n_par = 8;

		//Create matrix of determined variables
		Matrix <unsigned int> IX(2 * m, 1);
		Matrix <T> XMIN(1, n_par), X(1, n_par, 1), XMAX(1, n_par, 1), Y(2 * m, 1), W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);

		//Set iteration parameters
		const T eps = 1.0e-10, max_diff = 1.0e-10, k = 2.5, mult = 0.001;
		unsigned int iterations = 0;
		const unsigned int max_iter = 200;

		try
		{
			//Compute initial R value
			unsigned int total_samples_test = 0;
			Sample <T> sample_test;
			TAnalysisParameters <T> analysis_parameters_test(false);
			analysis_parameters_test.analysis_type.a_helt = true;

			//Radius mutliplier
			const T ks = 0.01;

			//Store projection properties before analysis
			const T R = (*i_projections)->getR();
			Point3DGeographic <T> cart_pole = (*i_projections)->getCartPole();
			const T lat0 = (*i_projections)->getLat0();
			const T lat1 = (*i_projections)->getLat1();
			const T lat2 = (*i_projections)->getLat2();
			const T lon0 = (*i_projections)->getLon0();
			const T dx = (*i_projections)->getDx();
			const T dy = (*i_projections)->getDy();
			const T c = (*i_projections)->getC();

			//Conditions
			const bool normal_aspect_enabled = (analysis_parameters.analyze_normal_aspect) && ((*i_projections)->getCartPole().getLat() == MAX_LAT);
			const bool transverse_aspect_enabled = (analysis_parameters.analyze_transverse_aspect) && ((*i_projections)->getCartPole().getLat() == 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);
			const bool oblique_aspect_enabled = (analysis_parameters.analyze_oblique_aspect) && ((*i_projections)->getCartPole().getLat() != 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);

			//Initialize matrices
			Matrix <T>  X0MIN(1, 5), X0MAX(1, 5), X0(1, 5);

			//Common for all apects
			X0MIN(0, 2) = lat0_min;     X0MAX(0, 2) = lat0_max;

			//Conic projections: c = lat2
			if (typeid(*(*i_projections)) == typeid(ProjectionConic<T>))
			{
				X0MIN(0, 4) = lat0_min;          X0MAX(0, 4) = lat0_max;;
			}

			//Otherwise: c is unspecified constant
			else
			{
				X0MIN(0, 4) = 0.0;          X0MAX(0, 4) = 1.0e3;
			}

			//Proces available aspects
			TProjectionAspect aspect = NormalAspect;
			for (unsigned int i = 0; i < 3; i++)
			{
				//At least one test is set
				if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
				{
					//Change analyzed aspect
					if (i == 1) aspect = TransverseAspect;
					else if (i == 2) aspect = ObliqueAspect;

					//Get current limits
					const T latp_init = 0.5 * (latp_interval_heur.min_val + latp_interval_heur.max_val);
					const T lonp_init = 0.5 * (lonp_interval_heur.min_val + lonp_interval_heur.max_val);

					//Create new meta-pole
					const Point3DGeographic <T> cart_pole_0(latp_init, lonp_init);

					//Set new position of the meta-pole for the transverse and oblique aspect
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole_0);

					//Use the cartometric analysis to set the initial value of the Earth radius: use the similarity transformation
					CartAnalysis::computeAnalysisForOneSample(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters_test, sample_test, false, total_samples_test, output);

					//Get initialize Earth radius using the similarity transformation
					const T R0 = (*i_projections)->getR() / sample_test.getScaleHelT();

					//Restore meta-pole position
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole);

					//Normal aspect
					if ((i == 0) && (normal_aspect_enabled))
					{

						//Set intervals for M7S method
						X0MIN(0, 0) = 90;           X0MAX(0, 0) = 90;
						X0MIN(0, 1) = 0;            X0MAX(0, 1) = 0;
						X0MIN(0, 3) = -MAX_LON;     X0MAX(0, 3) = MAX_LON;
					}

					else if ((i == 1) && (transverse_aspect_enabled))
					{
						//Set intervals for M7S method
						X0MIN(0, 0) = 0;		X0MAX(0, 0) = 0;
						X0MIN(0, 1) = lonp_min;		X0MAX(0, 1) = lonp_max;
						X0MIN(0, 3) = 0;		X0MAX(0, 3) = 0;
					}

					else if ((i == 2) && (oblique_aspect_enabled))
					{
						X0MIN(0, 0) = latp_min;		X0MAX(0, 0) = latp_max;
						X0MIN(0, 1) = lonp_min;		X0MAX(0, 1) = lonp_max;
						X0MIN(0, 3) = lonp_min;		X0MAX(0, 3) = lonp_max;
					}


					//X0MIN(0, 0) = 43;	X0MAX(0, 0) = 50;
					//X0MIN(0, 1) = 0;	X0MAX(0, 1) = 0;

					//X0MIN(0, 2) = 0;	X0MAX(0, 2) = 0;
					//X0MIN(0, 3) = 100;	X0MAX(0, 3) = 100;
					//X0MIN(0, 4) = 0;	X0MAX(0, 4) = 1000;

					//Variables
					unsigned int res_evaluation = 0;
					T x_mass_reference = 0.0, y_mass_reference = 0.0, min_cost = 0, alpha = 0, R_est = 1.0, scale = 1.0, q1 = 1, q2 = 1, k = 2.5;
					TMEstimatorsWeightFunction me_function = HuberFunction;

					//Matrices
					Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);
					Matrix <unsigned int> IX(2 * m, 1);

					Container <Node3DCartesianProjected <T> *> nl_projected;

					//Set initial solution and intervals depending ion the selected method
					if (analysis_parameters.analysis_method == SimplexMethod)
					{
						XMIN.replace(X0MIN, 0, 1); XMAX.replace(X0MAX, 0, 1);
						XMIN(0, 0) = ks * R0; XMAX(0, 0) = 1.0 / ks * R0;

						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters, aspect, best_sample,
							total_created_and_analyzed_samples_projection, res_evaluation, HuberFunction, k, IX, output), W, X, Y, V, XMIN, XMAX, iterations, mult * eps, max_iter, output);
					}

					else if (analysis_parameters.analysis_method == SimplexRotMethod)
					{
						XMIN.replace(X0MIN, 0, 1); XMAX.replace(X0MAX, 0, 1);
						XMIN(0, 0) = ks * R0; XMAX(0, 0) = 1.0 / ks * R0;
						XMIN(0, 6) = -MAX_LON; XMAX(0, 6) = MAX_LON;

						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters, aspect, best_sample,
							total_created_and_analyzed_samples_projection, res_evaluation, XMIN, XMAX, output), W, X, Y, V, XMIN, XMAX, iterations, mult * eps, max_iter, output);
					}

					else if (analysis_parameters.analysis_method == SimplexRot2Method)
					{
						//X0.print();
						//X0MIN.print();
						//X0MAX.print();
						/*
						q1 = 1.1226509817934028e-05;
						q2 = -4.4778786203080303e-07;
						R_est = 0.071682085811494392;
						*/
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, R_est, q1, q2, analysis_parameters, aspect, best_sample,
								total_created_and_analyzed_samples_projection, res_evaluation, HuberFunction, k, IX, output), W, X0, Y, V, X0MIN, X0MAX, iterations, mult * eps, max_iter, output);

						//X0.print();
						X(0, 0) = R_est;  X.replace(X0, 0, 1);

						scale = sqrt(q1 * q1 + q2 * q2); alpha = atan2(q2, q1) * 180 / M_PI;
					}

					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
					{
						XMIN.replace(X0MIN(0, 0, 0, 3), 1, 0); XMAX.replace(X0MAX(0, 0, 0, 3), 0, 1);
						XMIN(0, 0) = ks * R0;	XMAX(0, 0) = 1.0 / ks * X(0, 0);
						XMIN(0, 5) = -1.0e03;	XMAX(0, 5) = 1.0e3;
						XMIN(0, 6) = -1.0e03;	XMAX(0, 6) = 1.0e3;
						XMIN(0, 7) = 0;		XMAX(0, 7) = 1.0e8;

						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters, aspect, best_sample,
							total_created_and_analyzed_samples_projection, res_evaluation, output), W, X, Y, V, XMIN, XMAX, iterations, mult * eps, max_iter, output);
					}

					const T map_scale2 = (analysis_parameters.analysis_method != NonLinearLeastSquaresRot2Method ? R_def / X(0, 0) : R_def / R_est);

					//Matrix <T> WD = diag(W); WD.print(output);

					//Compute match ratio
					//TIndexList matched_points;
					//T match_ratio = getMatchRatioTissotIndicatrix(nl_test, pl_reference, trans(X), W, *i_projections, matched_points, CollectOn, map_scale, 30.0);
					//best_sample.setHomotheticTransformationPercMatch((unsigned int)match_ratio);
					//best_sample.setHomotheticTransformationMatchedPointsIndices ( matched_points );

					//Set residuals
					best_sample.setHomotheticTransformationRatio(min_cost);

					best_sample.setR(X(0, 0));
					best_sample.setLatP(X(0, 1));
					best_sample.setLonP(X(0, 2));
					best_sample.setLat0(X(0, 3));
					best_sample.setLat1(X(0, 3));
					best_sample.setLat2(X(0, 5));
					best_sample.setLon0(X(0, 4));
					best_sample.setC(X(0, 5));

					if (analysis_parameters.analysis_method == SimplexRotMethod)
						best_sample.setAlpha(X(0, 6));
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						best_sample.setAlpha(alpha);

					best_sample.setScaleHomT(map_scale2);
					best_sample.setIterations(iterations);
					best_sample.setResidualEval(res_evaluation);
					best_sample.setProj(*i_projections);

					//Add to the list
					if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
						sl.push_back(best_sample);

					//Restore projection properties after analysis
					(*i_projections)->setR(R);
					(*i_projections)->setCartPole(cart_pole);
					(*i_projections)->setLat0(lat0);
					(*i_projections)->setLat1(lat1);
					(*i_projections)->setLat2(lat2);
					(*i_projections)->setLon0(lon0);
					(*i_projections)->setDx(dx);
					(*i_projections)->setDy(dy);
					(*i_projections)->setC(c);

					//Print successfully analyzed samples for one cartographic projection
					std::cout << " [" << iterations << " it, f = " << min_cost << "] ";
					*output << " [" << iterations << " it, f = " << min_cost << "] ";
				}
			}

			std::cout << '\n';
			*output << '\n';
		}

		//Throw exception, analysis have not been computed
		catch (Exception & error)
		{
			if (analysis_parameters.print_exceptions)
				error.printException();
		}

		//Assign the amount of created samples
		total_created_or_thrown_samples += total_created_and_analyzed_samples_projection;
	}

}


template <typename T>
void CartAnalysis::computeAnalysisForAllSamplesDE(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
	unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
	//Find mimum using the differential evolution algorithm
	const unsigned int m = nl_test.size();

	//Total successfully computed analysis for one cartographic projection
	unsigned int total_created_and_analyzed_samples_projection = 0;

	//Find extreme values
	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	//Compute the approximate map scale: 1:1000 covers 15 x 15 sec
	const T map_scale1 = 1000 / (20.0 / 3600) * std::max(lat_interval.max_val - lat_interval.min_val,
		(lon_interval.max_val - lon_interval.min_val) * cos(0.5 * (lat_interval.max_val + lat_interval.min_val) * M_PI / 180));

	//Create sample for analyzed projection from command line and set flag for this sample
	if (analysis_parameters.analyzed_projections.size() > 0)
	{
		//Analyze all projections specified in command line
		for (typename TItemsList<Projection <T> *> ::Type::iterator i_projections = analysis_parameters.analyzed_projections.begin(); i_projections != analysis_parameters.analyzed_projections.end(); ++i_projections)
		{
			//Get analyzed projection
			Projection <T> *analyzed_proj = *i_projections;

			//List of points using new central meridian redefined in projection file
			Container <Point3DGeographic <T> *> pl_reference_red;

			//Reduce lon using a new central meridian redefined in projection file, if necessary
			if ((*i_projections)->getLon0() != 0.0) redLon(pl_reference, (*i_projections)->getLon0(), pl_reference_red);

			//Set pointer to processed file: reduced or non reduced
			Container <Point3DGeographic <T> *> * p_pl_reference = ((*i_projections)->getLon0() == 0.0 ? &pl_reference : &pl_reference_red);

			//Create temporary containers for non singular points
			Container <Node3DCartesian <T> *> nl_test_non_sing;
			Container <Point3DGeographic <T> *> pl_reference_non_sing;

			typename TDevIndexPairs <T>::Type non_singular_pairs;
			TIndexList non_singular_points;

			//Initialize non singular indices
			for (unsigned int i = 0; i < p_pl_reference->size(); i++) non_singular_points.push_back(i);

			//Set pointer to processed file: with or wihout singular points
			Container <Node3DCartesian <T> *> *p_nl_test = &nl_test;

			//Remove singular points to prevent throwing a sample
			bool singular_points_found = false;
			removeSingularPoints(*p_nl_test, *p_pl_reference, *i_projections, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

			//Some singular points have been found
			if (nl_test.size() != nl_test_non_sing.size())
			{
				//Set pointers to files without singular points
				p_nl_test = &nl_test_non_sing;	p_pl_reference = &pl_reference_non_sing;

				//Set flag to true, used for a sample using non-singular sets
				singular_points_found = true;

				//Correct meridians and parallels
				correctMeridiansAndParrallels <T>(meridians, parallels, non_singular_pairs);

				//Convert non singular pairs to index list: indices will be printed in output
				non_singular_points.clear();
				std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());
			}

			//Create sample
			Sample <T> analyzed_sample;

			//Compute analysis
			try
			{
				(void)computeAnalysisForOneSample(*p_nl_test, *p_pl_reference, meridians, parallels, faces_test, analyzed_proj, analysis_parameters, analyzed_sample, singular_points_found, total_created_and_analyzed_samples_projection, output);
			}

			//Throw exception
			catch (Exception & error)
			{
				if (analysis_parameters.print_exceptions) error.printException();
			}

			//Sample with analyzed projection has been successfully created (not thrown by the heuristic)
			if (total_created_and_analyzed_samples_projection > 0) sl[sl.size() - 1].setAnalyzedProjectionSample(true);
		}

		if (total_created_and_analyzed_samples_projection == 0) throw BadDataException("BadDataException: no analyzed projection has been used because of dissimilarity.", "Analysis has been stopped.");
	}

	//Process all cartographic projections from the list one by one
	for (typename TItemsList <Projection <T> *> ::Type::const_iterator i_projections = pl.begin(); i_projections != pl.end(); ++i_projections)
	{
		//Get limits of the cartographic pole latitude and longitude: some projections are defined only in normal position
		total_created_and_analyzed_samples_projection = 0;

		//Get defined radius
		const T R_def = (*i_projections)->getR(), R_est = R_def;

		//Print actual projection name to the log
		std::cout << (*i_projections)->getName() << ": ";
		*output << (*i_projections)->getName() << ": ";

		//Get both latp and lonp intervals: lonp intervals are set to the moved central meridian (further must be reduced)
		TMinMax <T> latp_interval_heur = (*i_projections)->getLatPIntervalH(lat_interval);
		TMinMax <T> lonp_interval_heur = (*i_projections)->getLonPIntervalH(lon_interval);
		TMinMax <T> lat0_interval = (*i_projections)->getLat0Interval();

		const T latp_min = (*i_projections)->getLatPInterval().min_val;
		const T latp_max = (*i_projections)->getLatPInterval().max_val;
		const T lonp_min = (*i_projections)->getLonPInterval().min_val;
		const T lonp_max = (*i_projections)->getLonPInterval().max_val;
		const T lat0_min = lat0_interval.min_val;
		const T lat0_max = lat0_interval.max_val;

		//Create sample
		Sample <T> best_sample;

		//Create matrix of determined variables
		unsigned short n_par = 6;
		if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
			n_par = 7;
		else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
			n_par = 8;

		//Create matrix of determined variables
		Matrix <unsigned int> IX(2 * m, 1);

		Matrix <T> XMIN(1, n_par), XMAX(1, n_par), XMIN2(1, n_par), XMAX2(1, n_par), X(1, n_par), XAVER(1, n_par), Y(2 * m, 1), W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);

		//Parameters of the genetic algorithm
		const T eps = 1.0e-10, CR = 0.8, k = 2.5;
		const unsigned int population = 2*n_par * XMIN.cols(), max_gen = 150;
		unsigned int iterations = 0;

		Matrix <T> F(1, 1); F(0, 0) = 0.5;

		try
		{
			//Compute initial R value
			unsigned int total_samples_test = 0;
			Sample <T> sample_test;
			TAnalysisParameters <T> analysis_parameters_test(false);
			analysis_parameters_test.analysis_type.a_helt = true;

			//Radius mutliplier
			const T ks = 0.01;

			//Store projection properties before analysis
			const T R = (*i_projections)->getR();
			Point3DGeographic <T> cart_pole = (*i_projections)->getCartPole();
			const T lat0 = (*i_projections)->getLat0();
			const T lat1 = (*i_projections)->getLat0();
			const T lat2 = (*i_projections)->getLat0();
			const T lon0 = (*i_projections)->getLon0();
			const T dx = (*i_projections)->getDx();
			const T dy = (*i_projections)->getDy();
			const T c = (*i_projections)->getC();

			//Conditions
			const bool normal_aspect_enabled = (analysis_parameters.analyze_normal_aspect) && ((*i_projections)->getCartPole().getLat() == MAX_LAT);
			const bool transverse_aspect_enabled = (analysis_parameters.analyze_transverse_aspect) && ((*i_projections)->getCartPole().getLat() == 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);
			const bool oblique_aspect_enabled = (analysis_parameters.analyze_oblique_aspect) && ((*i_projections)->getCartPole().getLat() != 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);

			//Initialize matrices
			Matrix <T>  X0MIN(1, 5), X0MAX(1, 5), X0(1, 5), XAVER0(1, 5);

			//Common for all apects
			X0MIN(0, 2) = lat0_min;     X0MAX(0, 2) = lat0_max;

			//Conic projections: c = lat2
			if (typeid(*(*i_projections)) == typeid(ProjectionConic<T>))
			{
				X0MIN(0, 4) = lat0_min;          X0MAX(0, 4) = lat0_max;;
			}

			//Otherwise: c is unspecified constant
			else
			{
				X0MIN(0, 4) = 0.0;          X0MAX(0, 4) = 1.0e3;
			}


			//Proces available aspects
			TProjectionAspect aspect = NormalAspect;
			for (unsigned int i = 0; i < 3; i++)
			{
				//At least one test is set
				if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
				{
					//Change analyzed aspect
					if (i == 1) aspect = TransverseAspect;
					else if (i == 2) aspect = ObliqueAspect;

					//Get current limits
					const T latp_init = 0.5 * (latp_interval_heur.min_val + latp_interval_heur.max_val);
					const T lonp_init = 0.5 * (lonp_interval_heur.min_val + lonp_interval_heur.max_val);

					//Create new meta-pole
					const Point3DGeographic <T> cart_pole_0(latp_init, lonp_init);

					//Set new position of the meta-pole for the transverse and oblique aspect
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole_0);

					//Use the cartometric analysis to set the initial value of the Earth radius: use the similarity transformation
					CartAnalysis::computeAnalysisForOneSample(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters_test, sample_test, false, total_samples_test, output);

					//Get initialize Earth radius using the similarity transformation
					const T R0 = (*i_projections)->getR() / sample_test.getScaleHelT();

					//Restore meta-pole position
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole);

					//Normal aspect
					if ((i == 0) && (normal_aspect_enabled))
					{

						//Set intervals for M7S method
						X0MIN(0, 0) = 90;           X0MAX(0, 0) = 90;
						X0MIN(0, 1) = 0;            X0MAX(0, 1) = 0;
						X0MIN(0, 3) = -MAX_LON;     X0MAX(0, 3) = MAX_LON;
					}

					else if ((i == 1) && (transverse_aspect_enabled))
					{
						//Set intervals for M7S method
						X0MIN(0, 0) = 0;            X0MAX(0, 0) = 0;
						X0MIN(0, 1) = lonp_min;     X0MAX(0, 1) = lonp_max;
						X0MIN(0, 3) = 0;		X0MAX(0, 3) = 0;
					}

					else if ((i == 2) && (oblique_aspect_enabled))
					{
						X0MIN(0, 0) = latp_min;	X0MAX(0, 0) = latp_max;
						X0MIN(0, 1) = lonp_min;	X0MAX(0, 1) = lonp_max;
						X0MIN(0, 3) = 0;		X0MAX(0, 3) = 0;
					}

					//Variables
					unsigned int res_evaluation = 0;
					T x_mass_reference = 0.0, y_mass_reference = 0.0, min_cost = 0, alpha = 0, R_est = 1.0, scale = 1.0, q1 = 1, q2 = 1, k = 2.5, res_aver, res_max;
					TMEstimatorsWeightFunction me_function = HuberFunction;

					//Matrices
					Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);
					Matrix <unsigned int> IX(2 * m, 1);

					Container <Node3DCartesianProjected <T> *> nl_projected;

					//Set initial solution and intervals depending ion the selected method
					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
					{
						XMIN.replace(X0MIN, 0, 1); XMAX.replace(X0MAX, 0, 1);
						XMIN(0, 0) = ks * R0; XMAX(0, 0) = 1.0 / ks * R0;

						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters, aspect, best_sample,
							total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, output), population, eps, max_gen, F, CR, DERandBest1Strategy, MFDE, W, X, Y, V, XMIN, XMAX, XAVER, res_aver, res_max, iterations);
					}

					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
					{
						XMIN.replace(X0MIN, 0, 1); XMAX.replace(X0MAX, 0, 1);
						XMIN(0, 0) = ks * R0; XMAX(0, 0) = 1.0 / ks * R0;
						XMIN(0, 6) = -MAX_LON; XMAX(0, 6) = MAX_LON;

						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters, aspect, best_sample,
							total_created_and_analyzed_samples_projection, res_evaluation, XMIN, XMAX, XAVER, output), population, eps, max_gen, F, CR, DERand1Strategy, MFDE, W, X, Y, V, XMIN, XMAX, XAVER, res_aver, res_max, iterations);
					}

					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
					{
						//X0.print();
						//X0MIN.print();
						//X0MAX.print();

						unsigned int eff = 0;
						for (int index = 0; index < 100; index++)
						{
							min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, R_est, q1, q2, analysis_parameters, aspect, best_sample,
								total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, output), population, eps, max_gen, F, CR, DERand1Strategy, MFDE, W, X0, Y, V, X0MIN, X0MAX, XAVER0, res_aver, res_max, iterations);
						
							if (min_cost < 1.0e-4) eff++;
						}

						std::cout << "eff = " << eff;

						X(0, 0) = R_est;  X.replace(X0, 0, 1);

						scale = sqrt(q1 * q1 + q2 * q2); alpha = atan2(q2, q1) * 180 / M_PI;
					}

					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
					{
						XMIN.replace(X0MIN(0, 0, 0, 3), 1, 0); XMAX.replace(X0MAX(0, 0, 0, 3), 0, 1);
						XMIN(0, 0) = ks * R0;	XMAX(0, 0) = 1.0 / ks * X(0, 0);
						XMIN(0, 5) = -1.0e03;	XMAX(0, 5) = 1.0e3;
						XMIN(0, 6) = -1.0e03;	XMAX(0, 6) = 1.0e3;
						XMIN(0, 7) = 0;		XMAX(0, 7) = 1.0e8;

						//min_cost = 					}
					}

					//Matrix <T> WD = diag(W); WD.print(output);

					const T map_scale2 = (analysis_parameters.analysis_method != NonLinearLeastSquaresRot2Method ? R_def / X(0, 0) : R_def / R_est);

					//Compute match ratio
					//TIndexList matched_points;
					//T match_ratio = getMatchRatioTissotIndicatrix(nl_test, pl_reference, trans(X), W, *i_projections, matched_points, CollectOn, map_scale, 30.0);
					//best_sample.setHomotheticTransformationPercMatch((unsigned int)match_ratio);
					//best_sample.setHomotheticTransformationMatchedPointsIndices ( matched_points );

					//Set residuals
					best_sample.setHomotheticTransformationRatio(min_cost);

					best_sample.setR(X(0, 0));
					best_sample.setLatP(X(0, 1));
					best_sample.setLonP(X(0, 2));
					best_sample.setLat0(X(0, 3));
					best_sample.setLat1(X(0, 3));
					best_sample.setLat2(X(0, 5));
					best_sample.setLon0(X(0, 4));
					best_sample.setC(X(0, 5));

					if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						best_sample.setAlpha(X(0, 6));
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						best_sample.setAlpha(alpha);

					best_sample.setScaleHomT(map_scale2);
					best_sample.setIterations(iterations);
					best_sample.setResidualEval(res_evaluation);
					best_sample.setProj(*i_projections);

					//Add to the list
					if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
						sl.push_back(best_sample);

					//Restore projection properties after analysis
					(*i_projections)->setR(R);
					(*i_projections)->setCartPole(cart_pole);
					(*i_projections)->setLat0(lat0);
					(*i_projections)->setLat1(lat1);
					(*i_projections)->setLat2(lat2);
					(*i_projections)->setLon0(lon0);
					(*i_projections)->setDx(dx);
					(*i_projections)->setDy(dy);
					(*i_projections)->setC(c);

					//Print successfully analyzed samples for one cartographic projection
					std::cout << " [" << iterations << " it, f = " << min_cost << "] ";
					*output << " [" << iterations << " it, f = " << min_cost << "] ";
				}
			}

			std::cout << '\n';
			*output << '\n';
		}

		//Throw exception
		catch (Exception & error)
		{
			if (analysis_parameters.print_exceptions)
				error.printException();
		}

		//Assign the amount of created samples
		total_created_or_thrown_samples += total_created_and_analyzed_samples_projection;
	}
}


template <typename T>
void CartAnalysis::computeAnalysisForAllSamplesNLS(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
	unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
	//Find mimum using the Minimum Least Squares
	const unsigned int m = nl_test.size();

	//Total successfully computed analysis for one cartographic projection
	unsigned int total_created_and_analyzed_samples_projection = 0;

	//Find extreme values
	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	//Compute the approximate map scale: 1:1000 covers 15 x 15 sec
	const T map_scale1 = 1000 / (15.0 / 3600) * std::max(lat_interval.max_val - lat_interval.min_val,
		(lon_interval.max_val - lon_interval.min_val) * cos(0.5 * (lat_interval.max_val + lat_interval.min_val) * M_PI / 180));

	//Create sample for analyzed projection from command line and set flag for this sample
	if (analysis_parameters.analyzed_projections.size() > 0)
	{
		//Analyze all projections specified in command line
		for (typename TItemsList<Projection <T> *> ::Type::iterator i_projections = analysis_parameters.analyzed_projections.begin(); i_projections != analysis_parameters.analyzed_projections.end(); ++i_projections)
		{
			//Get analyzed projection
			Projection <T> *analyzed_proj = *i_projections;

			//List of points using new central meridian redefined in projection file
			Container <Point3DGeographic <T> *> pl_reference_red;

			//Reduce lon using a new central meridian redefined in projection file, if necessary
			if ((*i_projections)->getLon0() != 0.0) redLon(pl_reference, (*i_projections)->getLon0(), pl_reference_red);

			//Set pointer to processed file: reduced or non reduced
			Container <Point3DGeographic <T> *> * p_pl_reference = ((*i_projections)->getLon0() == 0.0 ? &pl_reference : &pl_reference_red);

			//Create temporary containers for non singular points
			Container <Node3DCartesian <T> *> nl_test_non_sing;
			Container <Point3DGeographic <T> *> pl_reference_non_sing;

			typename TDevIndexPairs <T>::Type non_singular_pairs;
			TIndexList non_singular_points;

			//Initialize non singular indices
			for (unsigned int i = 0; i < p_pl_reference->size(); i++) non_singular_points.push_back(i);

			//Set pointer to processed file: with or wihout singular points
			Container <Node3DCartesian <T> *> *p_nl_test = &nl_test;

			//Remove singular points to prevent throwing a sample
			bool singular_points_found = false;
			removeSingularPoints(*p_nl_test, *p_pl_reference, *i_projections, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

			//Some singular points have been found
			if (nl_test.size() != nl_test_non_sing.size())
			{
				//Set pointers to files without singular points
				p_nl_test = &nl_test_non_sing;	p_pl_reference = &pl_reference_non_sing;

				//Set flag to true, used for a sample using non-singular sets
				singular_points_found = true;

				//Correct meridians and parallels
				correctMeridiansAndParrallels <T>(meridians, parallels, non_singular_pairs);

				//Convert non singular pairs to index list: indices will be printed in output
				non_singular_points.clear();
				std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());
			}

			//Create sample
			Sample <T> analyzed_sample;

			//Compute analysis
			try
			{
				(void)computeAnalysisForOneSample(*p_nl_test, *p_pl_reference, meridians, parallels, faces_test, analyzed_proj, analysis_parameters, analyzed_sample, singular_points_found, total_created_and_analyzed_samples_projection, output);
			}

			//Throw error
			catch (Exception & error)
			{
				if (analysis_parameters.print_exceptions) error.printException();
			}

			//Sample with analyzed projection has been successfully created (not thrown by the heuristic)
			if (total_created_and_analyzed_samples_projection > 0) sl[sl.size() - 1].setAnalyzedProjectionSample(true);
		}

		if (total_created_and_analyzed_samples_projection == 0) throw BadDataException("BadDataException: no analyzed projection has been used because of dissimilarity.", "Analysis has been stopped.");
	}

	//Process all cartographic projections from the list one by one
	for (typename TItemsList <Projection <T> *> ::Type::const_iterator i_projections = pl.begin(); i_projections != pl.end(); ++i_projections)
	{
		//Get limits of the cartographic pole latitude and longitude: some projections are defined only in normal position
		total_created_and_analyzed_samples_projection = 0;

		//Get defined radius
		T R_def = (*i_projections)->getR(), R_est = R_def;

		//Print actual projection name to the log
		std::cout << (*i_projections)->getName() << ": ";
		*output << (*i_projections)->getName() << ": ";

		//Get both latp and lonp intervals: lonp intervals are set to the moved central meridian (further must be reduced)
		TMinMax <T> latp_interval_heur = (*i_projections)->getLatPIntervalH(lat_interval);
		TMinMax <T> lonp_interval_heur = (*i_projections)->getLonPIntervalH(lon_interval);
		TMinMax <T> lat0_interval = (*i_projections)->getLat0Interval();

		const T latp_min = (*i_projections)->getLatPInterval().min_val;
		const T latp_max = (*i_projections)->getLatPInterval().max_val;
		const T lonp_min = (*i_projections)->getLonPInterval().min_val;
		const T lonp_max = (*i_projections)->getLonPInterval().max_val;
		const T lat0_min = lat0_interval.min_val;
		const T lat0_max = lat0_interval.max_val;

		//Create sample
		Sample <T> best_sample;

		//Create matrix of determined variables
		unsigned short n_par = 6;
		if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
			n_par = 7;
		else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
			n_par = 8;

		Matrix <unsigned int> IX = ones(2 * m, 1, 1.0);
		Matrix <T> X(n_par, 1), XMIN(n_par, 1), XMAX(n_par, 1), Y(2 * m, 1);

		//Set iteration parameters
		const T eps = 1.0e-10, max_diff = 1.0e-10, k = 2.5, ks = 0.01;
		T alpha_mult = 0.0001, nu = 0.0001;
		unsigned short iterations = 0;
		const unsigned int max_iter = 200;
		const T nu1 = 0.25, nu2 = 0.75, nu3 = 0.0001, gamma1 = 0.5, gamma2 = 2.0, lambda_min = 1.0e-6, lambda_max = 1.0e6;

		try
		{
			unsigned int total_samples_test = 0;
			Sample <T> sample_test;
			TAnalysisParameters <T> analysis_parameters_test(false);
			analysis_parameters_test.analysis_type.a_helt = true;

			//Store projection properties before analysis
			const T R = (*i_projections)->getR();
			const Point3DGeographic <T> cart_pole = (*i_projections)->getCartPole();
			const T lat0 = (*i_projections)->getLat0();
			const T lat1 = (*i_projections)->getLat1();
			const T lat2 = (*i_projections)->getLat2();
			const T lon0 = (*i_projections)->getLon0();
			const T dx = (*i_projections)->getDx();
			const T dy = (*i_projections)->getDy();
			const T c = (*i_projections)->getC();
			//const T alpha = ( *i_projections )->getAlpha();

			//Assign created samples amount
			const unsigned int total_created_and_analyzed_samples_projection_before = total_created_and_analyzed_samples_projection;

			//Compute mean of longitude to initialize lon0
			T lon_mean = 0;

			for (unsigned int i = 0; i < m; i++) lon_mean += pl_reference[i]->getLon();
			lon_mean /= m;

			//Conditions
			const bool normal_aspect_enabled = (analysis_parameters.analyze_normal_aspect) && ((*i_projections)->getCartPole().getLat() == MAX_LAT);
			const bool transverse_aspect_enabled = (analysis_parameters.analyze_transverse_aspect) && ((*i_projections)->getCartPole().getLat() == 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);
			const bool oblique_aspect_enabled = (analysis_parameters.analyze_oblique_aspect) && ((*i_projections)->getCartPole().getLat() != 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);

			//Compute lat0
			T lat0_init = 0.5 * (lat0_interval.min_val + lat0_interval.max_val);

			//Initialize matrices
			Matrix <T> X0(5, 1), X0MIN(5, 1), X0MAX(5, 1), XX(5, 1);

			//Common for all apects
			X0(2, 0) = lat0_init + 10;
			
			X0MIN(2, 0) = lat0_min;     X0MAX(2, 0) = lat0_max;
			
			//Conic projections: c = lat2
			if (typeid(*(*i_projections)) == typeid(ProjectionConic<T>))
			{
				X0MIN(4, 0) = 0;    X0MAX(4, 0) = lat0_max;
				X0(4, 0) = lat0_init + 20;
			}
			
			//Otherwise: c is unspecified constant
			else
			{
				X0MIN(4, 0) = 0.0;          X0MAX(4, 0) = 1.0e3;
				X0(4, 0) = 1;
			}
			
			//Proces available aspects
			TProjectionAspect aspect = NormalAspect;
			for (unsigned int i = 0; i < 3; i++)
			{
				//At least one test is set
				if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
				{
					//Change analyzed aspect
					if (i == 1) aspect = TransverseAspect;
					else if (i == 2) aspect = ObliqueAspect;

					//Get current limits
					const T latp_init = 0.5 * (latp_interval_heur.min_val + latp_interval_heur.max_val);
					const T lonp_init = 0.5 * (lonp_interval_heur.min_val + lonp_interval_heur.max_val);

					//Create new meta-pole
					const Point3DGeographic <T> cart_pole_0(latp_init, lonp_init);

					//Set new position of the meta-pole for the transverse and oblique aspect
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole_0);

					//Use the cartometric analysis to set the initial value of the Earth radius: use the similarity transformation
					CartAnalysis::computeAnalysisForOneSample(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters_test, sample_test, false, total_samples_test, output);

					//Get initialize Earth radius using the similarity transformation
					T R0 = (*i_projections)->getR() / sample_test.getScaleHelT();

					//Restore meta-pole position
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole);

					//Normal aspect
					if ((i == 0) && (normal_aspect_enabled))
					{
						//Initial solution for M7S method
						X0(0, 0) = MAX_LAT;
						X0(1, 0) = 0;
						X0(3, 0) = 10.0;

						//Set intervals for M7S method
						X0MIN(0, 0) = 90;           X0MAX(0, 0) = 90;
						X0MIN(1, 0) = 0;            X0MAX(1, 0) = 0;
						X0MIN(3, 0) = -MAX_LON;     X0MAX(3, 0) = MAX_LON;
					}

					else if ((i == 1) && (transverse_aspect_enabled))
					{
						//Initial solution for M7S method
						X0(0, 0) = 0;
						X0(1, 0) = lonp_init + 10;
						X0(3, 0) = 0;

						//Set intervals for M7S method
						X0MIN(0, 0) = 0;            X0MAX(0, 0) = 0;
						X0MIN(1, 0) = lonp_min;     X0MAX(1, 0) = lonp_max;
						X0MIN(3, 0) = 0;	    X0MAX(3, 0) = 0;
					}

					else if ((i == 2) && (oblique_aspect_enabled))
					{
						X0(0, 0) = latp_init + 20;
						X0(1, 0) = lonp_init + 10;
						X0(3, 0) = 10;
						/*
						X0(0, 0) = 90; X0(1, 0) = 17.6;  //Erich3
						X0(2, 0) = 51.2;
						R0 = 2453;
						*/
						X0MIN(0, 0) = latp_min;	X0MAX(0, 0) = latp_max;
						X0MIN(1, 0) = lonp_min;	X0MAX(1, 0) = lonp_max;
						X0MIN(3, 0) = 0;	X0MAX(3, 0) = 0; 
						X0MIN(3, 0) = - MAX_LON;	X0MAX(3, 0) = MAX_LON;
					}
					
					//Variables
					unsigned int res_evaluation = 0, jac_evaluation = 0;
					bool enable_additional_lon0_analysis = false;
					T x_mass_reference = 0.0, y_mass_reference = 0.0, min_cost = 0, alpha = 0, R_est = 1.0, scale = 1.0, q1 = 1, q2 = 1, k = 2.5;
					TMEstimatorsWeightFunction me_function = HuberFunction;

					//Matrices
					Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);
					Matrix <unsigned int> IX(2 * m, 1);

					Container <Node3DCartesianProjected <T> *> nl_projected;

					//Set initial solution and intervals depending on the selected method
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					{
						X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
						/*
						X(0, 0) = 0.1207;
						X(1, 0) = 90;
						X(2, 0) = 10;
						X(3, 0) = 2.9;
						X(4, 0) = -91.2;
						X(5, 0) = 0;
						X.print();
						*/
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test, pl_reference, (*i_projections), aspect, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
							analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, HuberFunction, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					{
						X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
						X(6, 0) = 30;  XMIN(6, 0) = -MAX_LON; XMAX(6, 0) = MAX_LON;

						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, (*i_projections), aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, *i_projections,
							x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					{
						XX = X0;

						//min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, (*i_projections), aspect, R_est, q1, q2, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
						//	R_est, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, XX, Y, V, X0MIN, X0MAX, iterations, alpha_mult, nu, 0.001 * eps, 2*max_iter, 0.001 * max_diff, output);
						/*
						R_est = 1, q1 = 1, q2 = 1;
						XX(0, 0) = 80;
						XX(1, 0) = 10;
						XX(2, 0) = 2.9;
						XX(3, 0) = -91.2;
						XX(4, 0) = 0;
						*/
						/*
						R_est = 0.018850637712149891, q1 = 1.2572863375365346e-08, q2 = 2.9546186587375037e-06;
						R_est = 0.016459097618021223, q1 = 1.7629759189707618e-06, q2 = -1.8834178011400873e-06;
						XX(0, 0) = 80;
						XX(1, 0) = 90;
						XX(2, 0) = 10;
						XX(3, 0) = 0;
						XX(4, 0) = 20;
						*/
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, (*i_projections), aspect, R_est, q1, q2, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
							R_est, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, XX, Y, V, X0MIN, X0MAX, iterations, alpha_mult, nu, 0.001 * eps,  max_iter, 0.0001 * max_diff, output);
							
						//min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ4 <T>(nl_test, pl_reference, (*i_projections), aspect, R_est, q1, q2, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
						//	R_est, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, XX, Y, V, X0MIN, X0MAX, iterations, alpha_mult, 0.001*eps, 2 * max_iter, 0.001*max_diff, output);


						X(0, 0) = R_est;  X.replace(XX, 1, 0);
						XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);

						scale = sqrt(q1 * q1 + q2 * q2); alpha = atan2(q2, q1) * 180 / M_PI;
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					{
						X.replace(X0(0, 3, 0, 0), 1, 0);  XMIN.replace(X0MIN(0, 3, 0, 0), 1, 0); XMAX.replace(X0MAX(0, 3, 0, 0), 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);

						XMIN(5, 0) = -1.0e03;	XMAX(5, 0) = 1.0e3;
						XMIN(6, 0) = -1.0e03;	XMAX(6, 0) = 1.0e3;
						XMIN(7, 0) = 0;		XMAX(7, 0) = 1.0e8;

			
						min_cost = min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test, pl_reference, (*i_projections), aspect, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
							analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);

					}

					//Matrix <T> WD = diag(W); WD.print(output);

					const T map_scale2 = (analysis_parameters.analysis_method != NonLinearLeastSquaresRot2Method ? R_def / X(0, 0) : R_def / R_est);

					//Compute match ratio
					//TIndexList matched_points;
					//T match_ratio = getMatchRatioTissotIndicatrix(nl_test, pl_reference, trans(X), W, *i_projections, matched_points, CollectOn, map_scale, 30.0);
					//best_sample.setHomotheticTransformationPercMatch((unsigned int)match_ratio);
					//best_sample.setHomotheticTransformationMatchedPointsIndices ( matched_points );

					//Set residuals
					best_sample.setHomotheticTransformationRatio(min_cost);

					best_sample.setR(X(0, 0));
					best_sample.setLatP(X(1, 0));
					best_sample.setLonP(X(2, 0));
					best_sample.setLat0(X(3, 0));
					best_sample.setLat1(X(3, 0));
					best_sample.setLat2(X(5, 0));
					best_sample.setLon0(X(4, 0));
					best_sample.setC(X(5, 0));

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						best_sample.setAlpha(X(6, 0));
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						best_sample.setAlpha(alpha);

					best_sample.setScaleHomT(map_scale2);
					best_sample.setIterations(iterations);
					best_sample.setResidualEval(res_evaluation);
					best_sample.setProj(*i_projections);

					//Add to the list
					if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
						sl.push_back(best_sample);

					//Restore projection properties after analysis
					(*i_projections)->setR(R);
					(*i_projections)->setCartPole(cart_pole);
					(*i_projections)->setLat0(lat0);
					(*i_projections)->setLat1(lat1);
					(*i_projections)->setLat2(lat2);
					(*i_projections)->setLon0(lon0);
					(*i_projections)->setDx(dx);
					(*i_projections)->setDy(dy);
					(*i_projections)->setC(c);

					//Print successfully analyzed samples for one cartographic projection
					std::cout << " [" << iterations << " it, f = " << min_cost  << "] ";
					*output << " [" << iterations << " it, f = " << min_cost << "] ";
				}
			}

			std::cout << '\n';
			*output << '\n';
		}

		//Throw exception, analysis have not been computed
		catch (Exception & error)
		{
			if (analysis_parameters.print_exceptions)
				error.printException();
		}

		//Assign the amount of created samples
		total_created_or_thrown_samples += total_created_and_analyzed_samples_projection;
	}
}

/*
template <typename T>
void CartAnalysis::analyzeProjectionIncr(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
	unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
	//Incremental algorithm for the map projection analysis
	//Find mimum using the Minimum Least Squares
	const unsigned int m = nl_test.size();

	//Total successfully computed analysis for one cartographic projection
	unsigned int total_created_and_analyzed_samples_projection = 0;

	//Find extreme values
	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	//Create list of samples
	std::list <TSampleProjection<T> > analyzed_projections;

	//Create initial solution: process all cartographic projections from the list one by one
	for (typename TItemsList <Projection <T> *> ::Type::const_iterator i_projections = pl.begin(); i_projections != pl.end(); ++i_projections)
	{
		//Get limits of the cartographic pole latitude and longitude: some projections are defined only in normal position
		total_created_and_analyzed_samples_projection = 0;

		//Get defined radius
		T R_def = (*i_projections)->getR(), R_est = R_def;

		//Print actual projection name to the log
		std::cout << (*i_projections)->getName() << ": ";
		*output << (*i_projections)->getName() << ": ";

		//Get both latp and lonp intervals: lonp intervals are set to the moved central meridian (further must be reduced)
		TMinMax <T> latp_interval_heur = (*i_projections)->getLatPIntervalH(lat_interval);
		TMinMax <T> lonp_interval_heur = (*i_projections)->getLonPIntervalH(lon_interval);
		TMinMax <T> lat0_interval = (*i_projections)->getLat0Interval();

		const T latp_min = (*i_projections)->getLatPInterval().min_val;
		const T latp_max = (*i_projections)->getLatPInterval().max_val;
		const T lonp_min = (*i_projections)->getLonPInterval().min_val;
		const T lonp_max = (*i_projections)->getLonPInterval().max_val;
		const T lat0_min = lat0_interval.min_val;
		const T lat0_max = lat0_interval.max_val;

		//Create sample
		Sample <T> best_sample;

		//Create matrix of determined variables
		unsigned short n_par = 6;
		if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
			n_par = 7;
		else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
			n_par = 8;

		Matrix <unsigned int> IX = ones(2 * m, 1, 1.0);
		Matrix <T> X(n_par, 1), XMIN(n_par, 1), XMAX(n_par, 1), Y(2 * m, 1);

		//Set iteration parameters
		const T eps = 1.0e-10, max_diff = 1.0e-10, k = 2.5, ks = 0.01;
		T alpha_mult = 0.0001, nu = 0.0001;
		unsigned short iterations = 0;
		const unsigned int max_iter = 200;
		const T nu1 = 0.25, nu2 = 0.75, nu3 = 0.0001, gamma1 = 0.5, gamma2 = 2.0, lambda_min = 1.0e-6, lambda_max = 1.0e6;

		try
		{
			unsigned int total_samples_test = 0;
			Sample <T> sample_test;
			TAnalysisParameters <T> analysis_parameters_test(false);
			analysis_parameters_test.analysis_type.a_helt = true;

			//Store projection properties before analysis
			const T R = (*i_projections)->getR();
			const Point3DGeographic <T> cart_pole = (*i_projections)->getCartPole();
			const T lat0 = (*i_projections)->getLat0();
			const T lat1 = (*i_projections)->getLat1();
			const T lat2 = (*i_projections)->getLat2();
			const T lon0 = (*i_projections)->getLon0();
			const T dx = (*i_projections)->getDx();
			const T dy = (*i_projections)->getDy();
			const T c = (*i_projections)->getC();
			//const T alpha = ( *i_projections )->getAlpha();

			//Assign created samples amount
			const unsigned int total_created_and_analyzed_samples_projection_before = total_created_and_analyzed_samples_projection;

			//Compute mean of longitude to initialize lon0
			T lon_mean = 0;

			for (unsigned int i = 0; i < m; i++) lon_mean += pl_reference[i]->getLon();
			lon_mean /= m;

			//Conditions
			const bool normal_aspect_enabled = (analysis_parameters.analyze_normal_aspect) && ((*i_projections)->getCartPole().getLat() == MAX_LAT);
			const bool transverse_aspect_enabled = (analysis_parameters.analyze_transverse_aspect) && ((*i_projections)->getCartPole().getLat() == 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);
			const bool oblique_aspect_enabled = (analysis_parameters.analyze_oblique_aspect) && ((*i_projections)->getCartPole().getLat() != 0.0 || (*i_projections)->getCartPole().getLat() == MAX_LAT) &&
				((*i_projections)->getLatPInterval().min_val != (*i_projections)->getLatPInterval().max_val);

			//Compute lat0
			T lat0_init = 0.5 * (lat0_interval.min_val + lat0_interval.max_val);

			//Initialize matrices
			Matrix <T> X0(5, 1), X0MIN(5, 1), X0MAX(5, 1), XX(5, 1);

			//Common for all apects
			X0(2, 0) = lat0_init + 10;

			X0MIN(2, 0) = lat0_min;     X0MAX(2, 0) = lat0_max;

			//Conic projections: c = lat2
			if (typeid(*(*i_projections)) == typeid(ProjectionConic<T>))
			{
				X0MIN(4, 0) = 0;    X0MAX(4, 0) = lat0_max;
				X0(4, 0) = lat0_init + 20;
			}

			//Otherwise: c is unspecified constant
			else
			{
				X0MIN(4, 0) = 0.0;          X0MAX(4, 0) = 1.0e3;
				X0(4, 0) = 1;
			}

			//Proces available aspects
			TProjectionAspect aspect = NormalAspect;
			for (unsigned int i = 0; i < 3; i++)
			{
				//At least one test is set
				if ((i == 0) && (normal_aspect_enabled) || (i == 1) && (transverse_aspect_enabled) || (i == 2) && (oblique_aspect_enabled))
				{
					//Change analyzed aspect
					if (i == 1) aspect = TransverseAspect;
					else if (i == 2) aspect = ObliqueAspect;

					//Get current limits
					const T latp_init = 0.5 * (latp_interval_heur.min_val + latp_interval_heur.max_val);
					const T lonp_init = 0.5 * (lonp_interval_heur.min_val + lonp_interval_heur.max_val);

					//Create new meta-pole
					const Point3DGeographic <T> cart_pole_0(latp_init, lonp_init);

					//Set new position of the meta-pole for the transverse and oblique aspect
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole_0);

					//Use the cartometric analysis to set the initial value of the Earth radius: use the similarity transformation
					CartAnalysis::computeAnalysisForOneSample(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections, analysis_parameters_test, sample_test, false, total_samples_test, output);

					//Get initialize Earth radius using the similarity transformation
					T R0 = (*i_projections)->getR() / sample_test.getScaleHelT();

					//Restore meta-pole position
					if (i != 0)
						(*i_projections)->setCartPole(cart_pole);

					//Normal aspect
					if ((i == 0) && (normal_aspect_enabled))
					{
						//Initial solution for M7S method
						X0(0, 0) = MAX_LAT;
						X0(1, 0) = 0;
						X0(3, 0) = 10.0;

						//Set intervals for M7S method
						X0MIN(0, 0) = 90;           X0MAX(0, 0) = 90;
						X0MIN(1, 0) = 0;            X0MAX(1, 0) = 0;
						X0MIN(3, 0) = -MAX_LON;     X0MAX(3, 0) = MAX_LON;
					}

					else if ((i == 1) && (transverse_aspect_enabled))
					{
						//Initial solution for M7S method
						X0(0, 0) = 0;
						X0(1, 0) = lonp_init + 10;
						X0(3, 0) = 0;

						//Set intervals for M7S method
						X0MIN(0, 0) = 0;            X0MAX(0, 0) = 0;
						X0MIN(1, 0) = lonp_min;     X0MAX(1, 0) = lonp_max;
						X0MIN(3, 0) = 0;	    X0MAX(3, 0) = 0;
					}

					else if ((i == 2) && (oblique_aspect_enabled))
					{
						X0(0, 0) = latp_init + 20;
						X0(1, 0) = lonp_init + 10;
						X0(3, 0) = 0;

						//X0(0, 0) = -67; X0(1, 0) = 120;
						//X0(2, 0) = 50;

						X0MIN(0, 0) = latp_min;	X0MAX(0, 0) = latp_max;
						X0MIN(1, 0) = lonp_min;	X0MAX(1, 0) = lonp_max;
						X0MIN(3, 0) = 0;	X0MAX(3, 0) = 0;
					}


					//Create initial vector
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					{
						X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					{
						X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
						X(6, 0) = 30;  XMIN(6, 0) = -MAX_LON; XMAX(6, 0) = MAX_LON;
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					{
						X(0, 0) = R_est;  X.replace(XX, 1, 0);
						XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					{
						X.replace(X0(0, 3, 0, 0), 1, 0);  XMIN.replace(X0MIN(0, 3, 0, 0), 1, 0); XMAX.replace(X0MAX(0, 3, 0, 0), 1, 0);
						X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
					}

					//Add to the stack
					TSampleProjection <T> sample_proj(n_par);

					sample_proj.weight = 0;
					sample_proj.projection = *i_projections;
					sample_proj.aspect = aspect;
					sample_proj.X = X;
					sample_proj.XMIN = XMIN;
					sample_proj.XMAX = XMAX;

					//Add to the list
					analyzed_projections.push_back(sample_proj);
				}
			}

		}

		//Throw exception, analysis have not been computed
		catch (Exception & error)
		{
			if (analysis_parameters.print_exceptions)
				error.printException();
		}

	}
		//Create data structures for the incremental algorithm
		Container <Node3DCartesian <T> *> nl_test_incr;
		Container <Point3DGeographic <T> *> pl_reference_incr;

		//Create sample
		Sample <T> best_sample;

		//Add first 5 points to the lists

		//Run incremental algorithm
		T average_weight = 0;
		for (unsigned int i = 0; i < m; i++)
		{
			unsigned int l = i;
			for (; (l< i + 5) && (l < m); l++)
			{
				nl_test_incr.push_back(nl_test[l]);
				pl_reference_incr.push_back(pl_reference[l]);
			}

			i = l;

			std::cout <<" ******** Points " << i << " \n";
			
			//Add point to the lists
			nl_test_incr.push_back(nl_test[i]);
			pl_reference_incr.push_back(pl_reference[i]);

			//Process lists of samples
			std::list <TSampleProjection<T> > ::iterator i_analyzed_projections = analyzed_projections.begin();

			//Create matrix of determined variables
			unsigned short n_par = 6;
			if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
				n_par = 7;
			else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
				n_par = 8;

			Matrix <unsigned int> IX = ones(2 * m, 1, 1.0);
			Matrix <T> X(n_par, 1), XMIN(n_par, 1), XMAX(n_par, 1), Y(2 * m, 1);

			//Set iteration parameters
			const T eps = 1.0e-10, max_diff = 1.0e-10, k = 2.5, ks = 0.01;
			T alpha_mult = 0.0001, nu = 0.0001;
			unsigned short iterations = 0;
			const unsigned int max_iter = 200;
			const T nu1 = 0.25, nu2 = 0.75, nu3 = 0.0001, gamma1 = 0.5, gamma2 = 2.0, lambda_min = 1.0e-6, lambda_max = 1.0e6;

			//Process all analyzed projections
			T average_weight_new = 0;
			for (unsigned int j = 0; i_analyzed_projections != analyzed_projections.end(); ++i_analyzed_projections)
			{
				std::cout << (*i_analyzed_projections).projection->getName();

				//Get actual sample and its properties
				T weight = (*i_analyzed_projections).weight;
				std::cout << " " << weight << " / " << average_weight;

				//Process sample
				if (weight <= average_weight)
				{
					Projection <T> *proj = (*i_analyzed_projections).projection;
					TProjectionAspect aspect = (*i_analyzed_projections).aspect;
					Matrix <T> X = (*i_analyzed_projections).X;
					Matrix <T> XMIN = (*i_analyzed_projections).XMIN;
					Matrix <T> XMAX = (*i_analyzed_projections).XMAX;

					//Variables
					unsigned int res_evaluation = 0, jac_evaluation = 0;
					bool enable_additional_lon0_analysis = false;
					T x_mass_reference = 0.0, y_mass_reference = 0.0, min_cost = 0, alpha = 0, R_est = 1.0, scale = 1.0, q1 = 1, q2 = 1, k = 2.5;
					TMEstimatorsWeightFunction me_function = HuberFunction;

					//Matrices
					Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);
					Matrix <unsigned int> IX(2 * m, 1);

					Container <Node3DCartesianProjected <T> *> nl_projected;

					//Set initial solution and intervals depending ion the selected method
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					{
						//X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						//X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);

						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_incr, pl_reference_incr, proj, aspect, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
							analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, HuberFunction, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, 0.05*max_iter, max_diff, output);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					{
						//X.replace(X0, 1, 0); XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);
						//X(0, 0) = R0; XMIN(0, 0) = ks * X(0, 0); XMAX(0, 0) = 1.0 / ks * X(0, 0);
						//X(6, 0) = 30;  XMIN(6, 0) = -MAX_LON; XMAX(6, 0) = MAX_LON;

						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_incr, pl_reference_incr, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
							x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, 10, max_diff, output);
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					{
						
						XX = X0;

						//min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, (*i_projections), aspect, R_est, q1, q2, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
						//	R_est, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, XX, Y, V, X0MIN, X0MAX, iterations, alpha_mult, nu, 0.001 * eps, 2*max_iter, 0.001 * max_diff, output);

						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_incr, pl_reference_incr, proj, aspect, R_est, q1, q2, enable_additional_lon0_analysis, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
						R_est, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, res_evaluation, me_function, k, IX, enable_additional_lon0_analysis, output), FAnalyzeProjC <double>(), W, XX, Y, V, X0MIN, X0MAX, iterations, alpha_mult, nu, 0.001 * eps, 1 * max_iter, 0.0001 * max_diff, output);

						X(0, 0) = R_est;  X.replace(XX, 1, 0);
						XMIN.replace(X0MIN, 1, 0); XMAX.replace(X0MAX, 1, 0);

						scale = sqrt(q1 * q1 + q2 * q2); alpha = atan2(q2, q1) * 180 / M_PI;
						
					}

					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					{
						min_cost = min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_incr, pl_reference_incr, proj, aspect, analysis_parameters.print_exceptions, jac_evaluation), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
							analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					}

					if (min_cost < 100)
					{
						//Update sample
						(*i_analyzed_projections).weight = min_cost;
						(*i_analyzed_projections).X = X;
						//X.print();

						std::cout << " / " << min_cost << '\n';

						//Compute new average weight
						if (j > 0) average_weight_new = (average_weight_new * j + min_cost) / (j + 1);

						//Increment j
						j++;
					}

					//Remove sample
					else
					{
						i_analyzed_projections = analyzed_projections.erase(i_analyzed_projections);
						++i_analyzed_projections;
					}
				}

				//Remove sample from the list
				else
				{
					std::cout << "*********** Deleted \n";
					i_analyzed_projections = analyzed_projections.erase(i_analyzed_projections);
					++i_analyzed_projections;
				}
			}


			//Assign new average cost
			average_weight = average_weight_new;
		}

		//Assign the amount of created samples
		total_created_or_thrown_samples += total_created_and_analyzed_samples_projection;
	

	//Create priority queue, add all analyzed samples
}
*/


template <typename T>
void CartAnalysis::createOptimalLatPLonPPositions(const Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, const TMinMax <T> &latp_interval_heur, const TMinMax <T> &lonp_interval_heur, const TAnalysisParameters <T> & analysis_parameters,
	const TProjectionAspect & proj_aspect, typename TItemsList <TProjectionPolePosition<T> >::Type & proj_pole_positions_list, std::ostream * output)
{
	//Create list of latp, lonp poitions with respect to complex criteriJ (used only first 20% of results)
	//TProjectionAspect = NormalAspect create latp, lonp positions for a normal aspect
	//TProjectionAspect = TransverseAspect = create latp, lonp positions for a transverse aspect
	//TProjectionAspect = ObliqueAspect = create latp, lonp positions for an oblique aspect
	T complex_crit_sum = 0.0;

	const bool lat0_set = (proj->getLat0() != 0.0) && (proj->getLat0() != 45.0);
	const bool latp_set = proj->getCartPole().getLat() != MAX_LAT;
	const bool lonp_set = proj->getCartPole().getLon() != 0.0;

	//Set latp for a projection aspect
	T latp = MAX_LAT;																							//Normal apect

	if (proj_aspect == TransverseAspect) latp = 0.0;																			//Transverse aspect
	else if (proj_aspect == ObliqueAspect)  latp = proj->getLatPInterval().min_val;															//Oblique aspect

	//Load from the configuration file
	if (latp_set) latp = proj->getCartPole().getLat();

	//Set latp_min, latp_max for a projection aspect
	T latp_min = MAX_LAT, latp_max = MAX_LAT;																				//Normal aspect

	if (proj_aspect == TransverseAspect)																					//Transverse aspect
	{
		latp_min = 0.0; latp_max = latp_min;
	}

	else if (proj_aspect == ObliqueAspect)																				//Oblique aspect
	{
		latp_min = proj->getLatPInterval().min_val; latp_max = proj->getLatPInterval().max_val;
	}

	//Set lonp_min, lonp_max for a projection aspect
	T lonp_min = 0.0, lonp_max = 0.0;																					//Normal aspect

	if (proj_aspect == TransverseAspect || proj_aspect == ObliqueAspect)																	//Transverse or oblique aspect
	{
		lonp_min = proj->getLonPInterval().min_val; lonp_max = proj->getLonPInterval().max_val;
	}

	//Process normal / transverse / oblique aspect of the map projection
	for (; (proj_aspect == TransverseAspect ? latp == 0.0 : (latp >= latp_min) && (latp <= latp_max)); latp += analysis_parameters.latp_step)
	{
		//Set lonp for a projection aspect
		T lonp = 0.0;																							//Normal aspect

		if ((latp != MAX_LAT) && (proj_aspect == TransverseAspect || proj_aspect == ObliqueAspect))												//Transverse or oblique aspect
			lonp = (lonp_set) ? proj->getCartPole().getLon() : proj->getLonPInterval().min_val;												//Load from configuration file or use default value

		//lonp = -90.0;
		for (; (latp == MAX_LAT ? lonp == 0.0 : (lonp >= lonp_min) && (lonp <= lonp_max)); lonp += analysis_parameters.lonp_step)
		{
			//Test, if a generated lonp value satisfies a condition of heuristic, for a normal aspect do no test
			if ((proj_aspect == NormalAspect) && (latp == MAX_LAT) ||																//Normal aspect
				(proj_aspect == TransverseAspect || (proj_aspect == ObliqueAspect) && ((latp != MAX_LAT) && (latp != 0.0))) &&						//Oblique aspect ( - normal and - transverse aspects ) or transverse aspect
				((lonp_interval_heur.min_val < lonp_interval_heur.max_val) && ((lonp >= lonp_interval_heur.min_val) && (lonp <= lonp_interval_heur.max_val)) ||		//Inside lonp interval given by heuristic
				(lonp_interval_heur.min_val > lonp_interval_heur.max_val) && ((lonp >= lonp_interval_heur.min_val) || (lonp <= lonp_interval_heur.max_val))) &&		//Outside lonp interval given by heuristic
				((latp >= latp_interval_heur.min_val) && (latp <= latp_interval_heur.max_val)))
			{
				//Set lat0: rememeber old value lat0
				const T lat0_old = proj->getLat0();

				//Process all undistorted meridians for latp, lonp positions
				for (T lat0 = (lat0_set ? lat0_old : proj->getLat0Interval().min_val); (lat0 <= proj->getLat0Interval().max_val); lat0 += analysis_parameters.lat0_step)
				{
					//Set lat0 to compute correct Tissot indicatrix
					proj->setLat0(lat0);

					//Compute coordinates of all geographic points points in sample's projection + complex criterium
					const unsigned int n_ref = pl_reference.size();
					T lat_lon_mbr[] = { MAX_LAT, MAX_LON, MIN_LAT, MIN_LON };

					for (unsigned int i = 0; i < n_ref; i++)
					{
						try
						{
							//Get type of the direction
							TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

							//Convert geographic point to oblique aspect
							const T lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), pl_reference[i]->getLon(), latp, lonp);
							const T lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), pl_reference[i]->getLon(), lat_trans, latp, lonp, trans_lon_dir);

							//Find extreme values
							if (lat_trans < lat_lon_mbr[0]) lat_lon_mbr[0] = lat_trans;
							else if (lat_trans > lat_lon_mbr[2]) lat_lon_mbr[2] = lat_trans;

							if (lon_trans < lat_lon_mbr[1]) lat_lon_mbr[1] = lon_trans;
							else if (lon_trans > lat_lon_mbr[3]) lat_lon_mbr[3] = lon_trans;
						}

						//Throw exception: a conversion error or bad data
						catch (Exception & error)
						{
							if (analysis_parameters.print_exceptions) error.printException(output);
						}
					}

					//Compute complex criterium
					T complex_crit = 0.0, weight_sum = 0.0;

					//Compute complex criteria in extreme points
					if (analysis_parameters.perform_heuristic)
					{
						for (unsigned int i = 0; i < 2; i += 2)
						{
							//Create temporary point
							Point3DGeographic <T> p_oblique_temp(lat_lon_mbr[i], lat_lon_mbr[i + 1]);

							T h = 1.0, k = 1.0;

							try
							{
								//Compute Tissot indiatrix parameters
								h = CartDistortion::H(NUM_DERIV_STEP, &p_oblique_temp, proj, analysis_parameters.print_exceptions);
								k = CartDistortion::K(NUM_DERIV_STEP, &p_oblique_temp, proj, analysis_parameters.print_exceptions);

								//Switch h <-> k
								if (h < k)
								{
									const T temp = h; h = k; k = temp;
								}
							}

							//Throw exception: can not compute Tissot indicatrix
							catch (Exception & error)
							{
								if (analysis_parameters.print_exceptions) 	error.printException(output);
							}

							//Compute complex criteria
							const T weight = cos(M_PI / 180.0 *  lat_lon_mbr[i]);
							complex_crit += (0.5 * (fabs(h - 1.0) + fabs(k - 1.0)) + h / k - 1.0) * weight;
							weight_sum += weight;
						}

						//Complex criteria
						complex_crit /= weight_sum;
						complex_crit_sum += complex_crit;
					}

					//Create new possible pole position
					TProjectionPolePosition <T> latp_lonp(latp, lonp, lat0, complex_crit);

					//Add to the list of possible pole positions
					proj_pole_positions_list.push_back(latp_lonp);

					//Projection has set lat0 in a definition file
					if (lat0_set) break;
				}

				//Set old lat0
				proj->setLat0(lat0_old);
			}

			//Projection has set lonp in a definition file
			if (lonp_set) break;
		}

		//Projection has set latp in a definition file
		if (latp_set) break;
	}

	//Remove inappropriate pole positions
	if ((analysis_parameters.perform_heuristic) && (proj_pole_positions_list.size() > 10))
	{
		//Erase all worse than 2 *mean( complex_crit )
		proj_pole_positions_list.erase(std::remove_if(proj_pole_positions_list.begin(), proj_pole_positions_list.end(), removeProjectionPolePositions <T>(2.0 * complex_crit_sum / proj_pole_positions_list.size())), proj_pole_positions_list.end());

		//Sort all created  positions according to cartographic pole values
		std::sort(proj_pole_positions_list.begin(), proj_pole_positions_list.end(), sortProjectionPolePositionsByLat());
	}
}


template <typename T>
void CartAnalysis::findLatPLonPIntervals(const Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TMinMax <T> &latp_interval_heur, TMinMax <T> &lonp_interval_heur)
{
	//Find optimal latp and lonp intervals with respect to geographic properties of the analyzed area and a type of the projection
	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	//Find, which intervals are covered by the analyzed territory: i1 = (-180, -90), i2 = (-90, 0), i3 = (0,90), i4 = (90, 180)
	bool i1 = false, i2 = false, i3 = false, i4 = false;

	for (unsigned int i = 0; (i < pl_reference.size()) && (!(i1 && i2 && i3 && i4)); i++)
	{
		if ((pl_reference[i]->getLon() > MIN_LON) && (pl_reference[i]->getLon() < -90.0)) i1 = true;
		else if ((pl_reference[i]->getLon() > -90.0) && (pl_reference[i]->getLon() < 0.0)) i2 = true;
		else if ((pl_reference[i]->getLon() > 0.0) && (pl_reference[i]->getLon() < 90.0)) i3 = true;
		else if ((pl_reference[i]->getLon() > 90.0) && (pl_reference[i]->getLon() < MAX_LON)) i4 = true;
	}

	//Test, if analyzed area is small. Up to 2 intervals are covered ( no three nor four intervals covered )
	if ((!(i1 && i2 && i3 && i4)) && (!(i1 && i2 && i3)) && (!(i2 && i3 && i4)) &&
		(!(i3 && i4 && i1)) && (!(i4 && i1 && i2)) && (latp_interval_heur.min_val != MAX_LAT))
	{
		//Compute heuristic intervals
		TMinMax <T> latp_interval_oblique = proj->getLatPIntervalH(lat_interval);
		TMinMax <T> lonp_interval_oblique = proj->getLonPIntervalH(lon_interval);

		//Switch min and max if itervals i1 and i4 are covered
		if (i1 && i4)
		{
			T temp = lonp_interval_oblique.min_val;
			lonp_interval_oblique.min_val = lonp_interval_oblique.max_val;
			lonp_interval_oblique.max_val = temp;
		}

		//Round values: min down, max up, to 10 deg
		lonp_interval_oblique.min_val = ((int)(lonp_interval_oblique.min_val / 10.0)) * 10.0;
		lonp_interval_oblique.max_val = ((int)(lonp_interval_oblique.max_val / 10.0 + 0.5)) * 10.0;

		//Assign intervals
		latp_interval_heur = latp_interval_oblique;
		lonp_interval_heur = lonp_interval_oblique;
	}
}



template <typename T>
T CartAnalysis::computeAnalysisForOneSample(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type & meridians, typename TParallelsList <T> ::Type & parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, const TAnalysisParameters <T> & analysis_parameters, Sample <T> &sample_res, bool singular_points_found, unsigned int & created_samples, std::ostream * output)
{
	//Compute all cartometric analysis for one sample
	unsigned int n_nsing = nl_test.size(), n_best = n_nsing;
	bool outliers_found = false;

	//Set pointers to processed non sets
	Container <Point3DGeographic <T> *> * p_pl_reference = &pl_reference;
	Container <Point3DGeographic <T> *> *p_pl_reference_non_sing = p_pl_reference;
	Container <Node3DCartesian <T> *> *p_nl_test_non_sing = &nl_test;

	//Create temporary containers for non singular points: containers of points
	Container <Point3DGeographic <T> *> pl_reference_red;
	Container <Node3DCartesian <T> *> nl_test_non_sing;
	Container <Point3DGeographic <T> *> pl_reference_non_sing;

	//Set pointers to meridians and parallels
	typename TMeridiansList <T> ::Type * p_meridians_non_sing = &meridians;
	typename TParallelsList <T> ::Type * p_parallels_non_sing = &parallels;

	//Create temporary meridians and parallels
	typename TMeridiansList <T> ::Type meridians_non_sing;
	typename TParallelsList <T> ::Type parallels_non_sing;

	//List of singular points
	TIndexList non_singular_points;

	for (unsigned int j = 0; j < n_nsing; j++) non_singular_points.push_back(j);

	//Reduce lon using a new central meridian redefined in projection file, if necessary
	if (proj->getLon0() != 0.0)
	{
		//Reduce lon
		redLon(pl_reference, proj->getLon0(), pl_reference_red);

		//Set pointer to processed file: reduced or non reduced
		p_pl_reference_non_sing = &pl_reference_red;
	}

	//Remove singular points, store non singular pairs
	typename TDevIndexPairs <T>::Type non_singular_pairs;
	removeSingularPoints(*p_nl_test_non_sing, *p_pl_reference_non_sing, proj, nl_test_non_sing, pl_reference_non_sing, non_singular_pairs);

	//Set singular point flag
	singular_points_found = false;

	//Singular points found: remove points, correct meridians and parallels
	n_nsing = nl_test_non_sing.size();

	if (nl_test.size() != n_nsing)
	{
		//Set flag to true, used for all samples using non-singular sets
		singular_points_found = true;

		//Create copy of meridians / parallels
		meridians_non_sing = meridians; parallels_non_sing = parallels;

		//Correct meridians and parallels
		correctMeridiansAndParrallels <T>(meridians_non_sing, parallels_non_sing, non_singular_pairs);

		//Convert non singular pairs to index list: indices will be printed in output
		non_singular_points.clear();
		std::transform(non_singular_pairs.begin(), non_singular_pairs.end(), std::back_inserter(non_singular_points), getSecondElementInPair());

		//Set pointers to non-singular meridians and parallels
		p_meridians_non_sing = &meridians_non_sing; p_parallels_non_sing = &parallels_non_sing;

		//Set pointers to newly created non-singular sets
		p_nl_test_non_sing = &nl_test_non_sing; p_pl_reference_non_sing = &pl_reference_non_sing;
	}

	//Set non-singular points
	sample_res.setNonSingularPointsIndices(non_singular_points);

	//Create empty list of projected points and transformed points
	Container <Node3DCartesianProjected <T> *> nl_projected;

	//Get rotation
	const T alpha = sample_res.getAlpha();

	//Compute coordinates of all geographic points points in sample's projection and add to the list
	for (unsigned int i = 0; i < n_nsing; i++)
	{
		try
		{
			//Coordinates of the current point
			T lat = (*p_pl_reference_non_sing)[i]->getLat();
			T lon = (*p_pl_reference_non_sing)[i]->getLon();

			//Cartographic pole
			const T latp = proj->getCartPole().getLat();
			const T lonp = proj->getCartPole().getLon();

			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Reduce lon
			T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), lon);

			T x = 0.0, y = 0.0, lat_trans = lat, lon_trans = lon;

			//More safe lat/lon conversion
			for (unsigned int j = 0; j < 3; j++)
			{
				try
				{
					//Convert geographic point to oblique aspect
					lat_trans = CartTransformation::latToLatTrans(lat, lon_red, latp, lonp);
					lon_trans = CartTransformation::lonToLonTrans(lat, lon_red, latp, lonp, trans_lon_dir);
				}

				//2 attempt to avoid the singularity
				catch (Exception &error)
				{
					//Move in latitude direction
					if (j == 0)
					{
						if (lat == MAX_LAT || lat == latp)
							lat -= GRATICULE_ANGLE_SHIFT;
						else
							lat += GRATICULE_ANGLE_SHIFT;
					}

					//Move in longitude direction
					else if (j == 1)
					{
						if (lon_red == MAX_LON || lon == lonp)
							lon_red -= GRATICULE_ANGLE_SHIFT;
						else
							lon_red += GRATICULE_ANGLE_SHIFT;
					}

					//Neither first nor the second shhifts do not bring improvement
					else if (j == 2)
					{
						std::cout << lat << "  " << lon_red;
						throw;
					}
				}
			}
			
			//Create temporary point
			Point3DGeographic <T> p_oblique_temp(lat_trans, lon_trans);

			for (unsigned int j = 0; j < 3; j++)
			{
				try
				{
					const T x_temp = CartTransformation::latLonToX(&p_oblique_temp, proj, analysis_parameters.print_exceptions);
					const T y_temp = CartTransformation::latLonToY(&p_oblique_temp, proj, analysis_parameters.print_exceptions);

					//Involved rotation
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					{
						x = x_temp * cos(alpha * M_PI / 180) - y_temp * sin(alpha * M_PI / 180);
						y = x_temp * sin(alpha * M_PI / 180) + y_temp * cos(alpha * M_PI / 180);
					}

					//Without rotation
					else
					{
						x = x_temp;
						y = y_temp;
					}
				}

				//2 attempt to avoid the singularity
				catch (MathException <T> &error)
				{
					//Move in latitude direction
					if (j == 0)
					{
						if (lat_trans == MAX_LAT) 
							lat_trans -= GRATICULE_ANGLE_SHIFT;
						else 
							lat_trans += GRATICULE_ANGLE_SHIFT;
					}

					//Move in longitude direction
					else if (j == 1)
					{
						if (lon_trans == MAX_LON) 
							lon_trans -= GRATICULE_ANGLE_SHIFT;
						else 
							lon_trans += GRATICULE_ANGLE_SHIFT;
					}

					//Neither first nor the second shhifts do not bring improvement
					else if (j == 2)
					{
						throw;
					}

					//Set new geographic coordinates of the point
					p_oblique_temp.setLat(lat_trans);
					p_oblique_temp.setLon(lon_trans);
				}
			}

			//Create new point
			Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y, 0.0, 0.0, 0.0, 0.0, 0.0, TTissotIndicatrix <T>(), 0.0);
			
			//Uncertainty regions: compute Tissot Indicatrix parameters only for perspective samples ( will be faster ) for homothetic and helmert transformations
			if ((analysis_parameters.analysis_type.a_homt || analysis_parameters.analysis_type.a_helt) && (analysis_parameters.match_method == MatchTissotIndicatrix))
			{
				//Create Tissot structure, initialize with default parameters independent on (lat, lon): a = 1, XMAX = 1, Ae = 0
				TTissotIndicatrix <T> tiss;

				//Compute Tissot indicatrix parameters for (lat, lon)
				try
				{
					tiss = CartDistortion::Tiss(NUM_DERIV_STEP, &p_oblique_temp, proj, analysis_parameters.print_exceptions);
				}

				//If impossible to compute Tissot indicatrix parameters for (lat, lon), use default parameters (circle a=1, XMAX=1)
				catch (Exception & error)
				{
					if (analysis_parameters.print_exceptions) error.printException(output);
				}

				//Set Tissot Indicatrix parameters for the point
				n_projected->setTiss(tiss);
			}
			
			//Create new cartographic point
			nl_projected.push_back(n_projected);
		}

		//Throw exception: unspecified exception
		catch (Exception & error)
		{
			if (analysis_parameters.print_exceptions) error.printException(output);
		}
	}
	
	//Remove duplicate elements from reference data set (projected points)
	///nl_projected.removeDuplicateElements(nl_projected.begin(), nl_projected.end(),
	//	sortPointsByX(), isEqualPointByPlanarCoordinates <Node3DCartesianProjected <T> *>());

	//Both datasets do not contain the samee amount of points
	if (n_nsing != nl_projected.size())
	{
		throw BadDataException("BadDataException: both datasets contain a different number of points. ", "Sample had been thrown...");
	}
	//Set pointers both to test and projected files with non_singular items
	Container <Node3DCartesian <T> *> * p_nl_test_best = p_nl_test_non_sing;
	Container <Node3DCartesianProjected <T> *> * p_nl_projected_best = &nl_projected;

	//Create temporary containers for k-best fit points: containers of points
	Container <Node3DCartesian <T> *> nl_test_best;
	Container <Node3DCartesianProjected <T> *> nl_projected_best;
	//List of k-best points
	TIndexList k_best_points;

	for (unsigned int i = 0; i < n_nsing; i++) k_best_points.push_back(i);

	//Set pointers to meridians and parallels
	typename TMeridiansList <T> ::Type * p_meridians_best = p_meridians_non_sing; typename TParallelsList <T> ::Type * p_parallels_best = p_parallels_non_sing;

	//Create temporary meridians and parallels
	typename TMeridiansList <T> ::Type meridians_best; typename TParallelsList <T> ::Type parallels_best;

	//List of k_best pairs and multiplication ratio
	typename TDevIndexPairs <T>::Type min_pairs;
	/*
	//Remove outliers
	if (analysis_parameters.remove_outliers)
	{
	T f_init = 0, f = 0;
	unsigned int iter = 0;
	const unsigned int max_iter = 20;
	const T k = 2.0, eps = 1.0e-10;
	TTransformationKeyHelmert2D <T> min_key;

	//Find k-best fit transformation key using Homothetic transformation
	Outliers::findOutliersME(*p_nl_test_non_sing, nl_projected, nl_test_best, nl_projected_best, k, eps, SimilarityScheme, YangFunction, max_iter, min_pairs, f_init, f, iter);
	//Outliers::findOutliersIRLS(*p_nl_test_non_sing, nl_projected, nl_test_best, nl_projected_best, min_key, min_pairs);

	//Any outlier found
	n_best = nl_projected_best.size();

	if (n_nsing != n_best)
	{
	//Set flag: outliers found
	outliers_found = true;

	//Create copy of meridians and parallels with no outliers
	meridians_best = *p_meridians_non_sing; parallels_best = *p_parallels_non_sing;

	//Correct meridians and parallels: remove inappropriate points from points indices of all meridians and parallels
	correctMeridiansAndParrallels <T>(meridians_best, parallels_best, min_pairs);

	//Convert k-best pairs to index list: indices will be printed in output
	k_best_points.clear();
	std::transform(min_pairs.begin(), min_pairs.end(), std::back_inserter(k_best_points), getSecondElementInPair());

	//Set pointers to meridians and parallels
	p_meridians_best = &meridians_best; p_parallels_best = &parallels_best;

	//Set pointers to newly created without outliers
	p_nl_test_best = &nl_test_best; p_nl_projected_best = &nl_projected_best;
	}
	}
	*/

	//Set  k-best points
	sample_res.setKBestPointsIndices(k_best_points);

	//Compare shape of equator, meridian and north / south pole using turning function, similarity transformation: this is a heuristic throwing unperspective samples
	T sample_cost = MAX_FLOAT;

	if ((analysis_parameters.perform_heuristic) && (checkSample(*p_meridians_best, *p_parallels_best, *p_nl_test_best, *p_nl_projected_best, analysis_parameters.heuristic_sensitivity_ratio)) ||
		(!analysis_parameters.perform_heuristic))
	{
		//Compute multiplication ratio
		float mult_ratio = 2.0 - n_best / (float)n_nsing;

		//Set properties for the resulting smple
		sample_res.setProj(proj);
		sample_res.setR(proj->getR());
		sample_res.setLatP(proj->getCartPole().getLat());
		sample_res.setLonP(proj->getCartPole().getLon());
		sample_res.setLat0(proj->getLat0());
		sample_res.setLon0(proj->getLon0());
		sample_res.setDx(proj->getDx());
		sample_res.setDy(proj->getDy());
		sample_res.setSingularPointsFound(singular_points_found);
		sample_res.setOutliersFound(outliers_found);
		sample_res.setNonSingularPointsIndices(non_singular_points);
		sample_res.setKBestPointsIndices(k_best_points);

		//Set pointer to a sample
		Sample <T> *p_sample = &sample_res;
		//2D Helmert transformation: results for rotated and unrotated sample
		if (analysis_parameters.analysis_type.a_helt)
			analyzeSampleHelmertTransformationDeviation(*p_sample, *p_nl_test_best, *p_nl_projected_best, analysis_parameters.match_method, mult_ratio);

		//Create list for rotated test points and rotated sample as the copy of the result sample
		Container <Node3DCartesian <T> *> nl_test_best_rot;
		Sample <T> sample_rot(sample_res);

		//Process unrotated and rotated sample
		bool rotated_sample = false;

		//for (unsigned int j = 0; j < 2; j++)
		{
			//2D homothetic transformation:
			if (analysis_parameters.analysis_type.a_homt)
				analyzeSampleHomotheticTransformationDeviation(*p_sample, *p_nl_test_best, *p_nl_projected_best, analysis_parameters.match_method, mult_ratio);

			//Cross nearest distance
			if (analysis_parameters.analysis_type.a_cnd)
				analyzeSampleCrossNearestNeighbourDistance(*p_sample, *p_nl_test_best, *p_nl_projected_best, mult_ratio);

			//Analysis of graticule: turning function
			if (analysis_parameters.analysis_type.a_gn_tf)
				analyzeSampleGeographicNetworkTurningFunctionRatio(*p_sample, *p_nl_test_best, *p_nl_projected_best, *p_meridians_best, *p_parallels_best, mult_ratio);

			//Analysis of Voronoi diagrams: turning function
			if (analysis_parameters.analysis_type.a_vd_tf)
				analyzeSampleUsingVoronoiDiagramTurningFunctionRatio2(*p_sample, *p_nl_test_best, *p_nl_projected_best, faces_test, analysis_parameters, mult_ratio);

			//Angle of rotation given by 2D Helmert transformation
			const T rot_angle = p_sample->getRotation();

			//Is there a significant improvement for rotated sample: homt / helt > IMPROVE_RATIO_STD_DEV
			//Is the angle of rotation in  interval (88, 92) + k*M_PI/2
			rotated_sample = (analysis_parameters.correct_rotation) && (IMPROVE_RATIO_STD_DEV * sample_res.getHelmertTransformationRatio() < sample_res.getHomotheticTransformationRatio()) &&
				((short)(fabs(rot_angle) + REM_DIV_ROT_ANGLE) % 90 < 2 * REM_DIV_ROT_ANGLE) && (fabs(rot_angle) > MAX_LAT - REM_DIV_ROT_ANGLE);
			/*
			//Rotate a sample by the -( angle )
			if ((rotated_sample) && (!p_sample->getRotatedSample()))
			{
			//Create copy of test points
			nl_test_best_rot = *p_nl_test_best;

			//Rotate test points by angle given by 2D Helmert transformation in interval (80,100) + 0.5*k*M_PI
			for (unsigned int k = 0; k < n_best; k++)
			{
			nl_test_best_rot[k]->setX((*p_nl_test_best)[k]->getX() * cos(rot_angle * M_PI / 180.0) - (*p_nl_test_best)[k]->getY() * sin(rot_angle * M_PI / 180.0));
			nl_test_best_rot[k]->setY((*p_nl_test_best)[k]->getX() * sin(rot_angle * M_PI / 180.0) + (*p_nl_test_best)[k]->getY() * cos(rot_angle * M_PI / 180.0));
			}

			//Set pointer to rotated sample
			p_sample = &sample_rot;
			p_nl_test_best = &nl_test_best_rot;

			//Set angle and rotation flag for a sample
			p_sample->setRotatedSample(true);
			}

			else break;
			*/
		}

		//Compute sample cost
		sample_cost = (rotated_sample ? sample_rot.getSampleCost(analysis_parameters.analysis_type) : sample_res.getSampleCost(analysis_parameters.analysis_type));

		//Increment successfully created samples for projection
		created_samples++;

		//Add also rotated sample to the list
		if (rotated_sample)
		{
			sample_res = sample_rot;
			created_samples++;
		}
	}

	//Return cost of the sample
	return sample_cost;
}


template <typename T>
bool CartAnalysis::checkSample(const typename TMeridiansList <T> ::Type &meridians, const typename TParallelsList <T> ::Type &parallels, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected <T> *> &nl_projected, const T heuristic_sensitivity_ratio)
{
	//Small heuristic for sample: compare shape or the prime meridian, equator (central meridian, central parallel) north pole and south pole for test and reference points
	//using turning function und analyze Helmert transformation (fast match using circles)
	bool prime_meridian_found = false, equator_found = false;

	//Create list of transformed points
	Container <Node3DCartesian <T> *> nl_transformed;

	//Analyze match ratio using Helmert transformation ( additional test if no meridian and parallel have been found )
	TTransformationKeyHelmert2D <T> key_helmert;
	HelmertTransformation2D::transformPoints(nl_projected, nl_test, nl_transformed, key_helmert);

	//Sample is too rotated
	const T rot_angle = atan2(key_helmert.c2, key_helmert.c1) * 180.0 / M_PI;

	if ((short)(fabs(rot_angle) + 3.0 * REM_DIV_ROT_ANGLE) % 90 > 6.0 * REM_DIV_ROT_ANGLE)
	{
		return false;
	}

	//Test match ratio: at least 75 percent of points matched (fast test with the circle)
	TIndexList matched_points;

	if (getMatchRatioCircle(nl_projected, nl_transformed, matched_points, CollectOff, MATCHING_FACTOR * heuristic_sensitivity_ratio) < 75)
	{
		return false;
	}

	//Process all meridians and find prime meridian
	for (typename TMeridiansList <T> ::Type::const_iterator i_meridians = meridians.begin(); i_meridians != meridians.end(); ++i_meridians)
	{
		//Find prime meridian
		if ((*i_meridians).getLon() == 0)
		{
			//Convert test meridian to Points 2D list
			Container <Point3DCartesian <T> > pl_meridian_test(nl_test, (*i_meridians).getPointsIndices());

			//Convert projected meridian to Points2D list
			Container <Point3DCartesian <T> > pl_meridian_projected(nl_projected, (*i_meridians).getPointsIndices());

			//Compute turning function difference for each test and projected meridian
			T turning_function_ratio_prime_meridian = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_meridian_test, pl_meridian_projected, RotationInvariant, ScaleInvariant);

			//Both prime meridians are not similar
			if (turning_function_ratio_prime_meridian > TURNING_FUNCTION_MAX_DIFFERENCE * pl_meridian_projected.size() * heuristic_sensitivity_ratio)
			{
				return false;
			}

			//Set prime meridian as found
			prime_meridian_found = true;
		}
	}

	//Process all parallels find equator, north pole, south pole
	for (typename TParallelsList <T> ::Type::const_iterator i_parallels = parallels.begin(); i_parallels != parallels.end(); ++i_parallels)
	{
		//Find equator
		if ((*i_parallels).getLat() == 0)
		{
			//Convert test equator to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_test(nl_test, (*i_parallels).getPointsIndices());

			//Convert projected equator to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_projected(nl_projected, (*i_parallels).getPointsIndices());

			//Compute turning function difference for each parallel
			T turning_function_ratio_equator = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_parallel_test, pl_parallel_projected, RotationInvariant, ScaleInvariant);

			//Both equators are not similar
			if (turning_function_ratio_equator > TURNING_FUNCTION_MAX_DIFFERENCE * pl_parallel_projected.size() * heuristic_sensitivity_ratio)
			{
				return false;
			}

			//Set equator as found
			equator_found = true;
		}

		//Find and analyze north pole
		if ((*i_parallels).getLat() == MAX_LAT)
		{
			//Convert test parallel to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_test(nl_test, (*i_parallels).getPointsIndices());

			//Convert projected parallel to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_projected(nl_projected, (*i_parallels).getPointsIndices());

			//Compute turning function difference for each parallel
			T turning_function_ratio_north_pole = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_parallel_test, pl_parallel_projected, RotationInvariant, ScaleInvariant);

			//Both north poles are not similar
			if (turning_function_ratio_north_pole > TURNING_FUNCTION_MAX_DIFFERENCE * pl_parallel_projected.size() * heuristic_sensitivity_ratio)
			{
				return false;
			}
		}

		//Find and analyze south pole
		if ((*i_parallels).getLat() == MIN_LAT)
		{
			//Convert test parallel to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_test(nl_test, (*i_parallels).getPointsIndices());

			//Convert projected parallel to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_projected(nl_projected, (*i_parallels).getPointsIndices());

			//Compute turning function difference for each parallel
			T turning_function_ratio_south_pole = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_parallel_test, pl_parallel_projected, RotationInvariant, ScaleInvariant);

			//Both south poles are not similar
			if (turning_function_ratio_south_pole > TURNING_FUNCTION_MAX_DIFFERENCE * pl_parallel_projected.size() * heuristic_sensitivity_ratio)
			{
				return false;
			}
		}
	}

	//If not prime meridian found, continue with the found central meridian of the analyzed area
	if ((!prime_meridian_found) && (meridians.size() > 0))
	{
		//Set central meridian of the dataset
		typename TMeridiansList <T> ::Type::const_iterator i_meridians = meridians.begin();

		for (unsigned int i = 0; i < meridians.size() / 2; ++i_meridians, ++i) {}

		//Convert test meridian to Points 2D list
		Container <Point3DCartesian <T> > pl_meridian_test(nl_test, (*i_meridians).getPointsIndices());

		//Convert projected meridian to Points2D list
		Container <Point3DCartesian <T> > pl_meridian_projected(nl_projected, (*i_meridians).getPointsIndices());

		//Compute turning function difference for each test and projected meridian
		T turning_function_ratio_meridian = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_meridian_test, pl_meridian_projected, RotationInvariant, ScaleInvariant);

		//Both prime meridians are not similar
		if (turning_function_ratio_meridian > TURNING_FUNCTION_MAX_DIFFERENCE * pl_meridian_projected.size() * heuristic_sensitivity_ratio)
		{
			return false;
		}
	}

	//If not equator found, continue with the found central parallel of the analyzed area
	if ((!equator_found) && (parallels.size() > 0))
	{
		//Set central parallel of the dataset
		typename TParallelsList <T> ::Type::const_iterator i_parallels = parallels.begin();

		for (unsigned int i = 0; i < parallels.size() / 2; ++i_parallels, ++i) {}

		//Convert test parallel to Points2D list
		Container <Point3DCartesian <T> > pl_parallel_test(nl_test, (*i_parallels).getPointsIndices());

		//Convert projected parallel to Points2D list
		Container <Point3DCartesian <T> > pl_parallel_projected(nl_projected, (*i_parallels).getPointsIndices());

		//Compute turning function difference for each parallel
		T turning_function_ratio_parallel = TurningFunction::compare2PolyLinesUsingTurningFunction(pl_parallel_test, pl_parallel_projected, RotationInvariant, ScaleInvariant);

		//Both equators are not similar
		if (turning_function_ratio_parallel >  TURNING_FUNCTION_MAX_DIFFERENCE * pl_parallel_projected.size() * heuristic_sensitivity_ratio)
		{
			return false;
		}
	}

	//All test successfull, sample may be perspective
	return true;
}


template <typename T>
void CartAnalysis::redLon(const Container <Point3DGeographic <T> *> & pl_source, const T lon0, Container <Point3DGeographic <T> *> & pl_destination)
{
	//Corect lon0: set new central meridian for all items of the destination container
	for (unsigned int i = 0; i < pl_source.size(); i++)
	{
		//Create copy of a point
		Point3DGeographic <T> *point = pl_source[i]->clone();

		//Set new longitude
		point->setLon(CartTransformation::redLon0(point->getLon(), lon0));

		//Add point with modified longitude to the list
		pl_destination.push_back(point);
	}
}


template <typename T>
void CartAnalysis::redLon(const Container <Point3DGeographic <T> *> & pl_source, const T lon0)
{
	//Corect lon0: set new central meridian for all items of the destination container
	for (unsigned int i = 0; i < pl_source.size(); i++)
	{
		//Create copy of a point
		pl_source[i]->setLon(CartTransformation::redLon0(pl_source[i]->getLon(), lon0));
	}
}



template <typename T>
T CartAnalysis::getInitialRadius(const Projection <T> *proj, const T scale, const Container <Point3DGeographic <T> *> & pl_reference)
{
	//Get initial radius of the of the  auxilliary sphere in the analyzed map
	T sum_lat = 0.0, sum_lon = 0.0, h_max = 0.0, k_max = 0.0;

	const Point3DGeographic<T> pole = proj->getCartPole();

	for (unsigned int i = 0; i < pl_reference.size(); i++)
	{
		//Get type of the direction
		TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

		//Reduce lon
		const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), proj->getLon0());

		//Convert geographic point to oblique aspect
		const T lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon());
		const T lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, lat_trans, proj->getCartPole().getLat(), proj->getCartPole().getLon(), trans_lon_dir);

		//Distortions
		Point3DGeographic <T> point_temp(lat_trans, lon_trans);
		const T h = CartDistortion::H(NUM_DERIV_STEP, &point_temp, proj, false);
		const T k = CartDistortion::K(NUM_DERIV_STEP, &point_temp, proj, false);

		if (h > h_max)
			h_max = h;

		if (k > k_max)
			k_max = k;
	}

	//Determine the initial radius
	const T R0 = proj->getR() / scale /* * sqrt( h_max * k_max )*/;

	return R0;
}


template <typename T>
void CartAnalysis::removeSingularPoints(const Container <Node3DCartesian <T> *> & nl_source, const Container <Point3DGeographic <T> *> & pl_source, const Projection <T>  *proj, Container <Node3DCartesian <T> *> &nl_destination, Container <Point3DGeographic <T> *> &pl_destination, typename TDevIndexPairs <T>::Type & non_singular_point_pairs)
{
	//Remove all singular points from computations
	for (unsigned int i = 0; i < nl_source.size(); i++)
	{
		if ((*pl_source[i] != proj->getCartPole()) && ((proj->getCartPole().getLon() >= 0) && (pl_source[i]->getLon() != proj->getCartPole().getLon() - 180.0) ||
			(proj->getCartPole().getLon() < 0) && (pl_source[i]->getLon() != proj->getCartPole().getLon() + 180.0)))
		{
			nl_destination.push_back(nl_source[i]->clone());
			pl_destination.push_back(pl_source[i]->clone());

			non_singular_point_pairs.push_back(std::make_pair(i + 1.0, i));
		}

		//else std::cout << i << '\n';
	}
}


template <typename T>
void CartAnalysis::removeSingularPoints(const Container <Node3DCartesian <T> *> &nl_source, const Container <Point3DGeographic <T> *> &pl_source, const Projection <T>  *proj, typename TDevIndexPairs <T>::Type & non_singular_point_pairs)
{
	//Remove all singular points from computations
	typename TItemsList <Node3DCartesian <T> *> ::Type::iterator i_nl = nl_source.begin();
	typename TItemsList <Point3DGeographic <T> *> ::Type::iterator i_pl = pl_source.begin();

	for (unsigned int i = 0; (i_nl != nl_source.end()) && (i_pl != pl_source.end()); i++)
	{
		if ((*pl_source[i] != proj->getCartPole()) && ((proj->getCartPole().getLon() >= 0) && (pl_source[i]->getLon() != proj->getCartPole().getLon() - 180.0) ||
			(proj->getCartPole().getLon() < 0) && (pl_source[i]->getLon() != proj->getCartPole().getLon() + 180.0)))
		{
			nl_source->erase(i_nl++); pl_source->erase(i_pl++);
			non_singular_point_pairs.push_back(std::make_pair(i + 1.0, i));
		}
		else
		{
			++i_nl; ++i_pl;
		}
	}
}


template <typename T>
void CartAnalysis::correctMeridiansAndParrallels(typename TMeridiansList <T> ::Type & meridians, typename TParallelsList <T> ::Type & parallels,
	typename TDevIndexPairs<T>::Type & point_pairs)
{
	//Process meridians and paralles: remove inappropriate points detected as outliers using IRLS

	//Process all meridians
	for (typename TMeridiansList <T> ::Type::iterator i_meridians = meridians.begin(); i_meridians != meridians.end();)
	{
		//Get points of the meridian: break invariant
		TIndexList &meridian_points = const_cast <TIndexList&> ((*i_meridians).getPointsIndices());

		//Compare meridian points with list of pairs: remove indices of points missing in the list of pairs
		meridian_points.erase(std::remove_if(meridian_points.begin(), meridian_points.end(),
			removeUnequalMeridianParallelPointIndices <T>(point_pairs)), meridian_points.end());

		//Transform all indices of k-best points to new index list
		std::transform(meridian_points.begin(), meridian_points.end(), meridian_points.begin(), findMeridianParallelPointIndices <T>(point_pairs));

		//Not enough points, erase meridian from set
		if (meridian_points.size() < RANSAC_MIN_LINE_POINTS)
		{
			meridians.erase(i_meridians++);
			continue;
		}
		else ++i_meridians;
	}

	//Process all parallels
	for (typename TParallelsList <T> ::Type::iterator i_parallels = parallels.begin(); i_parallels != parallels.end();)
	{
		//Get points of the parallel: break invariant
		TIndexList &parallel_points = const_cast <TIndexList&> ((*i_parallels).getPointsIndices());

		//Compare parallel points with list of pairs: remove indices of points missing in the list of pairs
		parallel_points.erase(std::remove_if(parallel_points.begin(), parallel_points.end(),
			removeUnequalMeridianParallelPointIndices <T>(point_pairs)), parallel_points.end());

		//Transform all indices of k-best points to new index list
		std::transform(parallel_points.begin(), parallel_points.end(), parallel_points.begin(), findMeridianParallelPointIndices <T>(point_pairs));

		//Not enough points, erase parallel from set
		if (parallel_points.size() < RANSAC_MIN_LINE_POINTS)
		{
			parallels.erase(i_parallels++);
			continue;
		}
		else ++i_parallels;
	}
}



template <typename T>
void CartAnalysis::correctPointsMeridiansAndParrallels(const Container <Node3DCartesian <T> *> &nl_test_corr, const Container <Point3DGeographic <T> *> &pl_reference_corr, const typename TMeridiansList <T> ::Type & meridians, const typename TParallelsList <T> ::Type & parallels, const unsigned int n,
	Container <Node3DCartesian <T> *> **p_nl_test, Container <Point3DGeographic <T> *> **p_pl_reference, typename TMeridiansList <T> ::Type ** p_meridians, typename TParallelsList <T> ::Type ** p_parallels, typename TMeridiansList <T> ::Type & meridians_corr,
	typename TParallelsList <T> ::Type & parallels_corr, typename TDevIndexPairs<T>::Type & point_pairs_corr, bool &uncorrect_points_found)
{
	//Process points, meridians and paralles: remove inappropriate points, correct meridians and parallels
	const unsigned int n_corr = nl_test_corr.size();

	//Correction of meridians, parallels and points must be done
	if (n != n_corr)
	{
		//Set flag to true, used for all samples using non-singular sets
		uncorrect_points_found = true;

		//Create copy of meridians / parallels
		meridians_corr = meridians; parallels_corr = parallels;

		//Process all meridians
		for (typename TMeridiansList <T> ::Type::iterator i_meridians = meridians.begin(); i_meridians != meridians.end();)
		{
			//Get points of the meridian: break invariant
			TIndexList &meridian_points = const_cast <TIndexList&> (i_meridians->getPointsIndices());

			//Compare meridian points with list of pairs: remove indices of points missing in the list of pairs
			meridian_points.erase(std::remove_if(meridian_points.begin(), meridian_points.end(),
				removeUnequalMeridianParallelPointIndices <T>(point_pairs_corr)), meridian_points.end());

			//Transform all indices of k-best points to new index list
			std::transform(meridian_points.begin(), meridian_points.end(), meridian_points.begin(), findMeridianParallelPointIndices <T>(point_pairs_corr));

			//Not enough points, erase meridian from set
			if (meridian_points.size() < RANSAC_MIN_LINE_POINTS)
			{
				meridians.erase(i_meridians++);
				continue;
			}
			else ++i_meridians;
		}

		//Process all parallels
		for (typename TParallelsList <T> ::Type::iterator i_parallels = parallels.begin(); i_parallels != parallels.end();)
		{
			//Get points of the parallel: break invariant
			TIndexList &parallel_points = const_cast <TIndexList&> (i_parallels->getPointsIndices());

			//Compare parallel points with list of pairs: remove indices of points missing in the list of pairs
			parallel_points.erase(std::remove_if(parallel_points.begin(), parallel_points.end(),
				removeUnequalMeridianParallelPointIndices <T>(point_pairs_corr)), parallel_points.end());

			//Transform all indices of k-best points to new index list
			std::transform(parallel_points.begin(), parallel_points.end(), parallel_points.begin(), findMeridianParallelPointIndices <T>(point_pairs_corr));

			//Not enough points, erase parallel from set
			if (parallel_points.size() < RANSAC_MIN_LINE_POINTS)
			{
				parallels.erase(i_parallels++);
				continue;
			}
			else ++i_parallels;
		}

		//Convert new pairs to index list: indices will be printed to the log file
		std::transform(point_pairs_corr.begin(), point_pairs_corr.end(), std::back_inserter(point_pairs_corr), getSecondElementInPair());

		//Set pointers to corrected meridians and parallels
		p_meridians = &meridians_corr; p_parallels = &parallels_corr;

		//Set pointers to corrected sets
		p_nl_test = &nl_test_corr; p_pl_reference = &pl_reference_corr;
	}

	//No correction neccessary
	else
	{
		//Set flag to true, used for all samples using non-singular sets
		uncorrect_points_found = false;

		//Initialize point indices: all indices are valid
		for (unsigned int j = 0; j < n_corr; j++) point_pairs_corr.push_back(j);
	}
}


//Each sample analysis
template <typename T>
void CartAnalysis::analyzeSampleCrossNearestNeighbourDistance(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const float mult_ratio)
{
	//Analyze all samples using cross nearest distance ratio
	try
	{
		const unsigned int n = nl_test.size();

		//Compute center of mass
		T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

		for (unsigned int i = 0; i < n; i++)
		{
			x_mass_test += nl_test[i]->getX();
			y_mass_test += nl_test[i]->getY();

			x_mass_reference += nl_projected[i]->getX();
			y_mass_reference += nl_projected[i]->getY();
		}

		//Centers of mass
		x_mass_test = x_mass_test / n;
		y_mass_test = y_mass_test / n;
		x_mass_reference = x_mass_reference / n;
		y_mass_reference = y_mass_reference / n;

		//Reduce coordinates
		Container <Node3DCartesian <T> *> nl_test_red;
		Container <Node3DCartesianProjected <T> *> nl_projected_red;

		for (unsigned int i = 0; i < n; i++)
		{
			Node3DCartesian <T> *point_test = new Node3DCartesianProjected <T>(nl_test[i]->getX() - x_mass_test, nl_test[i]->getY() - y_mass_test);
			Node3DCartesianProjected <T> *point_proj = new Node3DCartesianProjected <T>(nl_projected[i]->getX() - x_mass_reference, nl_projected[i]->getY() - y_mass_reference);

			nl_test_red.push_back(point_test);
			nl_projected_red.push_back(point_proj);
		}

		//Compute cross nearest distance ratio
		const T cross_nearest_neighbour_distance_ratio = NNDistance::getCrossNearestNeighbourDistance(nl_test_red, nl_projected_red);

		s.setCrossNearestNeighbourDistanceRatio(cross_nearest_neighbour_distance_ratio * mult_ratio);
		s.setCrossNearestNeighbourDistanceRatioPosition(1);
	}

	//Throw exception
	catch (Exception & error)
	{
		s.setCrossNearestNeighbourDistanceRatio(MAX_FLOAT);
		s.setCrossNearestNeighbourDistanceRatioPosition(-1);
	}
}



template <typename T>
void CartAnalysis::analyzeSampleHomotheticTransformationDeviation(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const TMatchPointsType match_type, const float mult_ratio)
{
	//Analyze sample using Homothetic transformation deviation
	try
	{
		//List of transformed points
		Container <Node3DCartesian <T> *> nl_transformed;

		//Compute homothetic transformation
		TTransformationKeyHomothetic2D <T> key_homothetic;
		HomotheticTransformation2D::transformPoints(nl_projected, nl_test, nl_transformed, key_homothetic);

		//Get transformation charactersitic
		TAccuracyCharacteristics <T> deviations = Transformation2D::getAccuracyCharacteristics(nl_projected, nl_test, nl_transformed, key_homothetic);
		T homothetic_transformation_ratio = deviations.std_dev;

		//Compute ratio and percentage match using uncertainty regions
		const unsigned int n = nl_projected.size();
		TIndexList matched_points;

		//Get extreme coordinates
		const T x_min = (*(std::min_element(nl_projected.begin(), nl_projected.end(), sortPointsByX())))->getX();
		const T y_min = (*(std::min_element(nl_projected.begin(), nl_projected.end(), sortPointsByY())))->getY();
		const T x_max = (*(std::max_element(nl_projected.begin(), nl_projected.end(), sortPointsByX())))->getX();
		const T y_max = (*(std::max_element(nl_projected.begin(), nl_projected.end(), sortPointsByY())))->getY();

		//Compute radius of the circle: modified equation given by Wamelen et al (2001)
		const T max_radius = 0.25 * std::max(x_max - x_min, y_max - y_min) / sqrt((float)n);

		//Compute matching
		T homothetic_transformation_perc_match = (match_type == MatchCircle ? getMatchRatioCircle(nl_projected, nl_transformed, matched_points, CollectOn, max_radius, 0.1) :
			getMatchRatioTissotIndicatrix(nl_projected, nl_transformed, matched_points, CollectOn, max_radius, 0.5));


		//Set ratio, percentage match, scale, tx and ty for the sample
		s.setHomotheticTransformationRatio(mult_ratio * homothetic_transformation_ratio);
		s.setHomotheticTransformationRatioPosition(1);
		s.setHomotheticTransformationPercMatch((unsigned int)homothetic_transformation_perc_match);
		s.setHomotheticTransformationMatchedPointsIndices(matched_points);
		s.setScaleHomT(key_homothetic.c);
		//s.setDx ( key_homothetic.x_mass_local - key_homothetic.x_mass_global / key_homothetic.c );
		//s.setDy ( key_homothetic.y_mass_local - key_homothetic.y_mass_global / key_homothetic.c );
	}

	//Throw exception
	catch (Exception & error)
	{
		s.setHomotheticTransformationRatio(MAX_FLOAT);
		s.setHomotheticTransformationPercMatch(0.0);
		s.setHomotheticTransformationRatioPosition(-1);
	}
}


template <typename T>
void CartAnalysis::analyzeSampleHelmertTransformationDeviation(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const TMatchPointsType match_type, const float mult_ratio)
{
	//Analyze sample using Helmert transformation deviation
	try
	{
		//List of transformed points
		Container <Node3DCartesian <T> *> nl_transformed;

		//Compute Helmert transformation
		TTransformationKeyHelmert2D <T> key_helmert;
		HelmertTransformation2D::transformPoints(nl_projected, nl_test, nl_transformed, key_helmert);
		//HelmertTransformation2D::transformPoints (nl_test, nl_projected, nl_transformed, key_helmert);
		//nl_transformed.print();

		//Compute ratio and percentage match using uncertainty regions
		TAccuracyCharacteristics <T> deviations = Transformation2D::getAccuracyCharacteristics(nl_projected, nl_test, nl_transformed, key_helmert);
		T helmert_transformation_ratio = deviations.std_dev;

		//Compute ratio and percentage match using uncertainty regions
		TIndexList matched_points;
		const unsigned int n = nl_projected.size();

		//Get extreme coordinates
		const T x_min = (*(std::min_element(nl_projected.begin(), nl_projected.end(), sortPointsByX())))->getX();
		const T y_min = (*(std::min_element(nl_projected.begin(), nl_projected.end(), sortPointsByY())))->getY();
		const T x_max = (*(std::max_element(nl_projected.begin(), nl_projected.end(), sortPointsByX())))->getX();
		const T y_max = (*(std::max_element(nl_projected.begin(), nl_projected.end(), sortPointsByY())))->getY();

		//Compute radius of the circle: modified equation given by Wamelen et al (2001)
		const T max_radius = 0.25 * std::max(x_max - x_min, y_max - y_min) / sqrt((float)n);

		//Compute matching
		T helmert_transformation_perc_match = (match_type == MatchCircle ? getMatchRatioCircle(nl_projected, nl_transformed, matched_points, CollectOn, max_radius, 0.1) :
			getMatchRatioTissotIndicatrix(nl_projected, nl_transformed, matched_points, CollectOn, max_radius, 0.5));

		//std::cout <<"match" << helmert_transformation_perc_match;

		//Set ratio and percentage match for the sample
		s.setHelmertTransformationRatio(mult_ratio * helmert_transformation_ratio);
		s.setHelmertTransformationRatioPosition(1);
		s.setHelmertTransformationPercMatch((unsigned int)helmert_transformation_perc_match);
		s.setHelmertTransformationMatchedPointsIndices(matched_points);
		s.setScaleHelT(sqrt(key_helmert.c1 * key_helmert.c1 + key_helmert.c2 * key_helmert.c2));
		s.setRotation(atan2(key_helmert.c2, key_helmert.c1) * 180.0 / M_PI);
		//s.setDx ( key_helmert.x_mass_local - key_helmert.x_mass_global / sqrt ( key_helmert.c1 * key_helmert.c1 + key_helmert.c2 * key_helmert.c2 ) );
		//s.setDy ( key_helmert.y_mass_local - key_helmert.y_mass_global / sqrt ( key_helmert.c1 * key_helmert.c1 + key_helmert.c2 * key_helmert.c2 ) );
	}

	//Throw exception
	catch (Exception & error)
	{
		s.setHelmertTransformationRatio(MAX_FLOAT);
		s.setHelmertTransformationPercMatch(0.0);
		s.setHelmertTransformationRatioPosition(-1);
	}
}


template <typename T>
void CartAnalysis::analyzeSampleGeographicNetworkTurningFunctionRatio(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected,
	const typename TMeridiansList <T> ::Type &meridians, const typename TParallelsList <T> ::Type &parallels, const float mult_ratio)
{
	//Analyze sample using angular turning function differences in geographic network
	T turning_function_ratio_meridians = 0, turning_function_ratio_parallels = 0;

	try
	{
		//No meridian or parallel
		if ((meridians.size() == 0) && (parallels.size() == 0))
		{
			throw BadDataException("BadDataException: no meridians and parallels. ", "Can not perform analysis of turning function.");
		}

		typename TMeridiansList <T> ::Type::const_iterator i_meridians = meridians.begin();

		//Process all meridians
		for (i_meridians = meridians.begin(); i_meridians != meridians.end(); ++i_meridians)
		{
			TIndexList il = (*i_meridians).getPointsIndices();

			//Convert test meridian to Points 2D list
			Container <Point3DCartesian <T> > pl_meridian_test(nl_test, (*i_meridians).getPointsIndices());

			//Convert projected meridian to Points2D list
			Container <Point3DCartesian <T> > pl_meridian_projected(nl_projected, (*i_meridians).getPointsIndices());

			//Compute turning function difference for each test and projected meridian
			turning_function_ratio_meridians += TurningFunction::compare2PolyLinesUsingTurningFunction(pl_meridian_test, pl_meridian_projected, RotationDependent, ScaleInvariant);
		}

		typename TParallelsList <T> ::Type::const_iterator i_parallels = parallels.begin();

		//Process all parallels
		for (i_parallels = parallels.begin(); i_parallels != parallels.end(); ++i_parallels)
		{
			//Convert test meridian to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_test(nl_test, (*i_parallels).getPointsIndices());

			//Convert projected meridian to Points2D list
			Container <Point3DCartesian <T> > pl_parallel_projected(nl_projected, (*i_parallels).getPointsIndices());

			//Compute turning function difference for each parallel
			turning_function_ratio_parallels += TurningFunction::compare2PolyLinesUsingTurningFunction(pl_parallel_test, pl_parallel_projected, RotationDependent, ScaleInvariant);
		}

		//Set turning function ratio
		s.setGNTurningFunctionRatio(mult_ratio * (turning_function_ratio_meridians + turning_function_ratio_parallels));
		s.setGNTurningFunctionRatioPosition(1);
	}

	//Throw exception
	catch (Exception & error)
	{
		s.setGNTurningFunctionRatio(MAX_FLOAT);
		s.setGNTurningFunctionRatioPosition(-1);
	}
}


/*
template <typename T, TDestructable destructable>
void CartAnalysis::analyzeSampleUsingVoronoiDiagramTurningFunctionRatio(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const Container <Face <T> *, destructable> &faces_test,
const TAnalysisParameters <T> & analysis_parameters, const float mult_ratio)
{
//Analyze sample using Voronoi diagram ratios
try
{
const unsigned int n_test_points = nl_test.size();

//Compute Voronoi diagram for test dataset after outliers removed
Container <HalfEdge <double> *> hl_dt_test, hl_vor_test, hl_merge_test;
Container <Node3DCartesian <double> *> nl_vor_test, intersections_test;
Container <VoronoiCell <double> *> vor_cells_test;

//Voronoi diagram of the reference dataset, data structures
Container <HalfEdge <T> *> hl_dt_reference, hl_vor_reference, hl_merge_reference;
Container <Node3DCartesian <T> *> nl_vor_reference, intersections_reference;
Container <VoronoiCell <T> *> vor_cells_list_reference;

//Create Voronoi diagram for the reference dataset in any case
Voronoi2D::VD((Container <Node3DCartesian <T> *> &) nl_projected, nl_vor_reference, hl_dt_reference, hl_vor_reference, vor_cells_list_reference, AppropriateBoundedCells, TopologicApproach, 0, analysis_parameters.print_exceptions);

//Modified dataset with found outliers, rotated dataset or dataset with found singular points: create new Voronoi diagram for the test dataset
//Do not use pre-generated test faces
if (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())
Voronoi2D::VD((Container <Node3DCartesian <double> *> &) nl_test, nl_vor_test, hl_dt_test, hl_vor_test, vor_cells_test, AppropriateBoundedCells, TopologicApproach, 0, analysis_parameters.print_exceptions);

//Get total unbounded cells for the test dataset
unsigned int total_bounded_pairs_of_cell = 0;

//Initialize ratio
T turning_function_difference = 0.0;

//Process all merged faces
for (unsigned int index_faces = 0; index_faces < n_test_points; index_faces++)
{
//An appropriate face in the non-modified test dataset exists
//For modified  datasets with found outliers, rotated or found singularitities: Construct a new face for test dataset
if ((s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound()) || faces_test[index_faces] != NULL)
{
//Get Voronoi cell of the reference dataset
VoronoiCell <T> *vor_cell_reference = dynamic_cast <VoronoiCell <T> *> (nl_projected[index_faces]->getFace());

//Get Voronoi cells of the modified test dataset from newly created Voronoi diagram: only if outliers found, rotated dataset or singular points found
VoronoiCell <T> *vor_cell_test = (s.getOutliersFound() || s.getRotatedSample() || s.getOutliersFound() ? dynamic_cast <VoronoiCell <T> *> (nl_test[index_faces]->getFace()) : NULL);

//Remove_poutliers = 0, rotated dataset = 0 and remove singular points = 0: bounded Voronoi cell of the reference dataset exists and corresponding face of the test dataset exists
//Remove_poutliers = 1, rotated dataset = 1 and remove singulat points = 1: bounded Voronoi cell of the reference datasets exist and corresponding bounded Voronoi cell exist
if ((vor_cell_reference != NULL) && (vor_cell_reference->getBounded()) && ((!s.getOutliersFound() && !s.getRotatedSample() && !s.getSingularPointsFound()) ||
(vor_cell_test != NULL) && (vor_cell_test->getBounded()) && (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())))
{
//Pointer to merged reference face
Face <T> * face_reference = NULL, *face_test = NULL;

try
{
//Merge reference Voronoi cell with adjacent cells
Voronoi2D::mergeVoronoiCellAndAdjacentCells(vor_cell_reference, &face_reference, intersections_reference, hl_merge_reference);

//Merge test Voronoi cell with adjacent cells: only if outliers found, dataset is rotated or singular points found
if (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())
Voronoi2D::mergeVoronoiCellAndAdjacentCells(vor_cell_test, &face_test, intersections_test, hl_merge_test);

//Count unbounded pairs of Voronoi cells
total_bounded_pairs_of_cell++;

//Compute turning function difference for both merged faces: for modified dataset use newly created face, otherwise use a precomputed face
turning_function_difference += s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound() ? TurningFunction::compare2FacesUsingTurningFunction(face_test, face_reference, RotationDependent, ScaleInvariant) :
TurningFunction::compare2FacesUsingTurningFunction(faces_test[index_faces], face_reference, RotationDependent, ScaleInvariant);

//Delete merged face
if (face_reference != NULL) delete face_reference;

if (face_test != NULL) delete face_test;
}

//Throw exception
catch (Exception & error)
{
if (face_reference != NULL) delete face_reference;

if (face_test != NULL) delete face_test;

throw error;
}
}
}
}

//Not enough corresponding pairs of bounded Voronoi cells suitable for analysis
if (total_bounded_pairs_of_cell < MIN_BOUNDED_VORONOI_CELLS)
{
throw BadDataException("BadDataException: not enough unbounded pairs, ", "set values");
}

//Set results of the analysis
s.setVoronoiCellTurningFunctionRatio(mult_ratio * sqrt(turning_function_difference / total_bounded_pairs_of_cell));
s.setVoronoiCellTurningFunctionRatioPosition(1);
}

//Error while sample processed or invalid configuration for cell by cell ratia
catch (Exception & error)
{
//Do not compute other analysis
if (analysis_parameters.analysis_type.a_vd_tf)
{
s.setVoronoiCellTurningFunctionRatio(MAX_FLOAT);
s.setVoronoiCellTurningFunctionRatioPosition(-1);
}
}
}
*/

template <typename T, TDestructable destructable>
void CartAnalysis::analyzeSampleUsingVoronoiDiagramTurningFunctionRatio2(Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const Container <Face <T> *, destructable> &faces_test,
	const TAnalysisParameters <T> & analysis_parameters, const float mult_ratio)
{
	//Analyze sample using Voronoi diagram ratios
	try
	{
		const unsigned int n_test_points = nl_test.size();

		//Compute Voronoi diagram for test dataset after outliers removed
		Container <HalfEdge <double> *> hl_dt_test, hl_vor_test, hl_merge_test;
		Container <Node3DCartesian <double> *> nl_vor_test, intersections_test;
		Container <VoronoiCell <double> *> vor_cells_test;

		//Voronoi diagram of the reference dataset, data structures
		Container <HalfEdge <T> *> hl_dt_reference, hl_vor_reference, hl_merge_reference;
		Container <Node3DCartesian <T> *> nl_vor_reference, intersections_reference;
		Container <VoronoiCell <T> *> vor_cells_list_reference;

		//Create Voronoi diagram for the reference dataset in any case
		Voronoi2D::VD((Container <Node3DCartesian <T> *> &) nl_projected, nl_vor_reference, hl_dt_reference, hl_vor_reference, vor_cells_list_reference, AppropriateBoundedCells, TopologicApproach, 0, analysis_parameters.print_exceptions);

		//Modified dataset with found outliers, rotated dataset or dataset with found singular points: create new Voronoi diagram for the test dataset
		//Do not use pre-generated test faces
		if (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())
			Voronoi2D::VD((Container <Node3DCartesian <double> *> &) nl_test, nl_vor_test, hl_dt_test, hl_vor_test, vor_cells_test, AppropriateBoundedCells, TopologicApproach, 0, analysis_parameters.print_exceptions);

		const unsigned int n_faces = faces_test.size();

		//Get total unbounded cells for the test dataset
		unsigned int total_bounded_pairs_of_cell = 0;

		//Initialize ratio
		T turning_function_difference = 0.0, inner_distance_difference = 0.0;

		//Process all merged faces
		for (unsigned int index_faces = 0; index_faces < n_test_points; index_faces++)
		{
			//An appropriate face in the non-modified test dataset exists
			//For modified  datasets with found outliers, rotated or found singularitities: Construct a new face for test dataset
			if ((s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound()) || ((n_faces > 0) && (faces_test[index_faces] != NULL)))
			{
				//Get Voronoi cell of the reference dataset
				VoronoiCell <T> *vor_cell_reference = dynamic_cast <VoronoiCell <T> *> (nl_projected[index_faces]->getFace());

				//Get Voronoi cells of the modified test dataset from newly created Voronoi diagram: only if outliers found, rotated dataset or singular points found
				VoronoiCell <T> *vor_cell_test = (s.getOutliersFound() || s.getRotatedSample() || s.getOutliersFound() ? dynamic_cast <VoronoiCell <T> *> (nl_test[index_faces]->getFace()) : NULL);

				//Remove_poutliers = 0, rotated dataset = 0 and remove singular points = 0: bounded Voronoi cell of the reference dataset exists and corresponding face of the test dataset exists
				//Remove_poutliers = 1, rotated dataset = 1 and remove singulat points = 1: bounded Voronoi cell of the reference datasets exist and corresponding bounded Voronoi cell exist
				if ((vor_cell_reference != NULL) && (vor_cell_reference->getBounded()) && ((!s.getOutliersFound() && !s.getRotatedSample() && !s.getSingularPointsFound()) ||
					(vor_cell_test != NULL) && (vor_cell_test->getBounded()) && (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())))
				{
					//Pointer to merged reference face
					Face <T> * face_reference = NULL, *face_test = NULL;

					try
					{
						//Merge reference Voronoi cell with adjacent cells
						Voronoi2D::mergeVoronoiCellAndAdjacentCells(vor_cell_reference, &face_reference, intersections_reference, hl_merge_reference);

						//Merge test Voronoi cell with adjacent cells: only if outliers found, dataset is rotated or singular points found
						if (s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound())
							Voronoi2D::mergeVoronoiCellAndAdjacentCells(vor_cell_test, &face_test, intersections_test, hl_merge_test);

						//Count unbounded pairs of Voronoi cells
						total_bounded_pairs_of_cell++;

						//Compute turning function difference for both merged faces: for modified dataset use newly created face, otherwise use a precomputed face
						turning_function_difference += s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound() ? TurningFunction::compare2FacesUsingTurningFunction(face_test, face_reference, RotationDependent, ScaleInvariant) :
							TurningFunction::compare2FacesUsingTurningFunction(faces_test[index_faces], face_reference, RotationDependent, ScaleInvariant);

						//Compute inner distance difference
						inner_distance_difference += s.getOutliersFound() || s.getRotatedSample() || s.getSingularPointsFound() ? InnerDistance::compare2FacesUsingInnerDistances(face_test, face_reference) :
							InnerDistance::compare2FacesUsingInnerDistances(faces_test[index_faces], face_reference);
						//inner_distance_difference += InnerDistance::compare2FacesUsingInnerDistances(faces_test[index_faces], face_reference);

						//Delete merged face
						if (face_reference != NULL) delete face_reference;

						if (face_test != NULL) delete face_test;
					}

					//Throw exception
					catch (Exception & error)
					{
						if (face_reference != NULL) delete face_reference;

						if (face_test != NULL) delete face_test;

						throw error;
					}
				}
			}
		}

		//Not enough corresponding pairs of bounded Voronoi cells suitable for analysis
		if (total_bounded_pairs_of_cell < MIN_BOUNDED_VORONOI_CELLS)
		{
			throw BadDataException("BadDataException: not enough unbounded pairs, ", "set values");
		}

		//Set results of the analysis
		s.setVoronoiCellTurningFunctionRatio(mult_ratio * sqrt(turning_function_difference / total_bounded_pairs_of_cell));
		s.setVoronoiCellInnerDistanceRatio(sqrt(inner_distance_difference / total_bounded_pairs_of_cell));
		s.setVoronoiCellTurningFunctionRatioPosition(1);
		s.setVoronoiCellInnerDistanceRatioPosition(1);
	}

	//Error while sample processed or invalid configuration for cell by cell ratia
	catch (Exception & error)
	{
		//faces_reference.clear();

		//Do not compute other analysis
		if (analysis_parameters.analysis_type.a_vd_tf)
		{
			s.setVoronoiCellTurningFunctionRatio(MAX_FLOAT);
			s.setVoronoiCellTurningFunctionRatioPosition(-1);
			s.setVoronoiCellInnerDistanceRatio(MAX_FLOAT);
			s.setVoronoiCellInnerDistanceRatioPosition(-1);
		}
	}
}


template <typename T>
void CartAnalysis::sortSamplesByComputedRatios(Container <Sample <T> > &sl, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type)
{

	//Sort results: Cross nearest distance
	if (analysis_type.a_cnd)
		std::sort(sl.begin(), sl.end(), sortSamplesByCrossNearestNeighbourDistanceRatio());

	typename TAnalysisParameters <T>::TAnalysisType a1(analysis_type.a_cnd, 0, 0, 0, 0, 0);
	setPositionForSortedSamples(sl, a1);

	//Sort results: Homothetic transformation, standard deviation + match function
	if (analysis_type.a_homt)
		std::sort(sl.begin(), sl.end(), sortSamplesByHomotheticTransformationRatio());

	typename TAnalysisParameters <T>::TAnalysisType a2(0, analysis_type.a_homt, 0, 0, 0, 0);
	setPositionForSortedSamples(sl, a2);

	//Sort results: Helmert transformation, standard deviation + match function
	if (analysis_type.a_helt)
		std::sort(sl.begin(), sl.end(), sortSamplesByHelmertTransformationRatio());

	typename TAnalysisParameters <T>::TAnalysisType a3(0, 0, analysis_type.a_helt, 0, 0, 0);
	setPositionForSortedSamples(sl, a3);

	//Sort results: geographic network turning function difference ratio
	if (analysis_type.a_gn_tf)
		std::sort(sl.begin(), sl.end(), sortSamplesByGNTurningFunctionRatio());

	typename TAnalysisParameters <T>::TAnalysisType a4(0, 0, 0, analysis_type.a_gn_tf, 0, 0);
	setPositionForSortedSamples(sl, a4);

	//Sort results: Voronoi cell Turning Function Ratio
	if (analysis_type.a_vd_tf)
		std::sort(sl.begin(), sl.end(), sortSamplesByVoronoiCellTurningFunctionRatio());

	typename TAnalysisParameters <T>::TAnalysisType a5(0, 0, 0, 0, analysis_type.a_vd_tf, 0);
	setPositionForSortedSamples(sl, a5);

	// Sort results : Voronoi cell Inner distance Ratio
	if (analysis_type.a_vd_id)
		std::sort(sl.begin(), sl.end(), sortSamplesByVoronoiCellInnerDistanceRatio());

	typename TAnalysisParameters <T>::TAnalysisType a6(0, 0, 0, 0, 0, analysis_type.a_vd_id);
	setPositionForSortedSamples(sl, a6);

	//Sort samples by all ratios: result = arithmetic mean
	std::sort(sl.begin(), sl.end(), sortSamplesByAllRatios <T>(analysis_type));
}


template <typename T>
void CartAnalysis::setPositionForSortedSamples(Container <Sample <T> > &sl, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type)
{
	//Set position for each sample after its sorting
	unsigned int n = sl.size();

	for (unsigned int i = 1; i < n; i++)
	{
		//Cross nearest distance criterion
		if (analysis_type.a_cnd)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getCrossNearestNeighbourDistanceRatio() -
				sl[i - 1].getCrossNearestNeighbourDistanceRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: first position
				if (sl[i - 1].getCrossNearestNeighbourDistanceRatioPosition() < 0)
					sl[i].setCrossNearestNeighbourDistanceRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getCrossNearestNeighbourDistanceRatioPosition() > 0)
				{
					if (sl[i].getCrossNearestNeighbourDistanceRatioPosition() > 0)
						sl[i].setCrossNearestNeighbourDistanceRatioPosition(sl[i - 1].getCrossNearestNeighbourDistanceRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setCrossNearestNeighbourDistanceRatioPosition(sl[i - 1].getCrossNearestNeighbourDistanceRatioPosition());
			}
		}

		//Homothetic transformation, standard deviation criterion + match factor
		if (analysis_type.a_homt)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getHomotheticTransformationRatio() -
				sl[i - 1].getHomotheticTransformationRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: first position
				if (sl[i - 1].getHomotheticTransformationRatioPosition() < 0)
					sl[i].setHomotheticTransformationRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getHomotheticTransformationRatioPosition() > 0)
				{
					if (sl[i].getHomotheticTransformationRatioPosition() > 0)
						sl[i].setHomotheticTransformationRatioPosition(sl[i - 1].getHomotheticTransformationRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setHomotheticTransformationRatioPosition(sl[i - 1].getHomotheticTransformationRatioPosition());
			}
		}

		//Helmert transformation, standard deviation criterion + match factor
		if (analysis_type.a_helt)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getHelmertTransformationRatio() -
				sl[i - 1].getHelmertTransformationRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: first position
				if (sl[i - 1].getHelmertTransformationRatioPosition() < 0)
					sl[i].setHelmertTransformationRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getHelmertTransformationRatioPosition() > 0)
				{
					if (sl[i].getHelmertTransformationRatioPosition() > 0)
						sl[i].setHelmertTransformationRatioPosition(sl[i - 1].getHelmertTransformationRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setHelmertTransformationRatioPosition(sl[i - 1].getHelmertTransformationRatioPosition());
			}
		}

		//Turning function criterion for meridians and parallels
		else if (analysis_type.a_gn_tf)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getGNTurningFunctionRatio() -
				sl[i - 1].getGNTurningFunctionRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: first position
				if (sl[i - 1].getGNTurningFunctionRatioPosition() < 0)
					sl[i].setGNTurningFunctionRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getGNTurningFunctionRatioPosition() > 0)
				{
					if (sl[i].getGNTurningFunctionRatioPosition() > 0)
						sl[i].setGNTurningFunctionRatioPosition(sl[i - 1].getGNTurningFunctionRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setGNTurningFunctionRatioPosition(sl[i - 1].getGNTurningFunctionRatioPosition());
			}
		}

		//Voronoi cell turning function criterion
		else if (analysis_type.a_vd_tf)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getVoronoiCellTurningFunctionRatio() -
				sl[i - 1].getVoronoiCellTurningFunctionRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: position = 1
				if (sl[i - 1].getVoronoiCellTurningFunctionRatioPosition() < 0)
					sl[i].setVoronoiCellTurningFunctionRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getVoronoiCellTurningFunctionRatioPosition() > 0)
				{
					if (sl[i].getVoronoiCellTurningFunctionRatioPosition() > 0)
						sl[i].setVoronoiCellTurningFunctionRatioPosition(sl[i - 1].getVoronoiCellTurningFunctionRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setVoronoiCellTurningFunctionRatioPosition(sl[i - 1].getVoronoiCellTurningFunctionRatioPosition());
			}
		}

		//Voronoi cell inner distance criterion
		else if (analysis_type.a_vd_id)
		{
			//Actual value differs from previous value: their difference > min value
			if (fabs(sl[i].getVoronoiCellInnerDistanceRatio() -
				sl[i - 1].getVoronoiCellInnerDistanceRatio()) > ARGUMENT_ROUND_ERROR)
			{
				//Previous value negative: position = 1
				if (sl[i - 1].getVoronoiCellInnerDistanceRatioPosition() < 0)
					sl[i].setVoronoiCellInnerDistanceRatioPosition(1);

				//Previous value positive: position = position + 1
				else if (sl[i - 1].getVoronoiCellInnerDistanceRatioPosition() > 0)
				{
					if (sl[i].getVoronoiCellInnerDistanceRatioPosition() > 0)
						sl[i].setVoronoiCellInnerDistanceRatioPosition(sl[i - 1].getVoronoiCellInnerDistanceRatioPosition() + 1);
				}
			}

			//Actual value jittere as the previous value: both values are equal
			else
			{
				sl[i].setVoronoiCellInnerDistanceRatioPosition(sl[i - 1].getVoronoiCellInnerDistanceRatioPosition());
			}
		}
	}
}


template <typename T>
void CartAnalysis::printResults(const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
	const TAnalysisParameters <T> & analysis_parameters, std::ostream * output)
{

	//Print first n items sorted by the similarity match ratio
	unsigned int items_printed = analysis_parameters.printed_results;
	const unsigned int n = sl.size(), n_test = nl_test.size();

	//Correct number of printed items
	if (items_printed > n) items_printed = n;

	//Some points were loaded
	if (n > 0)
	{

		//Table  1
		*output << "Results containg values of the criteria:" << std::endl << std::endl;

		//Set properties
		*output << std::showpoint << std::fixed << std::right;

		//Create header 1 : results of analalysis
		*output << std::setw(4) << "#"
			<< std::setw(8) << "Proj"
			<< std::setw(7) << "Categ"
			<< std::setw(6) << "latP"
			<< std::setw(7) << "lonP"
			<< std::setw(6) << "lat0"
			<< std::setw(7) << "lon0"
			<< std::setw(9) << "C"
			<< std::setw(6) << "BKEY"
			<< std::setw(9) << "CND[m]"
			<< std::setw(9) << "HOMT[m]";


		//Print type of match
		analysis_parameters.match_method == MatchCircle ? *output << std::setw(5) << "+ MC" : *output << std::setw(5) << "+ MT";

		*output << std::setw(9) << "HELT[m]";

		//Print type of match
		analysis_parameters.match_method == MatchCircle ? *output << std::setw(5) << "+ MC" : *output << std::setw(5) << "+ MT";

		*output << std::setw(9) << "GNTF"
			<< std::setw(9) << "VDTF" << std::endl;

		//Values of the criterion for each projection
		for (unsigned int i = 0; i < (analysis_parameters.analyzed_projections.size() == 0 ? items_printed : n); i++)
		{
			if (analysis_parameters.analyzed_projections.size() == 0)
			{
				sl[i].printSampleRatios(i + 1, analysis_parameters.analysis_type, n_test, output);
			}

			else if (sl[i].getAnalyzedProjectionSample())
			{
				sl[i].printSampleRatios(i + 1, analysis_parameters.analysis_type, n_test, output);
			}
		}

		//Table 2
		*output << std::endl << "Results containg positions of the criteria:" << std::endl << std::endl;

		//Create header 2: positions
		*output << std::setw(4) << "#"
			<< std::setw(8) << "Proj"
			<< std::setw(7) << "Categ"
			<< std::setw(6) << "latP"
			<< std::setw(7) << "lonP"
			<< std::setw(6) << "lat0"
			<< std::setw(7) << "lon0"
			<< std::setw(9) << "C"
			<< std::setw(6) << "CND"
			<< std::setw(6) << "HOMT"
			<< std::setw(6) << "HELT"
			<< std::setw(6) << "GNTF"
			<< std::setw(6) << "VDTF"
			<< std::setw(6) << "VDID" << std::endl;

		//Table 2: Positions of the criterion for each projection
		for (unsigned int i = 0; i < (analysis_parameters.analyzed_projections.size() == 0 ? items_printed : n); i++)
		{
			if (analysis_parameters.analyzed_projections.size() == 0)
			{
				sl[i].printSamplePositions(i + 1, analysis_parameters.analysis_type, output);
			}

			else if (sl[i].getAnalyzedProjectionSample())
			{
				sl[i].printSamplePositions(i + 1, analysis_parameters.analysis_type, output);
			}
		}

		*output << std::endl << "  ( * Sample with additionaly corrected rotation, -c is enabled. )" << std::endl << std::endl;

		//Table 3
		*output << std::endl << "Analyzed and reference points:" << std::endl << std::endl;

		//Create header
		*output << std::setw(3) << "#"
			<< std::setw(15) << "X_test"
			<< std::setw(15) << "Y_test"
			<< std::setw(13) << "Fi_ref"
			<< std::setw(13) << "La_ref" << std::endl;

		//Table 3: List of analyzed and reference points
		for (unsigned int i = 0; i < n_test; i++)
		{
			//Print all analyzed points
			*output << std::setw(3) << i << std::fixed << std::setprecision(3) << std::setw(15) << nl_test[i]->getX() << std::setw(15) <<
				nl_test[i]->getY() << std::fixed << std::setprecision(5) << std::setw(13) << nl_reference[i]->getLat() <<
				std::setw(13) << nl_reference[i]->getLon() << std::endl;
		}

		*output << std::endl << "Scale, rotation and matched points for each projection:" << std::endl << std::endl;

		//Table 4: Scale, rotation, all matched points and outliers
		for (unsigned int i = 0; i < (analysis_parameters.analyzed_projections.size() == 0 ? items_printed : n); i++)
		{
			if (analysis_parameters.analyzed_projections.size() == 0)
			{
				sl[i].printSampleMatchedPoints(nl_test, nl_reference, i + 1, analysis_parameters.analysis_type, output);
			}

			else if (sl[i].getAnalyzedProjectionSample())
			{
				sl[i].printSampleMatchedPoints(nl_test, nl_reference, i + 1, analysis_parameters.analysis_type, output);
			}
		}

		*output << std::endl;
		*output << std::endl;

	}

}

template <typename T>
void CartAnalysis::printResults2(const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
	const TAnalysisParameters <T> & analysis_parameters, std::ostream * output)
{

	//Print first n items sorted by the similarity match ratio
	unsigned int items_printed = analysis_parameters.printed_results;
	const unsigned int n = sl.size(), n_test = nl_test.size();

	//Correct number of printed items
	if (items_printed > n) items_printed = n;

	//Some points were loaded
	if (n > 0)
	{

		//Table  1
		*output << "Results containg values of the criteria:" << std::endl << std::endl;

		//Set properties
		*output << std::showpoint << std::fixed << std::right;

		//Create header 1 : results of analalysis
		*output << std::setw(4) << "#"
			<< std::setw(20) << "Category"
			<< std::setw(8) << "Proj"
			<< std::setw(10) << "R"
			<< std::setw(6) << "latP"
			<< std::setw(7) << "lonP"
			<< std::setw(6) << "lat0"
			<< std::setw(7) << "lon0"
			<< std::setw(9) << "C"
			<< std::setw(10) << "Rotation"
			<< std::setw(13) << "Scale"
			<< std::setw(9) << "Resid"
			<< std::setw(10) << "Iterat"
			<< std::setw(10) << "Evaluat" << std::endl;

		//Values of the criterion for each projection
		for (unsigned int i = 0; i < (analysis_parameters.analyzed_projections.size() == 0 ? items_printed : n); i++)
		{
			sl[i].printSampleRatios2(i + 1, analysis_parameters.analysis_type, n_test, output);
		}


		*output << std::endl;
		*output << std::endl;
	}

}


template <typename T>
void CartAnalysis::printResults3(const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
	const TAnalysisParameters <T> & analysis_parameters, unsigned int index, std::ostream * output)
{

	//Print selected item
	const unsigned int n = sl.size(), n_test = nl_test.size();

	//Some points were loaded
	if ( ( n > 0 ) && ( index < n) )
	{
		//Values of the criterion for each projection
		sl[index].printSampleRatios2(index + 1, analysis_parameters.analysis_type, n_test, output);
	}

}



template <typename T>
T CartAnalysis::getMatchRatioTissotIndicatrix(const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference, const Matrix <T> &X, Matrix <T> &W, const Projection <T> *proj, TIndexList & matched_points_ind,
	const TMatchedPointsCollectType collect_matched_points, const T map_scale, const T lambda)
{
	//Compute match ratio: how many points is inside the ucertainty region formed by Tissot a Indicatrix?
	Container <Node3DCartesianProjected <T> *> nl_projected_temp;
	const unsigned int n = nl_test.size();

	for (unsigned int i = 0; i < n; i++)
	{
		//Get type of the direction
		TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

		//Reduce lon
		const T lon_red = CartTransformation::redLon0(nl_reference[i]->getLon(), X(0, 4));

		//Convert geographic point to oblique aspect
		T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

		//Create Tissot structure, initialize with default parameters independent on (lat, lon): a = 1, XMAX = 1, Ae = 0
		TTissotIndicatrix <T> tiss;

		try
		{
			//Convert geographic point to oblique aspect
			lat_trans = CartTransformation::latToLatTrans(nl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2));
			lon_trans = CartTransformation::lonToLonTrans(nl_reference[i]->getLat(), lon_red, lat_trans, X(0, 1), X(0, 2), trans_lon_dir);

			//Compute x, y coordinates
			x = ArithmeticParser::parseEquation(proj->getXEquatPostfix(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), proj->getC(), X(0, 3), proj->getLat1(), proj->getLat2(), false);
			y = ArithmeticParser::parseEquation(proj->getYEquatPostfix(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), proj->getC(), X(0, 3), proj->getLat1(), proj->getLat2(), false);

			//Create oblique point
			Point3DGeographic <T> point_oblique_temp(lat_trans, lon_trans);
			tiss = CartDistortion::Tiss(NUM_DERIV_STEP, &point_oblique_temp, proj, true);
		}

		catch (Exception &error)
		{
			//Disable point from analysis: set weight to zero
			W(i, i) = 0; W(i + n, i + n) = 0;
		}

		//Create new cartographic point
		Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

		//Set Tissot Indicatrix parameters for the point
		n_projected->setTiss(tiss);

		//Add point to the list
		nl_projected_temp.push_back(n_projected);
	}

	//Computer centers of mass for both systems P, P'
	unsigned int n_points = 0;
	T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

	for (unsigned int i = 0; i < n; i++)
	{
		//Use only non singular points
		if (W(i, i) != 0.0)
		{
			x_mass_test += nl_test[i]->getX();
			y_mass_test += nl_test[i]->getY();

			x_mass_reference += nl_projected_temp[i]->getX();
			y_mass_reference += nl_projected_temp[i]->getY();

			n_points++;
		}
	}

	x_mass_test = x_mass_test / n_points;
	y_mass_test = y_mass_test / n_points;
	x_mass_reference = x_mass_reference / n_points;
	y_mass_reference = y_mass_reference / n_points;

	for (unsigned int i = 0; i < n; i++)
	{
		nl_projected_temp[i]->setX(nl_projected_temp[i]->getX());
		nl_projected_temp[i]->setY(nl_projected_temp[i]->getY());
	}

	//Compute reduced coordinates and add to the containers
	Container <Node3DCartesianProjected <T> * > nl_test_red, nl_projected_red;

	for (unsigned int i = 0; i < n; i++)
	{
		//Create reduced coordinates multiplied by the map scale
		const T x_test_red = (nl_test[i]->getX() - x_mass_test) * map_scale;
		const T y_test_red = (nl_test[i]->getY() - y_mass_test) * map_scale;
		const T x_reference_red = (nl_projected_temp[i]->getX() - x_mass_reference) * map_scale;
		const T y_reference_red = (nl_projected_temp[i]->getY() - y_mass_reference) * map_scale;

		//Creae temporary points
		Node3DCartesianProjected <T> *n1 = new Node3DCartesianProjected <T>(x_test_red, y_test_red);
		Node3DCartesianProjected <T> *n2 = new Node3DCartesianProjected <T>(x_reference_red, y_reference_red);

		//Add points to containers
		nl_test_red.push_back(n1);
		nl_projected_red.push_back(n2);
	}

	//Compute graphical accuracy of the map: uncertainty region = lambda * g.A (lambda = 20-30)
	const T graph_accuracy = 1000 * 1.0e-7 * map_scale;

	//Call Tissot match procedure
	return getMatchRatioTissotIndicatrix(nl_projected_red, nl_test_red, matched_points_ind, collect_matched_points, graph_accuracy, lambda);
}


template <typename Point1, typename Point2>
typename Point1::Type CartAnalysis::getMatchRatioTissotIndicatrix(const Container <Point1 *> &global_points, const Container <Point2 *> &matched_points, TIndexList & matched_points_ind,
	const TMatchedPointsCollectType collect_matched_points, const typename Point1::Type radius, const typename Point1::Type lambda)
{
	//Compute match ratio: how many points is inside the ucertainty region formed by Tissot a Indicatrix?
	unsigned int total_points_matched = 0;
	const unsigned int n = global_points.size();

	//Compute radius of the circle: modified equation given by Wamelen et al (2001)
	const typename Point1::Type max_radius = lambda * radius;

	//Find points satisfying the criterion d((x,y)trans, (x,y)) < min_error
	for (unsigned int i = 0; i < n; i++)
	{
		//Get Tissote Indicatrix parameters
		const TTissotIndicatrix <typename Point1::Type> tiss = global_points[i]->getTiss();

		//Matching condition
		if (PointEllipsePosition::getPointEllipsePosition(global_points[i]->getX(), global_points[i]->getY(), matched_points[i]->getX(), matched_points[i]->getY(), max_radius * tiss.a_tiss, max_radius * tiss.b_tiss, tiss.Ae_proj + tiss.b_mer - 90))
		{
			total_points_matched++;

			if (collect_matched_points == CollectOn) matched_points_ind.push_back(i);
		}
	}

	//Compute match ratio
	return (100.0 / n) * total_points_matched;
}


template <typename Point1, typename Point2>
typename Point1::Type CartAnalysis::getMatchRatioCircle(const Container <Point1 *> &global_points_source, const Container <Point2 *> &transformed_points_dest, TIndexList & matched_points,
	const TMatchedPointsCollectType collect_matched_points, const typename Point1::Type radius, const typename Point1::Type lambda)
{
	//Compute match ratio: how many points is inside the ucertainty region formed by a circle?
	unsigned int total_points_matched = 0;
	const unsigned int n = global_points_source.size();

	//Compute radius of the circle: modified equation given by Wamelen et al (2001)
	const typename Point1::Type max_radius = lambda * radius;

	//Find points satisfying the criterion d((x,y)trans, (x,y)) < min_error
	for (unsigned int i = 0; i < n; i++)
	{
		//Matching condition
		if (EuclDistance::getEuclDistance2D(global_points_source[i]->getX(), global_points_source[i]->getY(), transformed_points_dest[i]->getX(), transformed_points_dest[i]->getY()) < max_radius)
		{
			total_points_matched++;

			if (collect_matched_points == CollectOn) matched_points.push_back(i);
		}
	}

	//Compute match ratio
	return (100.0 / n) * total_points_matched;
}


template <typename T>
T CartAnalysis::solutionDiversity(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, const unsigned int n_items)
{
	//Compute solution diversity: variance of coordinate differences between samples
	unsigned int n_samples = sl.size(), n_points = pl_reference.size();
	T scale = 1.0;

	//The real amount of processed samples
	unsigned int n = std::min(n_samples, n_items);

	//Store coordinates
	Matrix <T> XY(2 * n * n_points, 1), W(2 * n * n_points, 2 * n * n_points, 0, 1),
		XYM(2 * n_points, 1), RXY(2 * n * n_points, 1);

	//Process all samples
	for (unsigned int i = 0; i < n; i++)
	{
		std::cout << "proj " << i << '\n';

		//Get projection and its properties
		Projection <T> *proj = sl[i].getProj();
		proj->setR(sl[i].getR());
		Point3DGeographic <T> cart_pole(sl[i].getLatP(), sl[i].getLonP());
		proj->setCartPole(cart_pole);
		proj->setLat0(sl[i].getLat0());
		proj->setLon0(sl[i].getLon0());
		proj->setDx(sl[i].getDx());
		proj->setDy(sl[i].getDy());
		scale = sl[i].getScaleHomT();

		//Process all points
		for (unsigned int j = 0; j < n_points; j++)
		{
			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Reduce lon
			const T lon_red = CartTransformation::redLon0(pl_reference[j]->getLon(), proj->getLon0());

			try
			{
				//Convert to lat, lon
				const T lat_trans = CartTransformation::latToLatTrans(pl_reference[j]->getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon());
				const T lon_trans = CartTransformation::lonToLonTrans(pl_reference[j]->getLat(), lon_red, proj->getCartPole().getLat(), proj->getCartPole().getLon(), trans_lon_dir);

				//Create temporary point
				Point3DGeographic <T> point_geo_temp(lat_trans, lon_trans);

				//Compute coordiantes
				RXY(i * n_points + j, 0) = CartTransformation::latLonToX(&point_geo_temp, proj, false) - nl_test[j]->getX();
				RXY(i * n_points + j + n * n_points, 0) = CartTransformation::latLonToY(&point_geo_temp, proj, false) - nl_test[j]->getY();
			}

			//Set zero weights
			catch (Exception &error)
			{
				W(i * n_points + j, i * n_points + j) = 0;
				W(i * n_points + j + n * n_points, i * n_points + j + n * n_points) = 0;
			}
		}
	}

	//Sum of squares of residuals = diversity of the solution
	Matrix <T> RES = trans(RXY) * W * RXY;

	return sqrt(RES(0, 0) / (n  * n_points)) * scale / 1000;
}


template <typename T>
void CartAnalysis::batchTestShifts(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
	unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> &best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	const unsigned int n_items = nl_test.size();
	const unsigned int n_tests = 500;

	//Create file names
	char output_file_text_m6[256], output_file_text_m8[256];
	strcpy(output_file_text_m6, analysis_parameters.test_file); strcpy(output_file_text_m8, analysis_parameters.test_file);
	strcat(output_file_text_m6, "m6.res"); strcat(output_file_text_m8, "m8.res");

	static std::ofstream output_file_m6, output_file_m8;

	//Process all shifts
	//1.0e-9
	//1.0e-7
	//1.0e-3 ---EQDC, S1
	//1.0e4
	//3.0e4
	const T res2 = 1.0e0;

	//Initialize random number generator
	srand((unsigned)time(0));

	for (unsigned int i = 0; i < 11; i++)
	{
		//Create variables for measurements
		T min_cost_a1 = 0, min_cost_a2 = 0, min_cost_c1 = 0, min_cost_c2 = 0;
		unsigned int eff_1 = 0, eff_2 = 0, iter_1 = 0, iter_2 = 0, j1 = 0, j2 = 0;

		//Compute shift
		const T shift = std::pow(10, i);

		//Open files
		output_file_m6.open(output_file_text_m6, std::ofstream::app); output_file_m8.open(output_file_text_m8, std::ofstream::app);

		for (unsigned int j = 0; j < n_tests; j++)
		{

			//Create matrices
			Matrix <T> X1(6, 1), Y1(2 * n_items, 1), W1(2 * n_items, 2 * n_items, 0.0, 1), V1(2 * n_items, 1);
			Matrix <T> X2(8, 1), Y2(2 * n_items, 1), W2(2 * n_items, 2 * n_items, 0.0, 1), V2(2 * n_items, 1);

			//Copy shifted points
			Container <Node3DCartesian <T> *> nl_test_temp;
			const unsigned int n_items = nl_test.size();

			for (unsigned int k = 0; k < n_items; k++)
			{
				Node3DCartesian <T> *node_temp = new Node3DCartesian <T>(nl_test[k]->getX() + shift, nl_test[k]->getY() - 1.0*shift);

				//Add shifted point to the list
				nl_test_temp.push_back(node_temp);
			}

			//Use the cartometric analysis to set the initial value of the Earth radius: use the similarity transformation
			unsigned int total_samples_test = 0;
			TAnalysisParameters <T> analysis_parameters_test(false);
			analysis_parameters_test.analysis_type.a_helt = true;

			Sample <T> sample_test;
			//CartAnalysis::computeAnalysisForOneSample(nl_test_temp, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters_test, sample_test, false, total_samples_test, output);

			//Get initialize Earth radius using the similarity transformation
			const T R_def = proj->getR() / sample_test.getScaleHelT();
			T dX0 = proj->getDx();
			T dy_init = proj->getDy();

			const T latp = -90 + 180 * rand() / (RAND_MAX + 1.0);
			const T lonp = -180 + 360 * rand() / (RAND_MAX + 1.0);
			const T lat0 = -90 + 180 * rand() / (RAND_MAX + 1.0);
			const T R = 1.0e8 * rand() / (RAND_MAX + 1.0);
			const T dx = -1.0e8 + 2.0e8 * rand() / (RAND_MAX + 1.0);
			const T dy = -1.0e8 + 2.0e8 * rand() / (RAND_MAX + 1.0);

			std::cout << "\n>> >> Shift = " << shift << "   test = " << j << '\n' << '\n';
			//std::cout << "DX = " << dX0 << "  DY = " << dy_init << '\n';

			//Initialize X1
			//Use the  values as the initial vector for the least squares solution
			//X1(0, 0) = 10;
			//X1(1, 0) = latp;
			//X1(2, 0) = lonp;
			//X1(3, 0) = 10;
			//X1(4, 0) = 0;
			//X1(5, 0) = 1;

			X1(0, 0) = R;
			X1(1, 0) = latp;
			X1(2, 0) = lonp;
			X1(3, 0) = lat0;
			X1(4, 0) = 0;
			X1(5, 0) = 1;

			Matrix <T> XI = X1;
			//T min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_temp, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j1), FAnalyzeProjV2 <T>(nl_test_temp, pl_reference, meridians, parallels, faces_test, proj,
			//	analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W1, X1, Y1, V1, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);

			T min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ2 <T>(nl_test_temp, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j1), FAnalyzeProjV2 <T>(nl_test_temp, pl_reference, meridians, parallels, faces_test, proj,
				analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W1, X1, Y1, V1, iterations, alpha_mult, eps, max_iter, max_diff, output);

			iter_1 += iterations;
			min_cost_a1 += min_cost;

			std::cout << "M6 = (" << min_cost << ")";

			//X1.print();

			if (min_cost < res2)
			{
				eff_1++;
				min_cost_c1 += min_cost;
			}

			else
			{
				std::cout << "M6: \n";
				XI.print();
				X1.print();
				//X1.print(output);
			}

			//Initialize X2
			//X2(0, 0) = 10;
			//X2(1, 0) = latp;
			//X2(2, 0) = lonp;
			//X2(3, 0) = 10;
			//X2(4, 0) = 0;
			//X2(5, 0) = dX0;
			//X2(6, 0) = dy_init;
			//X2(5, 0) = 1.0;
			//X2(6, 0) = 1.0;
			//X2(7, 0) = 1;

			X2(0, 0) = R;
			X2(1, 0) = latp;
			X2(2, 0) = lonp;
			X2(3, 0) = lat0;
			X2(4, 0) = 0;
			//X2(5, 0) = dX0;
			//X2(6, 0) = dy_init;
			X2(5, 0) = dx;
			X2(6, 0) = dy;
			X2(7, 0) = 1;

			//min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_temp, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j2), FAnalyzeProjV <T>(nl_test_temp, pl_reference, meridians, parallels, faces_test, proj,
			//	analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);

			min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ <T>(nl_test_temp, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j2), FAnalyzeProjV <T>(nl_test_temp, pl_reference, meridians, parallels, faces_test, proj,
				analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, iterations, alpha_mult, eps, max_iter, max_diff, output);

			iter_2 += iterations;
			min_cost_a2 += min_cost;

			//X2.print();

			if (min_cost < res2)
			{
				eff_2++;
				min_cost_c2 += min_cost;
			}

			else
			{
				std::cout << "M8: \n";
				X2.print();
				//X2.print(output);
			}

			std::cout << "   M8 = (" << min_cost << ") \n";
		}

		//Set output
		output = &output_file_m6;

		//Write result to file
		*output << shift << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\t' << iter_1 << '\t' << j1 << '\t' << eff_1 * 100.0 / n_tests << '\n';

		//Set output
		output = &output_file_m8;

		//Write results to file
		*output << shift << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_2 << '\t' << j2 << '\t' << eff_2 * 100.0 / n_tests << '\n';

		std::cout << "M6 = (" << min_cost_a1 << ") " << "   M8 = (" << min_cost_a2 << ") \n";

		//Close files
		output_file_m6.close(); output_file_m8.close();
	}
}




template <typename T>
void CartAnalysis::batchTestNLSP(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned short &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Perform batch tests
	unsigned short n_tests = 300;
	const T nu1 = 0.25, nu2 = 0.75, nu3 = 0.01, gamma1 = 0.25, gamma2 = 2.0, lambda_min = 1.0e-6, lambda_max = 1.0e6;

	const unsigned short n_items = nl_test.size();

	//Initialize random number generatorres
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_gnd[256], output_file_text_bfgs[256], output_file_text_bfgsh[256],
		output_file_text_lm[256], output_file_text_tr[256];

	strcpy(output_file_text_gnd, "gnd.log"); strcpy(output_file_text_bfgs, "bfgs.log"); strcpy(output_file_text_bfgsh, "bfgsh.log");
	strcpy(output_file_text_lm, "lm.log"); strcpy(output_file_text_tr, "tr.log");

	//Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 6380;

	//New output file variables
	static std::ofstream output_file_gnd, output_file_bfgs, output_file_bfgsh, output_file_lm, output_file_tr;

	//Perform for all noises
	for (int noise = 3; noise < 33; noise += 3)
	{
		std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

		T delta = 100;
		//delta = 100;
		//for (T delta = 1; delta < 100000; delta *= 10)
		{

			//Open output log files
			output_file_gnd.open(output_file_text_gnd, std::ofstream::app); output_file_bfgs.open(output_file_text_bfgs, std::ofstream::app); output_file_bfgsh.open(output_file_text_bfgsh, std::ofstream::app); output_file_lm.open(output_file_text_lm, std::ofstream::app); output_file_tr.open(output_file_text_tr, std::ofstream::app);

			//Create variables for measurements
			T min_cost_a1 = 0, min_cost_a2 = 0, min_cost_a3 = 0, min_cost_a4 = 0, min_cost_a5 = 0, min_cost_c1 = 0, min_cost_c2 = 0, min_cost_c3 = 0, min_cost_c4 = 0, min_cost_c5 = 0;
			unsigned int eff_1 = 0, eff_2 = 0, eff_3 = 0, eff_4 = 0, eff_5 = 0, iter_1 = 0, iter_2 = 0, iter_3 = 0, iter_4 = 0, iter_5 = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0;

			std::cout << "\n>> >> delta = " << delta << '\n' << '\n';
			//std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

			//Measure time difference
			time_t start, end;
			time(&start);

			//Perform all tests
			for (unsigned int j = 0; j < n_tests; j++)
			{
				//Lowest and highest values, range
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_a1 = (-noise / 100.0) * 360, high_a1 = (noise / 100.0) * 360;
				T low_a2 = (-noise / 100.0) * 180, high_a2 = (noise / 100.0) * 180;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_a1 = high_a1 - low_a1;
				T range_a2 = high_a2 - low_a2;
				T range_dx = high_dx - low_dx;

				//Random values
				T rand_R = low_R + range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp = std::min(std::max(-90.0, latp_init + low_a2 + range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp = std::min(std::max(-180.0, lonp_init + low_a1 + range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0 = std::min(std::max(lat0_min, lat0_init + low_a2 + range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx = low_dx + range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy = low_dx + range_dx*rand() / (RAND_MAX + 1.0);

				unsigned int n_par = 6;
				if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					n_par = 7;
				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					n_par = 5;
				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					n_par = 8;

				Matrix <T> X(n_par, 1), XMIN(n_par, 1), XMAX(n_par, 1), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1);

				//Store randomly generated values
				X(0, 0) = rand_R;
				X(1, 0) = rand_latp;
				X(2, 0) = rand_lonp;
				X(3, 0) = rand_lat0;

				//X(0, 0) = 5989.7376929;
				//X(1, 0) = -14.4607544;
				//X(2, 0) = 32.7027832;
				//X(3, 0) = 0;

				//Set intervals
				XMIN(0, 0) = low_R; XMAX(0, 0) = low_R + range_R;
				XMIN(1, 0) = latp_min; XMAX(1, 0) = latp_max;
				XMIN(2, 0) = lonp_min; XMAX(2, 0) = lonp_max;
				XMIN(3, 0) = lat0_min; XMAX(3, 0) = lat0_max;
				XMIN(4, 0) = 0.0; XMAX(4, 0) = 0.0;

				//Set intervals
				if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
				{
					XMIN(5, 0) = 0.0; XMAX(5, 0) = 10000000;
				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
				{
					X(5, 0) = 1;
					X(6, 0) = 80;

					XMIN(5, 0) = 0.0; XMAX(5, 0) = 10000000;
					XMIN(6, 0) = -MAX_LON; XMAX(6, 0) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
				{
					X(0, 0) = rand_latp;
					X(1, 0) = rand_lonp;
					X(2, 0) = rand_lat0;
					X(3, 0) = 0;
					X(4, 0) = 1;

					XMIN(0, 0) = latp_min; XMAX(0, 0) = latp_max;
					XMIN(1, 0) = lonp_min; XMAX(1, 0) = lonp_max;
					XMIN(2, 0) = lat0_min; XMAX(2, 0) = lat0_max;
					XMIN(3, 0) = 0.0; XMAX(3, 0) = 0.0;
					XMIN(4, 0) = 0.0; XMAX(4, 0) = 10000000;
				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
				{
					//X(5, 0) = rand_dx;
					//X(6, 0) = rand_dy;
					X(7, 0) = 1;

					XMIN(5, 0) = -1.0e09; XMAX(5, 0) = 1.0e9;
					XMIN(6, 0) = -1.0e09; XMAX(6, 0) = 1.0e9;
					XMIN(7, 0) = 0.0; XMAX(7, 0) = 10000000;
				}
				/*
				//*************Additionaly added ***************
				XX(0, 0) = 6201.4682959;
				XX(1, 0) = -8.6312988;
				XX(2, 0) = 99.8191406;
				XX(3, 0) = 55.8554932;
				*/
				/*
				XX(0, 0) = 6358.1232495;
				XX(1, 0) = -3.6805298;
				XX(2, 0) = 91.8839355;
				XX(3, 0) = 51.2771606;
				*/
				/*
				XX(0, 0) = 6056.8188062;
				XX(1, 0) = 7.0163086;
				XX(2, 0) = 69.5324707;
				XX(3, 0) = 41.8558350;

				XX(0, 0) = 6025.0066626;
				XX(0, 0) = -2.6611084   ;
				XX(0, 0) = 108.8235352  ;
				XX(0, 0) =  53.3862061;


				XX(0, 0) = 6684.3962549;
				XX(1, 0) = 12.8553223;
				XX(2, 0) = 117.4851562;
				XX(3, 0) = 49.8879395;
				*/
				/*
				XX(0, 0) =4937.3443359;
				XX(1, 0) =-50.2756348 ;
				XX(2, 0) =118.4370117 ;
				XX(3, 0) =87.5600586;
				*/
				//XX(0, 0) = 6078.7976807;
				//XX(1, 0) = -18.9316406;
				//XX(2, 0) = 160.9211426;
				//XX(3, 0) =  70.2994385;

				/*
				XX(0, 0) = 4734.9565063;
				XX(1, 0) = 4.6604004;
				XX(2, 0) = 105.7543945;
				XX(3, 0) = 57.6563721;

				X2(0, 0) = 5641.3223145;
				X2(1, 0) = 90.0000000;
				X2(2, 0) = 54.7185059;
				X2(3, 0) = 51.8251953;
				X2(5, 0) = 0;

				*/
				//**********************************************

				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << j << "/" << n_tests << '\n';
				std::cout << "R= " << rand_R << "   latp= " << rand_latp << "   lonp= " << rand_lonp << "   lato0= " << rand_lat0 << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 2.0e8;

				//Create copy for analysis
				Matrix<double> X2 = X, Y2 = Y, V2 = V, W2 = W;
				nl_projected.clear();

				//Gauss-Newton
				std::cout << "GND ";
				output = &output_file_gnd;
				T min_cost = 0, q1 = 1, q2 = 1, R_def = rand_R;

				bool test_gnd = false, test_bfgs = false, test_bfgsh = true, test_lm = false, test_tr = false;
				if (test_gnd)
				{
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j1), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, j1), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, aspect, R, q1, q2, analysis_parameters.print_exceptions, j1), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::GND(FAnalyzeProjJ <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j1), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);

					//T min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, (*i_projections), ObliqueAspect, q1, q2, analysis_parameters.print_exceptions), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, *i_projections,
					//	R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X0, Y, V, iterations, alpha_mult, nu, eps, 0.25 * max_iter, max_diff, output) : 0);


					min_cost_a1 += min_cost;
					iter_1 += iterations;
					std::cout << "(" << iterations << ") ";

					//Test results
					//bool correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
					//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
					//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));

					if (/*correct &&*/ min_cost < res2)
					{
						eff_1++;
						min_cost_c1 += min_cost;
					}

					else
					{
						X2.print();
					}
				}

				if (test_bfgs)
				{

					//BFGS
					std::cout << " BFGS ";
					q1 = 1, q2 = 1, R_def = rand_R;
					output = &output_file_bfgsh;
					X2 = X; Y2 = Y; V2 = V;  W2 = W;
					nl_projected.clear();
					/*
					R_def = 6445.2562903;
					X2(0, 0) = 90.0000000;
					X2(1, 0) = 6.5772949;
					X2(2, 0) = 44.0307251;
					*/

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGS(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j2), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGS(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, j2), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGS(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, aspect, R, q1, q2, analysis_parameters.print_exceptions, j2), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGS(FAnalyzeProjJ <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j2), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, eps, max_iter, max_diff, output);

					min_cost_a2 += min_cost;
					iter_2 += iterations;
					std::cout << "(" << iterations << ") ";

					//correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
					//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
					//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));

					if (/*correct &&*/ min_cost < res2)
					{
						eff_2++;
						min_cost_c2 += min_cost;
					}

					else
					{
						X2.print();
					}
				}

				if (test_bfgsh)
				{
					//BFGSH
					//T min_cost;
					std::cout << " BFGSH ";
					q1 = 1, q2 = 1, R_def = rand_R;
					output = &output_file_bfgsh;
					X2 = X; X2 = X; Y2 = Y; V2 = V;  W2 = W;
					nl_projected.clear();


					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j3), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, j3), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, aspect, R, q1, q2, analysis_parameters.print_exceptions, j3), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j3), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					min_cost_a3 += min_cost;
					iter_3 += iterations;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					//correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
					//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
					//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));


					if (/*correct &&*/ min_cost < res2)
					{
						eff_3++;
						min_cost_c3 += min_cost;
					}

					else
					{
						X2.print();
					}
				}

				if (test_lm)
				{
					//Levenberg-Marquardt
					std::cout << " LM ";
					q1 = 1, q2 = 1, R_def = rand_R;
					output = &output_file_lm;
					X2 = X; Y2 = Y; V2 = V;  W2 = W;
					nl_projected.clear();

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::LM(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j4), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::LM(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, j4), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::LM(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, aspect, R, q1, q2, analysis_parameters.print_exceptions, j4), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::LM(FAnalyzeProjJ <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions, j4), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);


					min_cost_a4 += min_cost;
					iter_4 += iterations;
					std::cout << "(" << iterations << ") ";

					//correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
					//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
					//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));

					if (/*correct &&*/ min_cost < res2)
					{
						eff_4++;
						min_cost_c4 += min_cost;
					}

					else
					{
						X2.print();
					}
				}

				if (test_tr)
				{
					//Trust region
					//T min_cost;
					std::cout << " TR ";
					q1 = 1, q2 = 1, R_def = rand_R;
					output = &output_file_tr;
					X2 = X; Y2 = Y; V2 = V;  W2 = W;
					nl_projected.clear();

					/*
					R_def = XXT(0, 0);
					XXT2(0, 0) = XXT(1,0);
					XXT2(1, 0) = XXT(2,0);
					XXT2(2, 0) = XXT(3,0);
					XXT2(3, 0) = XXT(4,0);
					XXT2(4, 0) = XXT(5,0);

					//R_def = 7843.7643677;
					//XXT2(0, 0) = 37.6457520;
					//XXT2(1, 0) = 133.6301270;
					//XXT2(2, 0) = 2.5852051;
					//XXT2(4, 0) = 1.0;
					//alpha_mult = 0.1;
					//nu = 0.01;

					//XXT(0, 0) = 5380;
					//XXT(1, 0) = 60;
					//XXT(2, 0) = 20;
					//XXT(4, 0) = 50.0;

					//XXT2(0, 0) = 10;
					//XXT2(1, 0) = 168;
					//XXT2(2, 0) = 50;
					//XXT2(4, 0) = 1.0;


					T min_cost = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod ? NonLinearLeastSquares::LMTR3(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, ObliqueAspect, q1, q2, analysis_parameters.print_exceptions, j5), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), WW, XXT2, YY, VV, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, output) : 0);

					// min_cost = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod ? NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, ObliqueAspect, q1, q2, analysis_parameters.print_exceptions), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
					//	R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), WW, XXT2, YY, VV, iterations, alpha_mult, nu, eps, max_iter, max_diff, output) : 0);


					//XXT2.print();
					*/

					//X2(0, 0) = 6387.4595581;
					//X2(1, 0) = 69.4913330;
					//X2(2, 0) = -54.1381836;
					//X2(3, 0) = -3.5195313;

					//X2(0, 0) = 5604.89;
					//X2(1, 0) = 85.9412;
					//X2(2, 0) = 56.7754;
					//X2(3, 0) = 81.5278;

					//X2(0, 0) = 6913.6912720;
					//X2(1, 0) = 2.5048828;
					//X2(2, 0) = -16.0620117;
					//X2(3, 0) = 80.7829590;

					//X2(0, 0) = 7839.9104736;
					//X2(1, 0) = -47.6553955;
					//X2(2, 0) = 163.4721680;
					//X2(3, 0) = 29.5053711;

					//R_def = 5158.8848511;
					//X2(0, 0) = 38.4433594;
					//X2(1, 0) = 162.2988281;
					//X2(2, 0) = 76.9761963;
					//X2(3, 0) = 0;
					//X2(4, 0) = 0;
					//X2(0, 0) = 8282.8747192;
					//X2(1, 0) = 90.0000000;
					//X2(2, 0) = 21.5222168;
					//X2(3, 0) = 61.7590332;

					//R_def = 6788.0309692;
					//X2(0, 0) = -39.6859131;
					//X2(1, 0) = 153.4724121;
					//X2(2, 0) = 12.7523193;


					////////Merc,T2
					//TR: M6, M8> tR0 = 100, add_step
					//TR3: M7S: tR0 = 20*ng, reflection or 3 000 000;
					//TR: M7, addd_step, tR0 = 100

					////////EQDC, T2
					// TR3: M7S : tR0 = 20 * ng;
					// TR2: M6
					// GND: add_step

					//////Eqdc, T1
					//TR: M7, add_Step, 100
					//TR3: M7S, add_step, 10
					//BFGS: M7S. add step

					///////Merc, T1
					//TR3: M7s, tR0 = 3 000 000, reflection
					//TR3: M7, tR0 = 100, reflection,
					/*
					R_def = 6264.2517334;
					R_def = 6380;
					X2(0, 0) = 90.0000000;
					X2(1, 0) = -95.4558105;
					X2(2, 0) = -5.9519043;
					X2(2, 0) = 45.9519043;
					*/
					/*
					R_def = 7896.2006836;
					X2(0, 0) = 90.0000000;
					X2(1, 0) = -74.6520996;
					X2(2, 0) = 76.4290771;
					*/

					//X2(0, 0) = 7062.5917236;
					//X2(1, 0) = -20.9652100;
					//X2(2, 0) = 168.0864258;
					//X2(3, 0) = 42.3297119;
					/*
									R_def = 7062.5917236;
									X2(0, 0) = -20.9652100;
									X2(1, 0) = 168.0864258;
									X2(2, 0) = 42.3297119;
									*/
					//R_def = 6380.5917236;
					//X2(0, 0) = 0.9652100;
					//X2(1, 0) = 88.0864258;
					//X2(2, 0) = 52.3297119;

					//T min_cost = 0;i
					R = 1;
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::TR(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, ObliqueAspect, analysis_parameters.print_exceptions, j5), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::TR(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, ObliqueAspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, j5), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::TR3(FAnalyzeProjJ4 <T>(nl_test, pl_reference, proj, aspect, R, q1, q2, analysis_parameters.print_exceptions, j5), FAnalyzeProjV4 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::TR(FAnalyzeProjJ <T>(nl_test, pl_reference, proj, ObliqueAspect, analysis_parameters.print_exceptions, j5), FAnalyzeProjV <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
						analysis_parameters, ObliqueAspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y2, V2, XMIN, XMAX, iterations, nu1, nu2, nu3, gamma1, gamma2, lambda_min, lambda_max, eps, max_iter, max_diff, delta, output);


					min_cost_a5 += min_cost;
					iter_5 += iterations;


					//bool correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
					//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
					//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));

					if (min_cost < res2)
					{
						eff_5++;
						min_cost_c5 += min_cost;
					}

					else
					{
						X2.print();
					}

					std::cout << 100.0 * eff_5 / j;
					std::cout << "(" << iterations << ")  cost=" << min_cost << "\n";

				}
			}

			//End time
			time(&end);
			float time_diff = difftime(end, start);

			//Write results to files
			output = &output_file_gnd;
			*output << noise << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\t' << iter_1 << '\t' << j1 << '\t' << eff_1 * 100.0 / n_tests << '\n';
			//*output << "Noise [%] = " << noise << '\n';
			//*output << "SUM_A = " << min_cost_a1 << '\n' << "SUMC = " << min_cost_c1 << '\n' << "IT = " << iter_1 << '\n' << "EFF = " << eff_1 * 100.0 / n_tests << '\n';

			output = &output_file_bfgs;
			*output << noise << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_2 << '\t' << j2 << '\t' << eff_2 * 100.0 / n_tests << '\n';
			//*output << "Noise [%] = " << noise << '\n';
			//*output << "SUM_A = " << min_cost_a2 << '\n' << "SUMC = " << min_cost_c2 << '\n' << "IT = " << iter_2 << '\n' << "EFF = " << eff_2 * 100.0 / n_tests << '\n';
			//*output << "TIME = " << time_diff << '\n';

			output = &output_file_bfgsh;
			*output << noise << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\t' << iter_3 << '\t' << j3 << '\t' << eff_3 * 100.0 / n_tests << '\n';
			//*output << "Noise [%] = " << noise << '\n';
			//*output << "SUM_A = " << min_cost_a3 << '\n' << "SUMC = " << min_cost_c3 << '\n' << "IT = " << iter_3 << '\n' << "EFF = " << eff_3 * 100.0 / n_tests << '\n';

			output = &output_file_lm;
			*output << noise << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\t' << iter_4 << '\t' << j4 << '\t' << eff_4 * 100.0 / n_tests << '\n';
			//*output << "Noise [%] = " << noise << '\n';
			//*output << "SUM_A = " << min_cost_a4 << '\n' << "SUMC = " << min_cost_c4 << '\n' << "IT = " << iter_4 << '\n' << "EFF = " << eff_4 * 100.0 / n_tests << '\n';

			output = &output_file_tr;
			*output << noise << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\t' << iter_5 << '\t' << j5 << '\t' << eff_5 * 100.0 / n_tests << '\n';
			//*output << "Noise [%] = " << noise << '\n';
			//*output << "SUM_A = " << min_cost_a5 << '\n' << "SUMC = " << min_cost_c5 << '\n' << "IT = " << iter_5 << '\n' << "EFF = " << eff_5 * 100.0 / n_tests << '\n';

			//Close files
			output_file_gnd.close(); output_file_bfgs.close(); output_file_bfgsh.close(); output_file_lm.close(); output_file_tr.close();
		}
	}
}


template <typename T>
void CartAnalysis::batchTestSimplex(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned short &iterations, const T eps, unsigned short max_iter, const T max_diff, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Perform batch tests
	unsigned short n_tests = 300;

	const unsigned short n_items = nl_test.size();

	//Initialize random number generator
	srand((unsigned)time(0));

	//Create file names
	char output_file_text[256];

	strcpy(output_file_text, "simplex.log");

	//Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 6380;

	//New output file variables
	static std::ofstream output_file;

	//Perform for all noises
	for (int noise = 100; noise < 150; noise += 50)
	{
		std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

		T delta = 100;
		//delta = 100;
		//for (T delta = 1; delta < 100000; delta *= 10)
		{

			//Open output log files
			output_file.open(output_file_text, std::ofstream::app);

			//Create variables for measurements
			T min_cost_a = 0, min_cost_c = 0;
			unsigned int eff = 0, iter = 0, it_res = 0;

			std::cout << "\n>> >> delta = " << delta << '\n' << '\n';
			//std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

			//Measure time difference
			time_t start, end;
			time(&start);

			//Perform all tests
			for (unsigned int j = 0; j < n_tests; j++)
			{

				//Lowest and highest values, range
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_a1 = (-noise / 100.0) * 360, high_a1 = (noise / 100.0) * 360;
				T low_a2 = (-noise / 100.0) * 180, high_a2 = (noise / 100.0) * 180;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_a1 = high_a1 - low_a1;
				T range_a2 = high_a2 - low_a2;
				T range_dx = high_dx - low_dx;

				//Min random values
				T rand_R_min = R0 - 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_min = std::min(std::max(-90.0, latp_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_min = std::min(std::max(-180.0, lonp_init - 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_min = std::min(std::max(lat0_min, lat0_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				//Max random values
				T rand_R_max = R0 + 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_max = std::min(std::max(-90.0, latp_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_max = std::min(std::max(-180.0, lonp_init + 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_max = std::min(std::max(lat0_min, lat0_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				unsigned int n_par = 6;
				if (analysis_parameters.analysis_method == SimplexRotMethod)
					n_par = 7;
				else if (analysis_parameters.analysis_method == SimplexRot2Method)
					n_par = 5;
				else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
					n_par = 8;

				Matrix <T> XMIN(1, n_par), X(1, n_par, 1), XMAX(1, n_par, 1), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1);

				//Minimum values
				XMIN(0, 0) = rand_R_min;
				XMIN(0, 1) = rand_latp_min;
				XMIN(0, 2) = rand_lonp_min;
				XMIN(0, 3) = rand_lat0_min;

				//Maximum values
				XMAX(0, 0) = rand_R_max;
				XMAX(0, 1) = rand_latp_max;
				XMAX(0, 2) = rand_lonp_max;
				XMAX(0, 3) = rand_lat0_max;
				XMAX(0, 4) = 0.0;

				//Set intervals
				if (analysis_parameters.analysis_method == SimplexMethod)
				{
					XMAX(0, 5) = 10000000;
				}

				else if (analysis_parameters.analysis_method == SimplexRotMethod)
				{
					XMIN(0, 5) = 0;
					XMIN(0, 6) = -MAX_LON;

					XMAX(0, 5) = 10000000;
					XMAX(0, 6) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == SimplexRot2Method)
				{
					XMIN(0, 0) = rand_latp_min;
					XMIN(0, 1) = rand_lonp_min;
					XMIN(0, 2) = rand_lat0_min;
					XMIN(0, 3) = 0;
					XMIN(0, 4) = 0;

					XMAX(0, 0) = rand_latp_max;
					XMAX(0, 1) = rand_lonp_max;
					XMAX(0, 2) = rand_lat0_max;
					XMAX(0, 3) = 0.0;
					XMAX(0, 4) = 10000000;
				}

				else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
				{
					XMIN(0, 5) = -1.0e03;
					XMIN(0, 6) = -1.0e03;
					XMIN(0, 7) = 0;

					XMAX(0, 5) = 1.0e3;
					XMAX(0, 6) = 1.0e3;
					XMAX(0, 7) = 10000000;
				}

				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << j << "/" << n_tests << '\n';

				XMIN.print(); XMAX.print();

				T min_cost = 0;
				T q1 = 1, q2 = 1, R_def = rand_R_min, k = 2.0;
				Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

				if (analysis_parameters.analysis_method == SimplexMethod)
					min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, X, Y, V, iterations, eps, max_iter, output);
				else if (analysis_parameters.analysis_method == SimplexRotMethod)
					min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, X, Y, V, iterations, eps, max_iter, output);
				else if (analysis_parameters.analysis_method == SimplexRot2Method)
					min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, W, X, Y, V, iterations, eps, max_iter, output);
				else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
					min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, X, Y, V, iterations, eps, max_iter, output);

				min_cost_a += min_cost;
				iter += iterations;


				//bool correct = (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod || analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
				//	((fabs(fabs(XXT(0, 0)) - fabs(XR(0, 0))) < XR(0, 0) / 50) && (fabs(fabs(XXT(1, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT(2, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT(3, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/) :
				//	((fabs(R_def - XR(0, 0)) < XR(0, 0) / 50) && (fabs(fabs(XXT2(0, 0)) - fabs(XR(1, 0))) < 180.0 / 50) && (fabs(fabs(XXT2(1, 0)) - fabs(XR(2, 0))) < 360.0 / 50) /*&& (fabs(fabs(XXT2(2, 0)) - fabs(XR(3, 0))) < 180.0 / 50)*/));

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 1.0e2;

				if (min_cost < res2)
				{
					eff++;
					min_cost_c += min_cost;
				}

				else
				{
					X.print();
				}

				std::cout << 100.0 * eff / j;
				std::cout << "(" << iterations << ")  cost=" << min_cost << "\n";

			}


			//End time
			time(&end);
			float time_diff = difftime(end, start);

			//Write results to files
			output = &output_file;
			*output << noise << '\t' << min_cost_a << '\t' << min_cost_c << '\t' << iter << '\t' << it_res << '\t' << eff * 100.0 / n_tests << '\n';

			//Close files
			output_file.close();
		}

	}
}


template <typename T>
void CartAnalysis::batchTestDiffEvolutionSchema(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{

	//Perform batch tests of the 6 differential evolution mutation schema
	unsigned short n_tests = 100;
	const unsigned short n_items = nl_test.size();
	const T CR = 0.8;

	unsigned int n_par = 6, m = nl_test.size();
	if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
		n_par = 7;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
		n_par = 5;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
		n_par = 8;

	//Create matrices
	//Matrix <double> XMI(1, n_par), XMAX(1, n_par), A2(1, n_par), B2(1, n_par), X(1, n_par), Y(2 * m, 1);
	Matrix <T> F(1, 1); F(0, 0) = 0.5;

	//Parameters of the genetic algorithm
	const unsigned int population = n_par * n_par, max_gen = 1000;

	//Initialize random number generator
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_de_rand1[256], output_file_text_de_rand2[256], output_file_text_de_rand_dir1[256], output_file_text_de_rand_dir2[256], output_file_text_de_rand_best1[256], output_file_text_de_rand_best2[256], output_file_text_de_rand_best_dir1[256];

	strcpy(output_file_text_de_rand1, "de_rand1.log"); strcpy(output_file_text_de_rand2, "de_rand2.log"); strcpy(output_file_text_de_rand_dir1, "de_rand_dir1.log"); strcpy(output_file_text_de_rand_dir2, "de_rand_dir2.log");
	strcpy(output_file_text_de_rand_best1, "de_rand_best1.log"); strcpy(output_file_text_de_rand_best2, "de_rand_best2.log"); strcpy(output_file_text_de_rand_best_dir1, "de_rand_best_dir1.log");

	//Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 6380.0;

	//New output file variables
	static std::ofstream output_file_de_rand1, output_file_de_rand2, output_file_de_rand_dir1, output_file_de_rand_dir2, output_file_de_rand_best1, output_file_de_rand_best2, output_file_de_rand_best_dir1;

	//Perform for all noises
	for (int noise = 100; noise < 150; noise += 50)
	{
		std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

		T delta = 100.0;
		//delta = 100;
		//for (T delta = 1; delta < 100000; delta *= 10)
		{

			//Open output log files
			output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app); output_file_de_rand2.open(output_file_text_de_rand2, std::ofstream::app);
			output_file_de_rand_dir1.open(output_file_text_de_rand_dir1, std::ofstream::app); output_file_de_rand_best1.open(output_file_text_de_rand_best1, std::ofstream::app);
			output_file_de_rand_best2.open(output_file_text_de_rand_best2, std::ofstream::app); output_file_de_rand_best_dir1.open(output_file_text_de_rand_best_dir1, std::ofstream::app);
			output_file_de_rand_dir2.open(output_file_text_de_rand_dir2, std::ofstream::app);

			//Create variables for measurements
			T min_cost_a1 = 0, min_cost_a2 = 0, min_cost_a3 = 0, min_cost_a4 = 0, min_cost_a5 = 0, min_cost_a6 = 0, min_cost_a7 = 0, min_cost_c1 = 0, min_cost_c2 = 0, min_cost_c3 = 0, min_cost_c4 = 0, min_cost_c5 = 0, min_cost_c6 = 0, min_cost_c7 = 0,
				res_aver1 = 0, res_aver2 = 0, res_aver3 = 0, res_aver4 = 0, res_aver5 = 0, res_aver6 = 0, res_aver7 = 0, res_max1 = 0, res_max2 = 0, res_max3 = 0, res_max4 = 0, res_max5 = 0, res_max6 = 0, res_max7 = 0,
				res_diff1 = 0, res_diff2 = 0, res_diff3 = 0, res_diff4 = 0, res_diff5 = 0, res_diff6 = 0, res_diff7 = 0;
			unsigned int eff_1 = 0, eff_2 = 0, eff_3 = 0, eff_4 = 0, eff_5 = 0, eff_6 = 0, eff_7 = 0, iter_1 = 0, iter_2 = 0, iter_3 = 0, iter_4 = 0, iter_5 = 0, iter_6 = 0, iter_7 = 0, it_res1 = 0, it_res2 = 0, it_res3 = 0, it_res4 = 0, it_res5 = 0, it_res6 = 0, it_res7 = 0;

			std::cout << "\n>> >> delta = " << delta << '\n' << '\n';
			//std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

			//Measure time difference
			time_t start, end;
			time(&start);

			//Perform all tests
			for (unsigned int j = 0; j < n_tests; j++)
			{

				//Lowest and highest values, range
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_a1 = (-noise / 100.0) * 360, high_a1 = (noise / 100.0) * 360;
				T low_a2 = (-noise / 100.0) * 180, high_a2 = (noise / 100.0) * 180;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_a1 = high_a1 - low_a1;
				T range_a2 = high_a2 - low_a2;
				T range_dx = high_dx - low_dx;

				//Min random values
				T rand_R_min = R0 - 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_min = std::min(std::max(-90.0, latp_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_min = std::min(std::max(-180.0, lonp_init - 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_min = std::min(std::max(lat0_min, lat0_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				//Max random values
				T rand_R_max = R0 + 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_max = std::min(std::max(-90.0, latp_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_max = std::min(std::max(-180.0, lonp_init + 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_max = std::min(std::max(lat0_min, lat0_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				Matrix <T> XMIN(1, n_par), X(1, n_par, 1), XMAX(1, n_par, 1), XAVER(1, n_par), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1);

				//Minimum values
				XMIN(0, 0) = rand_R_min;
				XMIN(0, 1) = rand_latp_min;
				XMIN(0, 2) = rand_lonp_min;
				XMIN(0, 3) = rand_lat0_min;
				XMIN(0, 4) = 0.0;

				//Maximum values
				XMAX(0, 0) = rand_R_max;
				XMAX(0, 1) = rand_latp_max;
				XMAX(0, 2) = rand_lonp_max;
				XMAX(0, 3) = rand_lat0_max;
				XMAX(0, 4) = 0.0;

				//Set intervals
				if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
				{
					XMAX(0, 5) = 10000000;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
				{
					XMIN(0, 5) = 0;
					XMIN(0, 6) = -MAX_LON;

					XMAX(0, 5) = 10000000;
					XMAX(0, 6) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
				{
					XMIN(0, 0) = rand_latp_min;
					XMIN(0, 1) = rand_lonp_min;
					XMIN(0, 2) = rand_lat0_min;
					XMIN(0, 3) = 0;
					XMIN(0, 4) = 0;

					XMAX(0, 0) = rand_latp_max;
					XMAX(0, 1) = rand_lonp_max;
					XMAX(0, 2) = rand_lat0_max;
					XMAX(0, 3) = 0.0;
					XMAX(0, 4) = 10000000;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
				{
					XMIN(0, 5) = -1.0e03;
					XMIN(0, 6) = -1.0e03;
					XMIN(0, 7) = 0;

					XMAX(0, 5) = 1.0e3;
					XMAX(0, 6) = 1.0e3;
					XMAX(0, 7) = 10000000;
				}

				XMIN.print(); XMAX.print();

				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << j << "/" << n_tests << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 2.0e8;

				std::cout << "\n DE_RAND1 \n";
				output = &output_file_de_rand1;

				T min_cost = 0, res_aver = 0, res_max = 0, res_diff = 0;
				T q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);

				bool test_de_rand1 = false, test_de_rand2 = false, test_de_rand_dir1 = false, test_de_rand_dir2 = false, test_de_rand_best1 = true, test_de_rand_best2 = false, test_de_rand_best_dir1 = false;

				if (test_de_rand1)
				{
					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res1, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res1, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res1, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res1, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a1 += min_cost;
					iter_1 += iterations;
					res_aver1 += res_aver;
					res_max1 += res_max;
					res_diff1 += res_diff;

					if (min_cost < res2)
					{
						eff_1++;
						min_cost_c1 += min_cost;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_rand2)
				{
					std::cout << "\n DE_RAND2 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand2;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res2, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res2, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res2, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res2, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a2 += min_cost;
					iter_2 += iterations;
					res_aver2 += res_aver;
					res_max2 += res_max;
					res_diff2 += res_diff;

					if (min_cost < res2)
					{
						eff_2++;
						min_cost_c2 += min_cost;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_rand_dir1)
				{
					std::cout << "\n DE_RAND_DIR1 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand_dir1;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res3, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res3, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res3, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res3, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a3 += min_cost;
					iter_3 += iterations;
					res_aver3 += res_aver;
					res_max3 += res_max;
					res_diff3 += res_diff;

					if (min_cost < res2)
					{
						eff_3++;
						min_cost_c3 += min_cost;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_rand_best1)
				{
					std::cout << "\n DE_RAND_BEST1 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand_best1;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res4, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res4, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res4, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res4, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a4 += min_cost;
					iter_4 += iterations;
					res_aver4 += res_aver;
					res_max4 += res_max;
					res_diff4 += res_diff;

					std::cout << "res_aver = " << res_aver4 << "  ";

					if (min_cost < res2)
					{
						eff_4++;
						min_cost_c4 += min_cost;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_rand_best2)
				{
					std::cout << "\n DE_RAND_BEST2 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand_best2;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res5, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res5, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res5, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res5, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBest2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a5 += min_cost;
					iter_5 += iterations;
					res_aver5 += res_aver;
					res_max5 += res_max;
					res_diff5 += res_diff;

					if (min_cost < res2)
					{
						eff_5++;
						min_cost_c5 += min_cost;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_rand_best_dir1)
				{
					std::cout << "\n DE_RAND_BEST_DIR1 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand_best_dir1;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res6, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBestDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res6, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBestDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res6, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBestDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res6, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandBestDir1Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a6 += min_cost;
					iter_6 += iterations;
					res_aver6 += res_aver;
					res_max6 += res_max;
					res_diff6 += res_diff;


					if (min_cost < res2)
					{
						eff_6++;
						min_cost_c6 += min_cost;
					}

					else
					{
						X.print();
					}
				}

				if (test_de_rand_dir2)
				{
					std::cout << "\n DE_RAND_DIR2 \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_rand_dir2;
					res_aver = 0; res_max = 0, res_diff = 0;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res7, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res7, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res7, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res7, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERandDir2Strategy, Fixed, W, X, Y, V, XAVER, res_aver, res_max, iterations, output);

					res_diff = res_max - min_cost;

					min_cost_a7 += min_cost;
					iter_7 += iterations;
					res_aver7 += res_aver;
					res_max7 += res_max;
					res_diff7 += res_diff;

					if (min_cost < res2)
					{
						eff_7++;
						min_cost_c7 += min_cost;
					}

					else
					{
						X.print();
					}
				}

				std::cout << 100.0 * eff_3 / j;
				std::cout << "(" << iterations << ")  cost=" << min_cost << "\n";

			}


			//End time
			time(&end);
			float time_diff = difftime(end, start);

			//Write results to files
			output = &output_file_de_rand1;
			*output << noise << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\t' << iter_1 << '\t' << it_res1 << '\t' << eff_1 * 100.0 / n_tests << '\t' << res_aver1 << '\t' << res_max1 << '\t' << res_diff1 << '\n';

			output = &output_file_de_rand2;
			*output << noise << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_2 << '\t' << it_res2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << res_aver2 << '\t' << res_max2 << '\t' << res_diff2 << '\n';

			output = &output_file_de_rand_dir1;
			*output << noise << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\t' << iter_3 << '\t' << it_res3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << res_aver3 << '\t' << res_max3 << '\t' << res_diff3 << '\n';

			output = &output_file_de_rand_best1;
			*output << noise << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\t' << iter_4 << '\t' << it_res4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << res_aver4 << '\t' << res_max4 << '\t' << res_diff4 << '\n';

			output = &output_file_de_rand_best2;
			*output << noise << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\t' << iter_5 << '\t' << it_res5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << res_aver5 << '\t' << res_max5 << '\t' << res_diff5 << '\n';

			output = &output_file_de_rand_best_dir1;
			*output << noise << '\t' << min_cost_a6 << '\t' << min_cost_c6 << '\t' << iter_6 << '\t' << it_res6 << '\t' << eff_6 * 100.0 / n_tests << '\t' << res_aver6 << '\t' << res_max6 << '\t' << res_diff6 << '\n';

			output = &output_file_de_rand_dir2;;
			*output << noise << '\t' << min_cost_a7 << '\t' << min_cost_c7 << '\t' << iter_7 << '\t' << it_res7 << '\t' << eff_7 * 100.0 / n_tests << '\t' << res_aver7 << '\t' << res_max7 << '\t' << res_diff7 << '\n';


			//Close files
			output_file_de_rand1.close(); output_file_de_rand2.close(); output_file_de_rand_dir1.close(); output_file_de_rand_dir2.close();
			output_file_de_rand_best1.close(); output_file_de_rand_best2.close(); output_file_de_rand_best_dir1.close();

		}

	}

}


template <typename T>
void CartAnalysis::batchTestDiffEvolutionAdaptiveControl(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, const TMutationStrategy strategy, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Perform batch tests of the 5 differential evolution adaptive control schema
	unsigned short n_tests = 100;
	const unsigned short n_items = nl_test.size();
	const T CR = 0.8;

	unsigned int n_par = 6, m = nl_test.size();
	if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
		n_par = 7;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
		n_par = 5;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
		n_par = 8;

	//Create matrices
	Matrix <double> /*XMIN(1, n_par), XMAX(1, n_par), A2(1, n_par), B2(1, n_par),*/ X(1, n_par), Y(2 * m, 1);

	//Parameters of the genetic algorithm
	const unsigned int population = n_par * n_par, max_gen = 1000;

	//Initialize random number generator
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_de_adapt_ran[256], output_file_text_de_adapt_decr[256], output_file_text_de_mfde[256], output_file_text_de_jitter[256], output_file_text_de_sacp[256];

	strcpy(output_file_text_de_adapt_ran, "de_adapt_ran.log"); strcpy(output_file_text_de_adapt_decr, "de_adapt_decr.log"); strcpy(output_file_text_de_mfde, "de_mfde.log");
	strcpy(output_file_text_de_jitter, "de_jitter.log"); strcpy(output_file_text_de_sacp, "de_sacp.log");

	//Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 6380;

	//New output file variables
	static std::ofstream output_file_de_adapt_decr, output_file_de_adapt_ran, output_file_de_mfde, output_file_de_jitter, output_file_de_sacp;

	//Perform for all noises
	for (int noise = 30; noise <= 50; noise += 20)
	{
		std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

		T delta = 100;
		//delta = 100;
		//for (T delta = 1; delta < 100000; delta *= 10)
		{

			//Open output log files
			output_file_de_adapt_ran.open(output_file_text_de_adapt_ran, std::ofstream::app); output_file_de_adapt_decr.open(output_file_text_de_adapt_decr, std::ofstream::app);
			output_file_de_mfde.open(output_file_text_de_mfde, std::ofstream::app); output_file_de_jitter.open(output_file_text_de_jitter, std::ofstream::app);
			output_file_de_sacp.open(output_file_text_de_sacp, std::ofstream::app);

			//Create variables for measurements
			T min_cost_a1 = 0, min_cost_a2 = 0, min_cost_a3 = 0, min_cost_a4 = 0, min_cost_a5 = 0,
				min_cost_c1 = 0, min_cost_c2 = 0, min_cost_c3 = 0, min_cost_c4 = 0, min_cost_c5 = 0,
				res_avera1 = 0, res_avera2 = 0, res_avera3 = 0, res_avera4 = 0, res_avera5 = 0,
				res_averc1 = 0, res_averc2 = 0, res_averc3 = 0, res_averc4 = 0, res_averc5 = 0,
				res_maxa1 = 0, res_maxa2 = 0, res_maxa3 = 0, res_maxa4 = 0, res_maxa5 = 0,
				res_maxc1 = 0, res_maxc2 = 0, res_maxc3 = 0, res_maxc4 = 0, res_maxc5 = 0,
				res_diffa1 = 0, res_diffa2 = 0, res_diffa3 = 0, res_diffa4 = 0, res_diffa5 = 0,
				res_diffc1 = 0, res_diffc2 = 0, res_diffc3 = 0, res_diffc4 = 0, res_diffc5 = 0;

			unsigned int eff_1 = 0, eff_2 = 0, eff_3 = 0, eff_4 = 0, eff_5 = 0,
				iter_a1 = 0, iter_a2 = 0, iter_a3 = 0, iter_a4 = 0, iter_a5 = 0,
				iter_c1 = 0, iter_c2 = 0, iter_c3 = 0, iter_c4 = 0, iter_c5 = 0,
				it_resa1 = 0, it_resa2 = 0, it_resa3 = 0, it_resa4 = 0, it_resa5 = 0,
				it_resc1 = 0, it_resc2 = 0, it_resc3 = 0, it_resc4 = 0, it_resc5 = 0;

			std::cout << "\n>> >> delta = " << delta << '\n' << '\n';
			//std::cout << "\n>> >> Noise = " << noise << '\n' << '\n';

			//Measure time difference
			time_t start, end;
			time(&start);

			//Perform all tests
			for (unsigned int j = 0; j < n_tests; j++)
			{

				//Lowest and highest values, range
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_a1 = (-noise / 100.0) * 360, high_a1 = (noise / 100.0) * 360;
				T low_a2 = (-noise / 100.0) * 180, high_a2 = (noise / 100.0) * 180;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_a1 = high_a1 - low_a1;
				T range_a2 = high_a2 - low_a2;
				T range_dx = high_dx - low_dx;

				//Min random values
				T rand_R_min = R0 - 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_min = std::min(std::max(-90.0, latp_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_min = std::min(std::max(-180.0, lonp_init - 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_min = std::min(std::max(lat0_min, lat0_init - 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_min = low_dx - 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				//Max random values
				T rand_R_max = R0 + 0.5 * range_R*rand() / (RAND_MAX + 1.0);
				T rand_latp_max = std::min(std::max(-90.0, latp_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), 90.0);
				T rand_lonp_max = std::min(std::max(-180.0, lonp_init + 0.5 * range_a1*rand() / (RAND_MAX + 1.0)), 180.0);
				T rand_lat0_max = std::min(std::max(lat0_min, lat0_init + 0.5 * range_a2*rand() / (RAND_MAX + 1.0)), lat0_max);
				T rand_dx_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy_max = low_dx + 0.5 * range_dx*rand() / (RAND_MAX + 1.0);

				Matrix <T> XMIN(1, n_par), X(1, n_par, 1), XMAX(1, n_par, 1), XAVER(1, n_par), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1);

				//Minimum values
				XMIN(0, 0) = rand_R_min;
				XMIN(0, 1) = rand_latp_min;
				XMIN(0, 2) = rand_lonp_min;
				XMIN(0, 3) = rand_lat0_min;
				XMIN(0, 4) = 0.0;

				//Maximum values
				XMAX(0, 0) = rand_R_max;
				XMAX(0, 1) = rand_latp_max;
				XMAX(0, 2) = rand_lonp_max;
				XMAX(0, 3) = rand_lat0_max;
				XMAX(0, 4) = 0.0;

				//Set intervals
				if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
				{
					XMAX(0, 5) = 10000000;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
				{
					XMIN(0, 5) = 0;
					XMIN(0, 6) = -MAX_LON;

					XMAX(0, 5) = 10000000;
					XMAX(0, 6) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
				{
					XMIN(0, 0) = rand_latp_min;
					XMIN(0, 1) = rand_lonp_min;
					XMIN(0, 2) = rand_lat0_min;
					XMIN(0, 3) = 0;
					XMIN(0, 4) = 0;

					XMAX(0, 0) = rand_latp_max;
					XMAX(0, 1) = rand_lonp_max;
					XMAX(0, 2) = rand_lat0_max;
					XMAX(0, 3) = 0.0;
					XMAX(0, 4) = 10000000;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
				{
					XMIN(0, 5) = -1.0e03;
					XMIN(0, 6) = -1.0e03;
					XMIN(0, 7) = 0;

					XMAX(0, 5) = 1.0e3;
					XMAX(0, 6) = 1.0e3;
					XMAX(0, 7) = 10000000;
				}

				XMIN.print(); XMAX.print();

				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << j << "/" << n_tests << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 1.0e2;

				std::cout << "\n Adapt_random \n";
				output = &output_file_de_adapt_ran;

				T min_cost = 0, res_aver = 0, res_max = 0, res_diff = 0;
				T q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
				unsigned int it_res = 0;

				bool test_de_adapt_ran = true, test_de_adapt_decr = true, test_de_mfde = true, test_de_jitter = true, test_de_sacp = true;
				//bool test_de_adapt_ran = false, test_de_adapt_decr = false, test_de_mfde = true, test_de_jitter = false, test_de_sacp = false;

				if (test_de_adapt_ran)
				{
					it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveRandom, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveRandom, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveRandom, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveRandom, W, X, Y, V, XAVER, res_aver, res_max, iterations);

					res_diff = res_max - min_cost;

					min_cost_a1 += min_cost;
					iter_a1 += iterations;
					it_resa1 += it_res;
					res_avera1 += res_aver;
					res_maxa1 += res_max;
					res_diffa1 += res_diff;

					if (min_cost < res2)
					{
						eff_1++;
						min_cost_c1 += min_cost;
						iter_c1 += iterations;
						it_resc1 += it_res;
						res_averc1 += res_aver;
						res_maxc1 += res_max;
						res_diffc1 += res_diff;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_adapt_decr)
				{
					std::cout << "\n Adapt_decr \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_adapt_decr;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveDecreasing, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveDecreasing, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveDecreasing, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, AdaptiveDecreasing, W, X, Y, V, XAVER, res_aver, res_max, iterations);

					res_diff = res_max - min_cost;

					min_cost_a2 += min_cost;
					iter_a2 += iterations;
					it_resa2 += it_res;
					res_avera2 += res_aver;
					res_maxa2 += res_max;
					res_diffa2 += res_diff;

					if (min_cost < res2)
					{
						eff_2++;
						min_cost_c2 += min_cost;
						iter_c2 += iterations;
						it_resc2 += it_res;
						res_averc2 += res_aver;
						res_maxc2 += res_max;
						res_diffc2 += res_diff;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_mfde)
				{
					std::cout << "\n MFDE \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_mfde;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, MFDE, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, MFDE, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, MFDE, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, MFDE, W, X, Y, V, XAVER, res_aver, res_max, iterations);

					res_diff = res_max - min_cost;

					min_cost_a3 += min_cost;
					iter_a3 += iterations;
					it_resa3 += it_res;
					res_avera3 += res_aver;
					res_maxa3 += res_max;
					res_diffa3 += res_diff;

					if (min_cost < res2)
					{
						eff_3++;
						min_cost_c3 += min_cost;
						iter_c3 += iterations;
						it_resc3 += it_res;
						res_averc3 += res_aver;
						res_maxc3 += res_max;
						res_diffc3 += res_diff;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_jitter)
				{
					std::cout << "\n JITTER \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_jitter;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, n_par);

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, Jitter, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, Jitter, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, Jitter, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, Jitter, W, X, Y, V, XAVER, res_aver, res_max, iterations);

					res_diff = res_max - min_cost;

					min_cost_a4 += min_cost;
					iter_a4 += iterations;
					it_resa4 += it_res;
					res_avera4 += res_aver;
					res_maxa4 += res_max;
					res_diffa4 += res_diff;


					if (min_cost < res2)
					{
						eff_4++;
						min_cost_c4 += min_cost;
						iter_c4 += iterations;
						it_resc4 += it_res;
						res_averc4 += res_aver;
						res_maxc4 += res_max;
						res_diffc4 += res_diff;
					}

					else
					{
						X.print();
					}
				}


				if (test_de_sacp)
				{
					std::cout << "\n SACP \n";
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);
					output = &output_file_de_sacp;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, SACP, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV3DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, XAVER, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, SACP, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, SACP, W, X, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjVDE <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, population, eps, max_gen, F, CR, strategy, SACP, W, X, Y, V, XAVER, res_aver, res_max, iterations);

					res_diff = res_max - min_cost;

					min_cost_a5 += min_cost;
					iter_a5 += iterations;
					it_resa5 += it_res;
					res_avera5 += res_aver;
					res_maxa5 += res_max;
					res_diffa5 += res_diff;


					if (min_cost < res2)
					{
						eff_5++;
						min_cost_c5 += min_cost;
						iter_c5 += iterations;
						it_resc5 += it_res;
						res_averc5 += res_aver;
						res_maxc5 += res_max;
						res_diffc5 += res_diff;
					}

					else
					{
						X.print();
					}
				}


				std::cout << 100.0 * eff_3 / j;
				std::cout << "(" << iterations << ")  cost=" << min_cost << "\n";

			}


			//End time
			time(&end);
			float time_diff = difftime(end, start);

			//Write results to files
			output = &output_file_de_adapt_ran;
			*output << noise << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\t' << iter_a1 << '\t' << iter_c1 << '\t' << it_resa1 << '\t' << it_resc1 << '\t' << eff_1 * 100.0 / n_tests << '\t' << res_avera1 << '\t' << res_averc1 << '\t' << res_maxa1 << '\t' << res_maxc1 << '\t' << res_diffa1 << '\t' << res_diffc1 << '\n';

			output = &output_file_de_adapt_decr;
			*output << noise << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_a2 << '\t' << iter_c2 << '\t' << it_resa2 << '\t' << it_resc2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << res_avera2 << '\t' << res_averc2 << '\t' << res_maxa2 << '\t' << res_maxc2 << '\t' << res_diffa2 << '\t' << res_diffc2 << '\n';

			output = &output_file_de_mfde;
			*output << noise << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\t' << iter_a3 << '\t' << iter_c3 << '\t' << it_resa3 << '\t' << it_resc3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << res_avera3 << '\t' << res_averc3 << '\t' << res_maxa3 << '\t' << res_maxc3 << '\t' << res_diffa3 << '\t' << res_diffc3 << '\n';

			output = &output_file_de_jitter;
			*output << noise << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\t' << iter_a4 << '\t' << iter_c4 << '\t' << it_resa4 << '\t' << it_resc4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << res_avera4 << '\t' << res_averc4 << '\t' << res_maxa4 << '\t' << res_maxc4 << '\t' << res_diffa4 << '\t' << res_diffc4 << '\n';;

			output = &output_file_de_sacp;
			*output << noise << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\t' << iter_a5 << '\t' << iter_c5 << '\t' << it_resa5 << '\t' << it_resc5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << res_avera5 << '\t' << res_averc5 << '\t' << res_maxa5 << '\t' << res_maxc5 << '\t' << res_diffa5 << '\t' << res_diffc5 << '\n';;

			//Close files
			output_file_de_adapt_ran.close(); output_file_de_adapt_decr.close(); output_file_de_mfde.close();
			output_file_de_jitter.close(); output_file_de_sacp.close();

		}

	}
}

template <typename T>
void CartAnalysis::batchTestDiffEvolutionOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init,
	const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, const TMutationStrategy strategy, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Outlier test for differential evolution
	unsigned short n_tests = 100;
	const unsigned short n_items = nl_test.size();

	//Divide P
	for (unsigned int i = 0; i < n_items; i++)
	{
		nl_test[i]->setX(nl_test[i]->getX() / 100000);
		nl_test[i]->setY(nl_test[i]->getY() / 100000);
	}

	//Max coordinate error
	const T sx = 5.0 / 3 / sqrt(2.0) / 1000;

	/*
	for (unsigned int i = 0; i < n_items; i++)
	{
	nl_test[i]->setX(nl_test[i]->getX());
	nl_test[i]->setY(nl_test[i]->getY());
	}

	//Max coordinate error
	const T sx = 5.0 / 3 / sqrt(2.0) * 100;
	*/
	const T sy = sx;
	const T slat = 4.5 / 3.0 / sqrt(2.0);
	const T slon = slat;

	unsigned int n_par = 6, m = nl_test.size();
	if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
		n_par = 7;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
		n_par = 5;
	else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
		n_par = 8;

	//Create matrices
	Matrix <double> /*XMIN(1, n_par), XMAX(1, n_par), A2(1, n_par), B2(1, n_par),*/ X(1, n_par), Y(2 * m, 1);

	//Parameters of the genetic algorithm
	const unsigned int population = n_par * n_par, max_gen = 500;
	const T CR = 0.8;

	//Initialize random number generator
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_de_rand1[256], output_file_text_de_rand2[256], output_file_text_de_rand3[256], output_file_text_de_rand4[256], output_file_text_de_rand5[256];

	strcpy(output_file_text_de_rand1, "de_rand_hub.log"); strcpy(output_file_text_de_rand2, "de_rand_andrew.log"); strcpy(output_file_text_de_rand3, "de_rand_danish.log"); strcpy(output_file_text_de_rand4, "de_rand_yang.log"); strcpy(output_file_text_de_rand5, "de_rand_tukey.log");

	//Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 6380;

	T mult1 = 1, mult2 = 0.5, mult3 = 0.0001;
	unsigned int mode = 2;

	//New output file variables
	static std::ofstream output_file_de_rand1, output_file_de_rand2, output_file_de_rand3, output_file_de_rand4, output_file_de_rand5;

	//Initialize random number generator
	srand(time(NULL));

	for (T k = 2.0; k <= 3; k += 0.5)
	{
		//Open output log file
		output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app);

		std::cout << " >>>>>>>>>>>>> k = " << k << " <<<<<<<<<<<<<<<<<<" << '\n';

		//Write to file
		output = &output_file_de_rand1;
		*output << "\n K = " << k << '\n';

		//Close files
		output_file_de_rand1.close();

		//Process for all noise levels
		unsigned int percs[] = { 10, 20, 30, 40, 45, 50 };
		for (unsigned int i_per = 0; i_per < 6; i_per++)
		{

			//Open output log files
			output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app), output_file_de_rand2.open(output_file_text_de_rand2, std::ofstream::app), output_file_de_rand3.open(output_file_text_de_rand3, std::ofstream::app), output_file_de_rand4.open(output_file_text_de_rand4, std::ofstream::app), output_file_de_rand5.open(output_file_text_de_rand5, std::ofstream::app);

			T min_cost_a_rand1 = 0, min_cost_a_rand2 = 0, min_cost_a_rand3 = 0, min_cost_a_rand4 = 0, min_cost_a_rand5 = 0,
				min_cost_c_rand1 = 0, min_cost_c_rand2 = 0, min_cost_c_rand3 = 0, min_cost_c_rand4 = 0, min_cost_c_rand5 = 0,
				res_aver_a_rand1 = 0, res_aver_a_rand2 = 0, res_aver_a_rand3 = 0, res_aver_a_rand4 = 0, res_aver_a_rand5 = 0,
				res_aver_c_rand1 = 0, res_aver_c_rand2 = 0, res_aver_c_rand3 = 0, res_aver_c_rand4 = 0, res_aver_c_rand5 = 0,
				res_max_a_rand1 = 0, res_max_a_rand2 = 0, res_max_a_rand3 = 0, res_max_a_rand4 = 0, res_max_a_rand5 = 0,
				res_max_c_rand1 = 0, res_max_c_rand2 = 0, res_max_c_rand3 = 0, res_max_c_rand4 = 0, res_max_c_rand5 = 0,
				res_diff_a_rand1 = 0, res_diff_a_rand2 = 0, res_diff_a_rand3 = 0, res_diff_a_rand4 = 0, res_diff_a_rand5 = 0,
				res_diff_c_rand1 = 0, res_diff_c_rand2 = 0, res_diff_c_rand3 = 0, res_diff_c_rand4 = 0, res_diff_c_rand5 = 0,
				total_efficiency_rand1 = 0, total_efficiency_rand2 = 0, total_efficiency_rand3 = 0, total_efficiency_rand4 = 0, total_efficiency_rand5 = 0,
				dx_rand_a1 = 0, dx_rand_a2 = 0, dx_rand_a3 = 0, dx_rand_a4 = 0, dx_rand_a5 = 0,
				dx_rand_c1 = 0, dx_rand_c2 = 0, dx_rand_c3 = 0, dx_rand_c4 = 0, dx_rand_c5 = 0,
				dx_rand_a11 = 0, dx_rand_a21 = 0, dx_rand_a31 = 0, dx_rand_a41 = 0, dx_rand_a51 = 0,
				dx_rand_c11 = 0, dx_rand_c21 = 0, dx_rand_c31 = 0, dx_rand_c41 = 0, dx_rand_c51 = 0,
				min_cost_a_rand11 = 0, min_cost_c_rand11 = 0, diff1 = 0, diff2 = 0;

			unsigned int eff_rand1 = 0, eff_rand2 = 0, eff_rand3 = 0, eff_rand4 = 0, eff_rand5 = 0,
				eff_rand11 = 0, eff_rand21 = 0, eff_rand31 = 0, eff_rand41 = 0, eff_rand51 = 0,
				iter_a_rand1 = 0, iter_a_rand2 = 0, iter_a_rand3 = 0, iter_a_rand4 = 0, iter_a_rand5 = 0,
				iter_c_rand1 = 0, iter_c_rand2 = 0, iter_c_rand3 = 0, iter_c_rand4 = 0, iter_c_rand5 = 0,
				it_res_a_rand1 = 0, it_res_a_rand2 = 0, it_res_a_rand3 = 0, it_res_a_rand4 = 0, it_res_a_rand5 = 0,
				it_res_c_rand1 = 0, it_res_c_rand2 = 0, it_res_c_rand3 = 0, it_res_c_rand4 = 0, it_res_c_rand5 = 0,
				iter_a_rand11 = 0, iter_c_rand11 = 0, it_res_a_rand11 = 0, it_res_c_rand11 = 0;

			unsigned int efficiency_rand1 = 0, efficiency_rand2 = 0, efficiency_rand3 = 0, efficiency_rand4 = 0, efficiency_rand5 = 0,
				false_out_rand1 = 0, false_out_rand2 = 0, false_out_rand3 = 0, false_out_rand4 = 0, false_out_rand5 = 0;

			const T perc = percs[i_per];
			//const T perc = 30;

			std::cout << " *********** perc = " << perc << "************* \n";

			//Repeat for all tests
			for (unsigned int i = 0; i < n_tests; i++)
			{
				Container <Node3DCartesian <T> *> nl_test_out = nl_test;
				Container <Point3DGeographic <T> *> pl_reference_out = pl_reference;

				std::cout << "Test " << i << "/" << n_tests << '\n';

				//Create random permutation
				unsigned int n_rand = (perc * n_items) / 100;
				//n_rand = 6;

				Matrix <unsigned int> I = RandomPermutation::randperm(n_items, n_rand);

				/*
				Matrix <unsigned int> I(1, 6);

				I(0, 0) = 0;
				I(0, 1) = 2;
				I(0, 2) = 4;
				I(0, 3) = 6;
				I(0, 4) = 8;
				I(0, 5) = 10;
				*/
				//Create initial weight matrix
				Matrix <T> WI = ones(2 * n_items, 2 * n_items, 1.0);

				for (unsigned int j = 0; j < n_rand; j++)
				{
					WI(I(0, j), I(0, j)) = 0.0;
					WI(I(0, j) + n_items, I(0, j) + n_items) = 0.0;
				}
				/*
				nl_test_out[I(0, 0)]->setX(nl_test[I(0, 0)]->getX() + 3 * sx);
				nl_test_out[I(0, 0)]->setY(nl_test[I(0, 0)]->getY() - 3 * sy);
				nl_test_out[I(0, 1)]->setX(nl_test[I(0, 1)]->getX() + 4 * sx);
				nl_test_out[I(0, 1)]->setY(nl_test[I(0, 1)]->getY() - 4 * sy);
				nl_test_out[I(0, 2)]->setX(nl_test[I(0, 2)]->getX() + 4 * sx);
				nl_test_out[I(0, 2)]->setY(nl_test[I(0, 2)]->getY() - 4 * sy);
				nl_test_out[I(0, 3)]->setX(nl_test[I(0, 3)]->getX() + 5 * sx);
				nl_test_out[I(0, 3)]->setY(nl_test[I(0, 3)]->getY() - 5 * sy);
				nl_test_out[I(0, 4)]->setX(nl_test[I(0, 4)]->getX() - 3 * sx);
				nl_test_out[I(0, 4)]->setY(nl_test[I(0, 4)]->getY() + 3 * sy);
				nl_test_out[I(0, 5)]->setX(nl_test[I(0, 5)]->getX() - 4 * sx);
				nl_test_out[I(0, 5)]->setY(nl_test[I(0, 5)]->getY() + 4 * sy);

				*/

				//Contaminated randomly selected points on P: outliers
				for (unsigned int j = 0; j < n_rand; j++)
				{
					//Get coordinates
					T x = nl_test[I(0, j)]->getX();
					T y = nl_test[I(0, j)]->getY();

					//Get random numbers
					const T r1 = ((T)rand() / (RAND_MAX));
					const T r2 = ((T)rand() / (RAND_MAX));

					//Get sign
					const T sig_x = copysign(1.0, r1 - 0.5);
					const T sig_y = copysign(1.0, r2 - 0.5);

					const T r3 = ((T)rand() / (RAND_MAX));
					const T r4 = ((T)rand() / (RAND_MAX));

					//Get random numbers
					const T dx = 3 + r3 * 6;
					const T dy = 3 + r4 * 6;

					//Errors
					const T deltax = sig_x * dx * sx;
					const T deltay = sig_y * dy * sy;

					//Add errors
					x = x + mult1 * deltax;
					y = y + mult1 * deltay;

					//Create new points
					Node3DCartesian <T> *p = new Node3DCartesian <T>(x, y);

					//Update coordinates
					nl_test_out[I(0, j)]->setX(x);
					nl_test_out[I(0, j)]->setY(y);
				}


				//Contaminated P: random errors, all items
				if (mode == 2)
				{
					for (unsigned int j = 0; j < n_items; j++)
					{
						//Get cooridnates
						T x1 = nl_test_out[j]->getX();
						T y1 = nl_test_out[j]->getY();
						T lat1 = pl_reference_out[j]->getLat();
						T lon1 = pl_reference_out[j]->getLon();

						//Get random numbers
						const T r1 = ((T)rand() / (RAND_MAX));
						const T r2 = ((T)rand() / (RAND_MAX));
						const T r3 = ((T)rand() / (RAND_MAX));
						const T r4 = ((T)rand() / (RAND_MAX));

						//Get signs
						const T sig_x = copysign(1.0, r1 - 0.5);
						const T sig_y = copysign(1.0, r2 - 0.5);
						const T sig_lat = copysign(1.0, r3 - 0.5);
						const T sig_lon = copysign(1.0, r4 - 0.5);

						//Get random numbers
						const T dx = ((T)rand() / (RAND_MAX));
						const T dy = ((T)rand() / (RAND_MAX));
						const T dlat = ((T)rand() / (RAND_MAX));
						const T dlon = ((T)rand() / (RAND_MAX));

						//Errors
						const T deltax = sig_x * dx * sx;
						const T deltay = sig_y * dy * sy;
						const T deltalat = sig_lat * dlat * slat; //
						const T deltalon = sig_lon * dlon * slon; //

						//slon multiplier
						const T m = 1.0 / cos(lat1 * M_PI / 180);

						//Add errors
						x1 = x1 + mult2 * deltax;
						y1 = y1 + mult2 * deltay;
						lat1 = lat1 + mult2 * deltalat;
						lon1 = lon1 + m * mult2 * deltalon;

						//Update first point
						nl_test_out[j]->setX(x1);
						nl_test_out[j]->setY(y1);

						//Update second point
						pl_reference_out[j]->setLat(lat1);
						pl_reference_out[j]->setLon(lon1);
					}
				}


				// Lowest and highest values, range
				T noise = 50;
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_dx = high_dx - low_dx;

				//Random values
				T rand_R = (low_R + range_R*rand() / (RAND_MAX + 1.0));
				T rand_latp = latp_min + (latp_max - latp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lonp = lonp_min + (lonp_max - lonp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lat0 = lat0_min + (lat0_max - lat0_min) * rand() / (RAND_MAX + 1.0);
				T rand_dx = low_dx + range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy = low_dx + range_dx*rand() / (RAND_MAX + 1.0);

				Matrix <T> XMIN(1, n_par), X(1, n_par, 1), XMAX(1, n_par, 1), XAVER(1, n_par), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1),
					XEQDC(1, n_par), XLAEA(1, n_par), XMERC(1, n_par), XOK(1, n_par, 1), DXT(1, n_par);

				//Store randomly generated values
				X(0, 0) = rand_R;
				X(0, 1) = rand_latp;
				X(0, 2) = rand_lonp;
				X(0, 3) = rand_lat0;

				//Set intervals
				XMIN(0, 0) = 0; XMAX(0, 0) = low_R + range_R;
				XMIN(0, 1) = latp_min; XMAX(0, 1) = latp_max;
				XMIN(0, 2) = lonp_min; XMAX(0, 2) = lonp_max;
				XMIN(0, 3) = lat0_min; XMAX(0, 3) = lat0_max;
				XMIN(0, 4) = 0.0; XMAX(0, 4) = 0.0;

				//Correct values
				XEQDC(0, 0) = 0.06380;
				XEQDC(0, 1) = 90;
				XEQDC(0, 2) = 0;
				XEQDC(0, 3) = 50;
				XEQDC(0, 4) = 0.0;

				XLAEA(0, 0) = 0.06380;
				XLAEA(0, 1) = 52;
				XLAEA(0, 2) = 10;
				XLAEA(0, 3) = 0.0;
				XLAEA(0, 4) = 0.0;

				XMERC(0, 0) = 0.06380;
				XMERC(0, 1) = 0;
				XMERC(0, 2) = 90;
				XMERC(0, 3) = 47.0;
				XMERC(0, 4) = 0.0;

				//Threshold
				DXT(0, 0) = 6380;
				DXT(0, 1) = 180;
				DXT(0, 2) = 360;
				DXT(0, 3) = 180;

				//Set intervals
				if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
				{
					XMAX(0, 5) = 0;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod)
				{
					XMIN(0, 5) = 0; XMAX(0, 5) = 0;
					XMIN(0, 6) = -MAX_LON; XMAX(0, 6) = MAX_LON;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
				{
					X(0, 0) = rand_latp;
					X(0, 1) = rand_lonp;
					X(0, 2) = rand_lat0;
					X(0, 3) = 0;
					X(0, 4) = 1;

					XMIN(0, 0) = latp_min; XMAX(0, 0) = latp_max;
					XMIN(0, 1) = lonp_min; XMAX(0, 1) = lonp_max;
					XMIN(0, 2) = lat0_min; XMAX(0, 2) = lat0_max;
					XMIN(0, 3) = 0.0; XMAX(0, 3) = 0.0;
					XMIN(0, 4) = 0.0; XMAX(0, 4) = 0.0;

					//Correct values
					XEQDC(0, 0) = 90;
					XEQDC(0, 1) = 0;
					XEQDC(0, 2) = 50;
					XEQDC(0, 3) = 0.0;

					XLAEA(0, 0) = 52;
					XLAEA(0, 1) = 10;
					XLAEA(0, 2) = 0.0;
					XLAEA(0, 3) = 0.0;

					XMERC(0, 0) = 0;
					XMERC(0, 1) = 90;
					XMERC(0, 2) = 47.0;
					XMERC(0, 3) = 0.0;

					DXT(0, 0) = 180;
					DXT(0, 1) = 360;
					DXT(0, 2) = 180;
				}

				else if (analysis_parameters.analysis_method == DifferentialEvolutionShiftsMethod)
				{
					//XMIN(0, 5) = -1.0e03;
					//XMIN(0, 6) = -1.0e03;
					XMIN(0, 7) = 1;

					XMIN(0, 5) = -1.0e09; XMAX(0, 5) = 1.0e9;
					XMIN(0, 6) = -1.0e09; XMAX(0, 6) = 1.0e9;
					XMIN(0, 7) = 0.0; XMAX(0, 7) = 10000000;
				}

				//Update threshold
				DXT = DXT * 1.0 / 200;

				//XMIN.print(); XMAX.print();

				//Print actual values
				//std::cout << "Noise = " << noise << '\t' << "Test " << j << "/" << n_tests << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				//const T res2 = 1.0e2;

				const T res2 = 5.0e-4;

				//Get projection ID
				const char * proj_text = proj->getName();
				if (strcmp(proj_text, "eqdc") == 0)
				{
					XOK = XEQDC;
				}

				else if (strcmp(proj_text, "laea") == 0)
				{
					XOK = XLAEA;
				}

				else if (strcmp(proj_text, "merc") == 0)
				{
					XOK = XMERC;
				}

				//T res2 = 0;
				//Rotated case
				//if (analysis_parameters.analysis_method == DifferentialEvolutionRotMethod || analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)

				T min_cost = 0, res_aver = 0, res_max = 0, res_diff = 0;
				unsigned int it_res = 0;

				bool test_huber = true, test_andrew = true, test_danish = true, test_yang = true, test_tukey = true;

				//DE / RAND 1, HUBER function
				//**********************************************************************************************************************************
				if (test_huber)
				{
					//DE\RAND\BEST1
					std::cout << "\n Huber  \n";
					T q1 = 1, q2 = 1, R_def = 1;
					output = &output_file_de_rand1;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					analysis_parameters.remove_outliers = true;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand1 += efficiency;
					false_out_rand1 += false_out;
					total_efficiency_rand1 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					Matrix <T> WID = diag(WI);
					Matrix <T> WD = diag(W);
					//WD.print();

					res_diff = res_max - min_cost;

					min_cost_a_rand1 += min_cost;
					iter_a_rand1 += iterations;
					it_res_a_rand1 += it_res;
					res_aver_a_rand1 += res_aver;
					res_max_a_rand1 += res_max;
					res_diff_a_rand1 += res_diff;

					//Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a1 += norm(DX * trans(DX));

					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{

						eff_rand1++;
						min_cost_c_rand1 += min_cost;
						iter_c_rand1 += iterations;
						it_res_c_rand1 += it_res;
						res_aver_c_rand1 += res_aver;
						res_max_c_rand1 += res_max;
						res_diff_c_rand1 += res_diff;
						dx_rand_c1 += norm(DX * trans(DX));

					}

					else
					{
						XX.print();
					}


					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; it_res = 0;

					analysis_parameters.remove_outliers = false;


					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);


					if (n_par == 5) X2(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						X2(0, 1) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					diff2 = norm(DX2 * trans(DX2));
					dx_rand_a11 += norm(DX2 * trans(DX2));
					min_cost_a_rand11 += min_cost;
					iter_a_rand11 += iterations;
					it_res_a_rand11 += it_res;

					diff1 = norm(trans(DX) * DX);

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					if (detected)
					{
						eff_rand11++;
						min_cost_c_rand11 += min_cost;
						dx_rand_c11 += norm(DX2 * trans(DX2));
						iter_c_rand11 += iterations;
						it_res_c_rand11 += it_res;
					}

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a1 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c1 << "\t DNWSC = " << dx_rand_c11 << "\t EFFW = " << eff_rand1 * 100.0 / (i + 1) << "\t EFFNW = " << eff_rand11 * 100.0 / (i + 1) << '\n';

					///*********************************************************
					XX.print(); X2.print();

				}


				//DE / RAND 1, ANDREW function
				//**********************************************************************************************************************************
				if (test_andrew)
				{
					//DE\RAND\BEST1
					std::cout << "\n Andrew 1 \n";
					T q1 = 1, q2 = 1, R_def = 1;
					output = &output_file_de_rand1;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					analysis_parameters.remove_outliers = true;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, IX, output), XMIN, XMAX, population, mult3  * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand2 += efficiency;
					false_out_rand2 += false_out;
					total_efficiency_rand2 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					Matrix <T> WID = diag(WI);
					Matrix <T> WD = diag(W);
					//WD.print();

					res_diff = res_max - min_cost;

					min_cost_a_rand2 += min_cost;
					iter_a_rand2 += iterations;
					it_res_a_rand2 += it_res;
					res_aver_a_rand2 += res_aver;
					res_max_a_rand2 += res_max;
					res_diff_a_rand2 += res_diff;

					//Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a2 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{

						eff_rand2++;
						min_cost_c_rand2 += min_cost;
						iter_c_rand2 += iterations;
						it_res_c_rand2 += it_res;
						res_aver_c_rand2 += res_aver;
						res_max_c_rand2 += res_max;
						res_diff_c_rand2 += res_diff;
						dx_rand_c2 += norm(DX * trans(DX));

					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = ones(1, n_par, 1.0);
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, IX, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);

					if (n_par == 5) X2(0, 4) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(0, 1) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a21 += norm(DX2 * trans(DX2));

					T diff1 = norm(trans(DX) * DX);

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					//if (min_cost < res2)
					if (detected)
					{
					eff_rand21++;
					dx_rand_c21 += norm(DX2 * trans(DX2));
					}
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a2 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c2 << "\t DNWSC = " << dx_rand_c11 << '\n';

					//********************************************************
					XX.print();
					//X2.print();
				}


				//DE / RAND 1, DANISH function
				//**********************************************************************************************************************************
				if (test_danish)
				{
					//DE\RAND\BEST1
					std::cout << "\n Danish \n";
					T q1 = 1, q2 = 1, R_def = 1;
					output = &output_file_de_rand1;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					analysis_parameters.remove_outliers = true;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand3 += efficiency;
					false_out_rand3 += false_out;
					total_efficiency_rand3 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					Matrix <T> WID = diag(WI);
					Matrix <T> WD = diag(W);
					//WD.print();

					res_diff = res_max - min_cost;

					min_cost_a_rand3 += min_cost;
					iter_a_rand3 += iterations;
					it_res_a_rand3 += it_res;
					res_aver_a_rand3 += res_aver;
					res_max_a_rand3 += res_max;
					res_diff_a_rand3 += res_diff;

					//Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a3 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{

						eff_rand3++;
						min_cost_c_rand3 += min_cost;
						iter_c_rand3 += iterations;
						it_res_c_rand3 += it_res;
						res_aver_c_rand3 += res_aver;
						res_max_c_rand3 += res_max;
						res_diff_c_rand3 += res_diff;
						dx_rand_c3 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = ones(1, n_par, 1.0);
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, IX, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);

					if (n_par == 5) X2(0, 4) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(0, 1) = 0;
					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a31 += norm(DX2 * trans(DX2));

					T diff1 = norm(trans(DX) * DX);

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					//if (min_cost < res2)
					if (detected)
					{
					eff_rand31++;
					dx_rand_c31 += norm(DX2 * trans(DX2));
					}
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a3 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c3 << "\t DNWSC = " << dx_rand_c11 << '\n';
					//********************************************************

					XX.print();
					//X2.print();
				}


				//DE / RAND 1, YANG function
				//**********************************************************************************************************************************
				if (test_yang)
				{
					//DE\RAND\BEST1
					std::cout << "\n Yang \n";
					T q1 = 1, q2 = 1, R_def = 1;
					output = &output_file_de_rand1;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					analysis_parameters.remove_outliers = true;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, YangFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, YangFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand4 += efficiency;
					false_out_rand4 += false_out;
					total_efficiency_rand4 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					Matrix <T> WID = diag(WI);
					Matrix <T> WD = diag(W);
					//WD.print();

					res_diff = res_max - min_cost;

					min_cost_a_rand4 += min_cost;
					iter_a_rand4 += iterations;
					it_res_a_rand4 += it_res;
					res_aver_a_rand4 += res_aver;
					res_max_a_rand4 += res_max;
					res_diff_a_rand4 += res_diff;

					//Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a4 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_rand4++;
						min_cost_c_rand4 += min_cost;
						iter_c_rand4 += iterations;
						it_res_c_rand4 += it_res;
						res_aver_c_rand4 += res_aver;
						res_max_c_rand4 += res_max;
						res_diff_c_rand4 += res_diff;
						dx_rand_c4 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = ones(1, n_par, 1.0);
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, YangFunction, k, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, YangFunction, k, IX, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);

					if (n_par == 5) X2(0, 4) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(0, 1) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a41 += norm(DX2 * trans(DX2));

					T diff1 = norm(trans(DX) * DX);

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					//if (min_cost < res2)
					if (detected)
					{
					eff_rand41++;
					dx_rand_c41 += norm(DX2 * trans(DX2));
					}
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a4 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c4 << "\t DNWSC = " << dx_rand_c11 << '\n';

					//********************************************************

					XX.print();
					//X2.print();
				}


				//DE / RAND 1, TUKEY function
				//**********************************************************************************************************************************
				if (test_tukey)
				{
					//DE\RAND\BEST1
					std::cout << "\n Tukey \n";
					T q1 = 1, q2 = 1, R_def = 1;
					output = &output_file_de_rand1;
					res_aver = 0; res_max = 0, res_diff = 0, it_res = 0;
					Matrix <T> F(1, 1); F(0, 0) = 0.5;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					analysis_parameters.remove_outliers = true;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
						min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, IX, output), XMIN, XMAX, population, mult3 * eps, max_gen, F, CR, DERand1Strategy, Fixed, W, XX, Y, V, XAVER, res_aver, res_max, iterations);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand5 += efficiency;
					false_out_rand5 += false_out;
					total_efficiency_rand5 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					Matrix <T> WID = diag(WI);
					Matrix <T> WD = diag(W);
					//WD.print();

					res_diff = res_max - min_cost;

					min_cost_a_rand5 += min_cost;
					iter_a_rand5 += iterations;
					it_res_a_rand5 += it_res;
					res_aver_a_rand5 += res_aver;
					res_max_a_rand5 += res_max;
					res_diff_a_rand5 += res_diff;

					//Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a5 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));
					//if (min_cost < res2)
					if (detected)
					{

						eff_rand5++;
						min_cost_c_rand5 += min_cost;
						iter_c_rand5 += iterations;
						it_res_c_rand5 += it_res;
						res_aver_c_rand5 += res_aver;
						res_max_c_rand5 += res_max;
						res_diff_c_rand5 += res_diff;
						dx_rand_c5 += norm(DX * trans(DX));

					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = ones(1, n_par, 1.0);
					q1 = 1, q2 = 1, R_def = 0.5 * (rand_R_min + rand_R_max);

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == DifferentialEvolutionMethod)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV2DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);
					else if (analysis_parameters.analysis_method == DifferentialEvolutionRot2Method)
					min_cost = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DE <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
					total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, IX, output), XMIN, XMAX, population, eps, max_gen, F, CR, DERand1Strategy, Fixed, W2, X2, Y, V, XAVER, res_aver, res_max, iterations);

					if (n_par == 5) X2(0, 4) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(0, 1) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a51 += norm(DX2 * trans(DX2));

					T diff1 = norm(trans(DX) * DX);

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					//if (min_cost < res2)
					if (detected)
					{
					eff_rand51++;
					dx_rand_c51 += norm(DX2 * trans(DX2));
					}

					*/
					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a5 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c5 << "\t DNWSC = " << dx_rand_c11 << '\n';

					//********************************************************

					XX.print();
					//X2.print();
				}
			}

			//End time
			//time(&end);
			//float time_diff = difftime(end, start);

			//Write results to files
			output = &output_file_de_rand1;
			*output << perc << '\t' << efficiency_rand1 / n_tests << '\t' << false_out_rand1 / n_tests << '\t' << total_efficiency_rand1 / n_tests << '\t' << dx_rand_a1 << '\t' << dx_rand_c1 << '\t' << iter_a_rand1 << '\t' << iter_c_rand1 << '\t' << it_res_a_rand1 << '\t' << it_res_c_rand1 << '\t' << eff_rand1 * 100.0 / n_tests << '\t' << min_cost_a_rand1 << '\t' << min_cost_c_rand1 << '\n';
			*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << iter_a_rand2 << '\t' << iter_c_rand2 << '\t' << it_res_a_rand2 << '\t' << it_res_c_rand2 << '\t' << eff_rand2 * 100.0 / n_tests << '\t' << min_cost_a_rand2 << '\t' << min_cost_c_rand2 << '\n';
			*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << iter_a_rand3 << '\t' << iter_c_rand3 << '\t' << it_res_a_rand3 << '\t' << it_res_c_rand3 << '\t' << eff_rand3 * 100.0 / n_tests << '\t' << min_cost_a_rand3 << '\t' << min_cost_c_rand3 << '\n';
			*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << iter_a_rand4 << '\t' << iter_c_rand4 << '\t' << it_res_a_rand4 << '\t' << it_res_c_rand4 << '\t' << eff_rand4 * 100.0 / n_tests << '\t' << min_cost_a_rand4 << '\t' << min_cost_c_rand4 << '\n';
			*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << iter_a_rand5 << '\t' << iter_c_rand5 << '\t' << it_res_a_rand5 << '\t' << it_res_c_rand5 << '\t' << eff_rand5 * 100.0 / n_tests << '\t' << min_cost_a_rand5 << '\t' << min_cost_c_rand5 << '\n';
			*output << perc << '\t' << "0" << '\t' << "0" << '\t' << "0" << '\t' << dx_rand_a11 << '\t' << dx_rand_c11 << '\t' << iter_a_rand11 << '\t' << iter_c_rand11 << '\t' << it_res_a_rand11 << '\t' << it_res_c_rand11 << '\t' << eff_rand11 * 100.0 / n_tests << '\t' << min_cost_a_rand11 << '\t' << min_cost_c_rand11 << '\n';

			//output = &output_file_de_rand2;
			//*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << dx_rand_a21 << '\t' << dx_rand_c21 << '\t' << iter_a_rand2 << '\t' << iter_c_rand2 << '\t' << it_res_a_rand2 << '\t' << it_res_c_rand2 << '\t' << eff_rand2 * 100.0 / n_tests << '\t' << eff_rand21 * 100.0 / n_tests << '\t' << res_aver_a_rand2 << '\t' << res_aver_c_rand2 << '\t' << res_max_a_rand2 << '\t' << res_max_c_rand2 << '\t' << res_diff_a_rand2 << '\t' << res_diff_c_rand2 << '\n';

			//output = &output_file_de_rand3;
			//*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << dx_rand_a31 << '\t' << dx_rand_c31 << '\t' << iter_a_rand3 << '\t' << iter_c_rand3 << '\t' << it_res_a_rand3 << '\t' << it_res_c_rand3 << '\t' << eff_rand3 * 100.0 / n_tests << '\t' << eff_rand31 * 100.0 / n_tests << '\t' << res_aver_a_rand3 << '\t' << res_aver_c_rand3 << '\t' << res_max_a_rand3 << '\t' << res_max_c_rand3 << '\t' << res_diff_a_rand3 << '\t' << res_diff_c_rand3 << '\n';

			//output = &output_file_de_rand4;
			//*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << dx_rand_a41 << '\t' << dx_rand_c41 << '\t' << iter_a_rand4 << '\t' << iter_c_rand4 << '\t' << it_res_a_rand4 << '\t' << it_res_c_rand4 << '\t' << eff_rand4 * 100.0 / n_tests << '\t' << eff_rand41 * 100.0 / n_tests << '\t' << res_aver_a_rand4 << '\t' << res_aver_c_rand4 << '\t' << res_max_a_rand4 << '\t' << res_max_c_rand4 << '\t' << res_diff_a_rand4 << '\t' << res_diff_c_rand4 << '\n';

			//output = &output_file_de_rand5;
			//*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << dx_rand_a51 << '\t' << dx_rand_c51 << '\t' << iter_a_rand5 << '\t' << iter_c_rand5 << '\t' << it_res_a_rand5 << '\t' << it_res_c_rand5 << '\t' << eff_rand5 * 100.0 / n_tests << '\t' << eff_rand51 * 100.0 / n_tests << '\t' << res_aver_a_rand5 << '\t' << res_aver_c_rand5 << '\t' << res_max_a_rand5 << '\t' << res_max_c_rand5 << '\t' << res_diff_a_rand5 << '\t' << res_diff_c_rand5 << '\n';

			//Close files
			output_file_de_rand1.close(); output_file_de_rand2.close(); output_file_de_rand3.close(); output_file_de_rand4.close(); output_file_de_rand5.close();
		}
	}
}


template <typename T>
void CartAnalysis::batchTestNLSPOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
	const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Outlier test for hybrid BFGS
	//Perform batch tests
	unsigned short n_tests = 100;
	const T nu1 = 0.25, nu2 = 0.75, nu3 = 0.01, gamma1 = 0.25, gamma2 = 2.0, lambda_min = 1.0e-6, lambda_max = 1.0e6;

	const unsigned short n_items = nl_test.size();

	//Initialize random number generatorres
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_de_rand1[256];

	strcpy(output_file_text_de_rand1, "bfgsh.log");

	//New output file variables
	static std::ofstream  output_file_de_rand1;

	unsigned int n_par = 6;
	if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
		n_par = 7;
	else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
		n_par = 5;
	else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
		n_par = 8;

	//Divide P
	for (unsigned int i = 0; i < n_items; i++)
	{
		nl_test[i]->setX(nl_test[i]->getX() / 100000);
		nl_test[i]->setY(nl_test[i]->getY() / 100000);

		//nl_test[i]->setX(nl_test[i]->getX());
		//nl_test[i]->setY(nl_test[i]->getY());
	}

	//Max coordinate errors
	const T sx = 5.0 / 3 / sqrt(2.0) / 1000;
	//const T sx = 5.0 / 3 / sqrt(2.0) * 100;
	const T sy = sx;
	const T slat = 4.5 / 3.0 / sqrt(2.0);
	const T slon = slat;

	// Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 0.1;

	T mult1 = 1, mult2 = 0.5, mult3 = 0.01;
	unsigned int mode = 2;

	//Perform for all noises
	for (T k = 2.0; k <= 3; k += 0.5)
	{
		//Open output log file
		output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app);

		std::cout << " >>>>>>>>>>>>> k = " << k << " <<<<<<<<<<<<<<<<<<" << '\n';

		//Write to file
		output = &output_file_de_rand1;
		*output << "\n K = " << k << '\n';

		//Close files
		output_file_de_rand1.close();

		//Process for all noise levels
		unsigned int percs[] = { 10, 20, 30, 40, 45, 50 };
		//unsigned int percs[] = { 45, 50 };
		for (unsigned int i_per = 0; i_per < 6; i_per++)
		{

			//Open output log files
			output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app);

			//Create variables for measurements
			T  min_cost_a1 = 0, min_cost_a2 = 0, min_cost_a3 = 0, min_cost_a4 = 0, min_cost_a5 = 0,
				min_cost_a11 = 0, min_cost_a21 = 0, min_cost_a31 = 0, min_cost_a41 = 0, min_cost_a51 = 0,
				min_cost_c1 = 0, min_cost_c2 = 0, min_cost_c3 = 0, min_cost_c4 = 0, min_cost_c5 = 0,
				min_cost_c11 = 0, min_cost_c21 = 0, min_cost_c31 = 0, min_cost_c41 = 0, min_cost_c51 = 0,
				total_efficiency_rand1 = 0, total_efficiency_rand2 = 0, total_efficiency_rand3 = 0, total_efficiency_rand4 = 0, total_efficiency_rand5 = 0,
				dx_rand_a1 = 0, dx_rand_a2 = 0, dx_rand_a3 = 0, dx_rand_a4 = 0, dx_rand_a5 = 0,
				dx_rand_c1 = 0, dx_rand_c2 = 0, dx_rand_c3 = 0, dx_rand_c4 = 0, dx_rand_c5 = 0,
				dx_rand_a11 = 0, dx_rand_a21 = 0, dx_rand_a31 = 0, dx_rand_a41 = 0, dx_rand_a51 = 0,
				dx_rand_c11 = 0, dx_rand_c21 = 0, dx_rand_c31 = 0, dx_rand_c41 = 0, dx_rand_c51 = 0,
				diff1 = 0, diff2 = 0, jac_a11 = 0, jac_c11 = 0;
			unsigned int eff_1 = 0, eff_2 = 0, eff_3 = 0, eff_4 = 0, eff_5 = 0,
				eff_11 = 0, eff_21 = 0, eff_31 = 0, eff_41 = 0, eff_51 = 0,
				iter_a1 = 0, iter_a2 = 0, iter_a3 = 0, iter_a4 = 0, iter_a5 = 0,
				iter_a11 = 0, iter_a21 = 0, iter_a31 = 0, iter_a41 = 0, iter_a51 = 0,
				jac_a1 = 0, jac_a2 = 0, jac_a3 = 0, jac_a4 = 0, jac_a5 = 0,
				iter_c1 = 0, iter_c2 = 0, iter_c3 = 0, iter_c4 = 0, iter_c5 = 0,
				iter_c11 = 0, iter_c21 = 0, iter_c31 = 0, iter_c41 = 0, iter_c51 = 0,
				jac_c1 = 0, jac_c2 = 0, jac_c3 = 0, jac_c4 = 0, jac_c5 = 0, jac11 = 0,
				efficiency_rand1 = 0, efficiency_rand2 = 0, efficiency_rand3 = 0, efficiency_rand4 = 0, efficiency_rand5 = 0,
				false_out_rand1 = 0, false_out_rand2 = 0, false_out_rand3 = 0, false_out_rand4 = 0, false_out_rand5 = 0;

			const T perc = percs[i_per];
			//const T perc = 30;

			std::cout << " *********** perc = " << perc << "************* \n";

			//Repeat for all tests
			for (unsigned int i = 0; i < n_tests; i++)
			{
				Container <Node3DCartesian <T> *> nl_test_out = nl_test;
				Container <Point3DGeographic <T> *> pl_reference_out = pl_reference;

				std::cout << "Test " << i << "/" << n_tests << '\n';

				//Create random permutation
				unsigned int n_rand = (perc * n_items) / 100;
				//n_rand = 6;

				Matrix <unsigned int> I = RandomPermutation::randperm(n_items, n_rand);

				/*
				Matrix <unsigned int> I(1, 6);

				I(0, 0) = 0;
				I(0, 1) = 2;
				I(0, 2) = 4;
				I(0, 3) = 6;
				I(0, 4) = 8;
				I(0, 5) = 10;
				*/
				//Create initial weight matrix
				Matrix <T> WI = ones(2 * n_items, 2 * n_items, 1.0);

				for (unsigned int j = 0; j < n_rand; j++)
				{
					WI(I(0, j), I(0, j)) = 0.0;
					WI(I(0, j) + n_items, I(0, j) + n_items) = 0.0;
				}
				/*
				nl_test_out[I(0, 0)]->setX(nl_test[I(0, 0)]->getX() + 3 * sx);
				nl_test_out[I(0, 0)]->setY(nl_test[I(0, 0)]->getY() - 3 * sy);
				nl_test_out[I(0, 1)]->setX(nl_test[I(0, 1)]->getX() + 4 * sx);
				nl_test_out[I(0, 1)]->setY(nl_test[I(0, 1)]->getY() - 4 * sy);
				nl_test_out[I(0, 2)]->setX(nl_test[I(0, 2)]->getX() + 4 * sx);
				nl_test_out[I(0, 2)]->setY(nl_test[I(0, 2)]->getY() - 4 * sy);
				nl_test_out[I(0, 3)]->setX(nl_test[I(0, 3)]->getX() + 5 * sx);
				nl_test_out[I(0, 3)]->setY(nl_test[I(0, 3)]->getY() - 5 * sy);
				nl_test_out[I(0, 4)]->setX(nl_test[I(0, 4)]->getX() - 3 * sx);
				nl_test_out[I(0, 4)]->setY(nl_test[I(0, 4)]->getY() + 3 * sy);
				nl_test_out[I(0, 5)]->setX(nl_test[I(0, 5)]->getX() - 4 * sx);
				nl_test_out[I(0, 5)]->setY(nl_test[I(0, 5)]->getY() + 4 * sy);

				*/

				//Contaminated randomly selected points on P: outliers
				for (unsigned int j = 0; j < n_rand; j++)
				{
					//Get coordinates
					T x = nl_test[I(0, j)]->getX();
					T y = nl_test[I(0, j)]->getY();

					//Get random numbers
					const T r1 = ((T)rand() / (RAND_MAX));
					const T r2 = ((T)rand() / (RAND_MAX));

					//Get sign
					const T sig_x = copysign(1.0, r1 - 0.5);
					const T sig_y = copysign(1.0, r2 - 0.5);

					const T r3 = ((T)rand() / (RAND_MAX));
					const T r4 = ((T)rand() / (RAND_MAX));

					//Get random numbers
					const T dx = 3 + r3 * 6;
					const T dy = 3 + r4 * 6;

					//Errors
					const T deltax = sig_x * dx * sx;
					const T deltay = sig_y * dy * sy;

					//Add errors
					x = x + mult1 * deltax;
					y = y + mult1 * deltay;

					//Create new points
					Node3DCartesian <T> *p = new Node3DCartesian <T>(x, y);

					//Update coordinates
					nl_test_out[I(0, j)]->setX(x);
					nl_test_out[I(0, j)]->setY(y);
				}


				//Contaminated P: random errors, all items
				if (mode == 2)
				{
					for (unsigned int j = 0; j < n_items; j++)
					{
						//Get cooridnates
						T x1 = nl_test_out[j]->getX();
						T y1 = nl_test_out[j]->getY();
						T lat1 = pl_reference_out[j]->getLat();
						T lon1 = pl_reference_out[j]->getLon();

						//Get random numbers
						const T r1 = ((T)rand() / (RAND_MAX));
						const T r2 = ((T)rand() / (RAND_MAX));
						const T r3 = ((T)rand() / (RAND_MAX));
						const T r4 = ((T)rand() / (RAND_MAX));

						//Get signs
						const T sig_x = copysign(1.0, r1 - 0.5);
						const T sig_y = copysign(1.0, r2 - 0.5);
						const T sig_lat = copysign(1.0, r3 - 0.5);
						const T sig_lon = copysign(1.0, r4 - 0.5);

						//Get random numbers
						const T dx = ((T)rand() / (RAND_MAX));
						const T dy = ((T)rand() / (RAND_MAX));
						const T dlat = ((T)rand() / (RAND_MAX));
						const T dlon = ((T)rand() / (RAND_MAX));

						//Errors
						const T deltax = sig_x * dx * sx;
						const T deltay = sig_y * dy * sy;
						const T deltalat = sig_lat * dlat * slat; //
						const T deltalon = sig_lon * dlon * slon; //

						//slon multiplier
						const T m = 1.0 / cos(lat1 * M_PI / 180);

						//Add errors
						x1 = x1 + mult2 * deltax;
						y1 = y1 + mult2 * deltay;
						lat1 = lat1 + mult2 * deltalat;
						lon1 = lon1 + m * mult2 * deltalon;

						//Update first point
						nl_test_out[j]->setX(x1);
						nl_test_out[j]->setY(y1);

						//Update second point
						pl_reference_out[j]->setLat(lat1);
						pl_reference_out[j]->setLon(lon1);
					}
				}

				//Lowest and highest values, range
				T noise = 50;
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_dx = high_dx - low_dx;

				//Random values
				T rand_R = (low_R + range_R*rand() / (RAND_MAX + 1.0));
				T rand_latp = latp_min + (latp_max - latp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lonp = lonp_min + (lonp_max - lonp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lat0 = lat0_min + (lat0_max - lat0_min) * rand() / (RAND_MAX + 1.0);
				T rand_dx = low_dx + range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy = low_dx + range_dx*rand() / (RAND_MAX + 1.0);

				Matrix <T> X(n_par, 1), XMIN(n_par, 1), XMAX(n_par, 1), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1),
					XEQDC(n_par, 1), XLAEA(n_par, 1), XMERC(n_par, 1), XOK(n_par, 1, 1), DXT(n_par, 1);

				//Store randomly generated values
				X(0, 0) = rand_R;
				X(1, 0) = rand_latp;
				X(2, 0) = rand_lonp;
				X(3, 0) = rand_lat0;

				//Set intervals
				XMIN(0, 0) = 0; XMAX(0, 0) = low_R + range_R;
				XMIN(1, 0) = latp_min; XMAX(1, 0) = latp_max;
				XMIN(2, 0) = lonp_min; XMAX(2, 0) = lonp_max;
				XMIN(3, 0) = lat0_min; XMAX(3, 0) = lat0_max;
				XMIN(4, 0) = 0.0; XMAX(4, 0) = 0.0;

				//Correct values
				XEQDC(0, 0) = 0.06380;
				XEQDC(1, 0) = 90;
				XEQDC(2, 0) = 0;
				XEQDC(3, 0) = 50;
				XEQDC(4, 0) = 0.0;

				XLAEA(0, 0) = 0.06380;
				XLAEA(1, 0) = 52;
				XLAEA(2, 0) = 10;
				XLAEA(3, 0) = 0.0;
				XLAEA(4, 0) = 0.0;

				XMERC(0, 0) = 0.06380;
				XMERC(1, 0) = 0;
				XMERC(2, 0) = 90;
				XMERC(3, 0) = 47.0;
				XMERC(4, 0) = 0.0;

				//Threshold
				DXT(0, 0) = 6380;
				DXT(1, 0) = 180;
				DXT(2, 0) = 360;
				DXT(3, 0) = 180;

				//Set intervals
				if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
				{
					XMIN(5, 0) = 0.0; XMAX(5, 0) = 10000000;
				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
				{
					X(5, 0) = 1;
					X(6, 0) = 80;

					XMIN(5, 0) = 0.0; XMAX(5, 0) = 10000000;
					XMIN(6, 0) = -MAX_LON; XMAX(6, 0) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
				{
					X(0, 0) = rand_latp;
					X(1, 0) = rand_lonp;
					X(2, 0) = rand_lat0;
					X(3, 0) = 0;
					X(4, 0) = 1;

					XMIN(0, 0) = latp_min; XMAX(0, 0) = latp_max;
					XMIN(1, 0) = lonp_min; XMAX(1, 0) = lonp_max;
					XMIN(2, 0) = lat0_min; XMAX(2, 0) = lat0_max;
					XMIN(3, 0) = 0.0; XMAX(3, 0) = 0.0;
					XMIN(4, 0) = 0.0; XMAX(4, 0) = 10000000;

					//Correct values
					XEQDC(0, 0) = 90;
					XEQDC(1, 0) = 0;
					XEQDC(2, 0) = 50;
					XEQDC(3, 0) = 0.0;

					XLAEA(0, 0) = 52;
					XLAEA(1, 0) = 10;
					XLAEA(2, 0) = 0.0;
					XLAEA(3, 0) = 0.0;

					XMERC(0, 0) = 0;
					XMERC(1, 0) = 90;
					XMERC(2, 0) = 47.0;
					XMERC(3, 0) = 0.0;

					DXT(0, 0) = 180;
					DXT(1, 0) = 360;
					DXT(2, 0) = 180;
				}

				else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
				{
					//X(5, 0) = rand_dx;
					//X(6, 0) = rand_dy;
					X(7, 0) = 1;

					XMIN(5, 0) = -1.0e09; XMAX(5, 0) = 1.0e9;
					XMIN(6, 0) = -1.0e09; XMAX(6, 0) = 1.0e9;
					XMIN(7, 0) = 0.0; XMAX(7, 0) = 10000000;
				}


				//Update threshold
				DXT = DXT * 1.0 / 200;

				//Get projection ID
				const char * proj_text = proj->getName();
				if (strcmp(proj_text, "eqdc") == 0)
				{
					XOK = XEQDC;
				}

				else if (strcmp(proj_text, "laea") == 0)
				{
					XOK = XLAEA;
				}

				else if (strcmp(proj_text, "merc") == 0)
				{
					XOK = XMERC;
				}


				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << i << "/" << n_tests << '\n';
				std::cout << "R= " << rand_R << "   latp= " << rand_latp << "   lonp= " << rand_lonp << "   lato0= " << rand_lat0 << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 2.0e8;

				bool test_huber = true, test_andrew = false, test_danish = false, test_yang = false, test_tukey = false;
				/*
				//******************* Iterative solution
				{
				unsigned int jac = 0;
				T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;

				Matrix <T> XX = X, XX2 = X;
				W = eye(2 * n_items, 2 * n_items, 1.0);
				Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

				analysis_parameters.remove_outliers = true;

				for (unsigned int j = 0; j < 10; j++)
				{
				//q1 = 1, q2 = 1, R_def = 1, XX = X;


				//Matrix <T> WW = diag(W);
				//WW.print();

				T min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
				R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, HuberFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, 0.00001 * max_diff, output);

				if (j == 0)
				{
				XX2 = XX;
				std::cout << "XX2\n";
				XX2.print();
				}
				/*
				Point3DGeographic <T> cart_pole(XX(0, 0), XX(1, 0));
				proj->setR(R_def);
				proj->setCartPole(cart_pole);
				proj->setLat0(XX(2, 0));
				proj->setLon0(XX(3, 0));
				proj->setDx(0.0);
				proj->setDy(0.0);
				proj->setC(XX(4, 0));

				T mini = norm(trans(V) * W * V);
				std::cout << mini << '\n';
				//W = eye(2 * n_items, 2 * n_items, 1.0);

				//Compute coordinate differences (residuals): items of V matrix
				Container <Node3DCartesianProjected <T> *> nl_projected_temp;

				for (unsigned int i = 0; i < n_items; i++)
				{
				//Get type of the direction
				TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

				//Reduce lon
				const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), XX(3, 0));

				T lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0;

				try
				{
				//Convert geographic point to oblique aspect
				lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, XX(0, 0), XX(1, 0));
				lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, lat_trans, XX(0, 0), XX(1, 0), trans_lon_dir);

				//Compute x, y coordinates
				x = ArithmeticParser::parseEquation(proj->getXEquatPostfix(), lat_trans, lon_trans, R_def, proj->getA(), proj->getB(), XX(4, 0), XX(2, 0), proj->getLat1(), proj->getLat2(), false);
				y = ArithmeticParser::parseEquation(proj->getYEquatPostfix(), lat_trans, lon_trans, R_def, proj->getA(), proj->getB(), XX(4, 0), XX(2, 0), proj->getLat1(), proj->getLat2(), false);
				}

				catch (Exception &error)
				{
				//Disable point from analysis: set weight to zero
				W(i, i) = 0; W(i + n_items, i + n_items) = 0;
				}

				//Create new cartographic point
				Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

				//Add point to the list
				nl_projected_temp.push_back(n_projected);
				}

				// Weighted Helmert transformation
				Matrix <T> C(2, 2), beta(4, 1), P(n_items, 2), Q(n_items, 2), XMIN(2 * n_items, 4);

				for (unsigned int i = 0; i < n_items; i++)
				{
				P(i, 0) = nl_test[i]->getX(); P(i, 1) = nl_test[i]->getY();
				Q(i, 0) = nl_projected_temp[i]->getX(); Q(i, 1) = nl_projected_temp[i]->getY();
				}

				HelmertTransformation2D::getTransformKey2(P, Q, W, XMIN, beta, Y, C);

				//Get centers of gravity
				T x_mass_reference = C(1, 0);
				T y_mass_reference = C(1, 1);
				T x_mass_test = beta(2, 0);		//Determined from the transformation
				T y_mass_test = beta(3, 0);		//Determined from the trandformation

				//Get coefficients
				q1 = beta(0, 0);
				q2 = beta(1, 0);
				T alpha = atan2(q2, q1) * 180.0 / M_PI;

				//Remove outliers
				Matrix <T> PR(n_items, 2), QR(n_items, 2), Eps(2 * n_items, 1);
				for (unsigned int i = 0; i < n_items; i++)
				{
				PR(i, 0) = (nl_test[i]->getX() - x_mass_test);
				PR(i, 1) = (nl_test[i]->getY() - y_mass_test);

				QR(i, 0) = (q1 * (nl_projected_temp[i]->getX() - x_mass_reference) - q2 * (nl_projected_temp[i]->getY() - y_mass_reference));
				QR(i, 1) = (q2 * (nl_projected_temp[i]->getX() - x_mass_reference) + q1 * (nl_projected_temp[i]->getY() - y_mass_reference));
				}

				//Remove  outliers
				T eps_init = 0, eps = 0;
				unsigned int iter = 0;
				Outliers::findOutliersME(PR, QR, k, 1.0e-10, SimilarityScheme, DanishFunction2, 30, W, IX, Eps, eps_init, eps, iter);

				HelmertTransformation2D::getTransformKey2(P, Q, W, XMIN, beta, Y, C);

				//Get centers of gravity
				x_mass_reference = C(1, 0);
				y_mass_reference = C(1, 1);
				x_mass_test = beta(2, 0);		//Determined from the transformation
				y_mass_test = beta(3, 0);		//Determined from the trandformation

				//Get coefficients
				q1 = beta(0, 0);
				q2 = beta(1, 0);
				alpha = atan2(q2, q1) * 180.0 / M_PI;

				*/
				//Matrix <T> WWW = diag(W);
				//WWW.print();
				/*
				XX.print();
				XX2.print();
				}

				//Efficiency
				T efficiency = 0, false_out = 0, total_efficiency;
				for (unsigned int j = 0; j < 2 * n_items; j++)
				{
				if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
				if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
				}

				efficiency = 100.0 * efficiency / (2 * n_rand);
				false_out = 100.0 * false_out / (2 * n_items);
				total_efficiency = efficiency * (1 - false_out / 100.0);

				//All cases
				efficiency_rand1 += efficiency;
				false_out_rand1 += false_out;
				total_efficiency_rand1 += total_efficiency;
				std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

				min_cost_a1 += min_cost;
				iter_a1 += iterations;
				jac_a1 += jac;
				std::cout << "(" << iterations << ") ";
				std::cout << min_cost;

				// Succcess case
				if (n_par == 5) XX(4, 0) = 0;
				if (strcmp(proj_text, "eqdc") == 0)
				XX(1, 0) = 0;

				Matrix <T> DX = abs(abs(XX) - XOK);
				dx_rand_a1 += norm(DX * trans(DX));
				bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

				if (detected)
				{
				eff_1++;
				min_cost_c1 += min_cost;
				iter_c1 += iterations;
				jac_c1 += jac;
				dx_rand_c1 += norm(DX * trans(DX));
				}

				else
				{
				//XX.print();
				}

				XX.print();
				XX2.print();
				}
				//********************
				*/

				//HUBER function
				//**********************************************************************************************************************************
				if (test_huber)
				{
					//X(0, 0) = 32.4206543;
					//X(1, 0) = 27.4108887;
					//X(2, 0) = 8.1736755;
					/*
					X(0, 0) = -7;
					X(1, 0) = -103;
					X(2, 0) = 39;
					*/

					//X(0, 0) = 86.61;
					//X(1, 0) = 125.28;
					//X(2, 0) = 47.21;
					//BFGSH
					std::cout << " Huber ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;

					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);
					/*
					XX(0, 0) = -45.9613037;
					XX(1, 0) =57.9748535;
					XX(2, 0) = 21.2915039;
					XX(3, 0) = 0;
					XX(4, 0) = 1;
					*/

					//XX(0, 0) = 90;
					//XX(1, 0) = 0;
					//XX(2, 0) = 10;
					//XX.print();

					//nl_test_out.print();
					//R = 4524.8219604 latp = 32.4206543 lonp = 27.4108887 lato0 = 8.1736755
					//R = 4640.0494995 latp = -55.5084229 lonp = -149.0185547 lato0 = 60.1884460
					//X(0, 0) = 32.4206543;
					//X(1, 0) = 27.4108887;
					//X(2, 0) = 8.1736755;
					/*
					X(0, 0) = -7;
					X(1, 0) = -103;
					X(2, 0) = 39;
					*/
					/*
					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, HuberFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
					x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, HuberFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					*/
					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand1 += efficiency;
					false_out_rand1 += false_out;
					total_efficiency_rand1 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a1 += min_cost;
					iter_a1 += iterations;
					jac_a1 += jac;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(1, 0) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a1 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

					if (detected)
					{
						eff_1++;
						min_cost_c1 += min_cost;
						iter_c1 += iterations;
						jac_c1 += jac;
						dx_rand_c1 += norm(DX * trans(DX));
					}

					else
					{
						//XX.print();
					}

					XX.print();

					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;
					X2 = X;

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, HuberFunction, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 *  max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, HuberFunction, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);

					if (n_par == 5) X2(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						X2(1, 0) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					diff2 = norm(DX2 * trans(DX2));
					dx_rand_a11 += norm(DX2 * trans(DX2));
					iter_a11 += iterations;
					jac_a11 += jac;
					min_cost_a11 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(1, 0) < DXT(1, 0)) && (DX2(2, 0) < DXT(2, 0)) && (DX2(3, 0) < DXT(3, 0));
					if (detected)
					{
						eff_11++;
						min_cost_c11 += min_cost;
						iter_c11 += iterations;
						jac_c11 += jac;
						dx_rand_c11 += norm(DX2 * trans(DX2));
					}


					diff1 = norm(trans(DX) * DX);

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a1 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c1 << "\t DNWSC = " << dx_rand_c11 << "\t EFFW = " << eff_1 * 100.0 / (i + 1) << "\t EFFNW = " << eff_11 * 100.0 / (i + 1) << '\n';
					X2.print();
				}

				//ANDREW function
				//**********************************************************************************************************************************
				if (test_andrew)
				{
					//BFGSH
					std::cout << " Andrew ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, AndrewFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, AndrewFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand2 += efficiency;
					false_out_rand2 += false_out;
					total_efficiency_rand2 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a2 += min_cost;
					iter_a2 += iterations;
					jac_a2 += jac;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(1, 0) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a2 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

					if (detected)
					{
						eff_2++;
						min_cost_c2 += min_cost;
						iter_c2 += iterations;
						jac_c2 += jac;
						dx_rand_c2 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
					x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, AndrewFunction, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);

					if (n_par == 5) X2(4, 0) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(1, 0) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a21 += norm(DX2 * trans(DX2));
					iter_a21 += iterations;
					min_cost_a21 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(1, 0) < DXT(1, 0)) && (DX2(2, 0) < DXT(2, 0)) && (DX2(3, 0) < DXT(3, 0));
					if (detected)
					{
					eff_21++;
					min_cost_c21 += min_cost;
					dx_rand_c21 += norm(DX2 * trans(DX2));
					iter_c21 += iterations;
					}

					T diff1 = norm(trans(DX) * DX);
					*/
					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a2 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c2 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}

				//DANISH function
				//**********************************************************************************************************************************
				if (test_danish)
				{
					//BFGSH
					std::cout << " Danish ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, DanishFunction2, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, DanishFunction2, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand3 += efficiency;
					false_out_rand3 += false_out;
					total_efficiency_rand3 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a3 += min_cost;
					iter_a3 += iterations;
					jac_a3 += jac;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(1, 0) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a3 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

					if (detected)
					{
						eff_3++;
						min_cost_c3 += min_cost;
						iter_c3 += iterations;
						jac_c3 += jac;
						dx_rand_c3 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
					x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, DanishFunction2, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, 0.0001 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					if (n_par == 5) X2(4, 0) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(1, 0) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a31 += norm(DX2 * trans(DX2));

					iter_a31 += iterations;
					min_cost_a31 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(1, 0) < DXT(1, 0)) && (DX2(2, 0) < DXT(2, 0)) && (DX2(3, 0) < DXT(3, 0));
					if (detected)
					{
					eff_31++;
					min_cost_c31 += min_cost;
					dx_rand_c31 += norm(DX2 * trans(DX2));
					iter_c31 += iterations;
					}

					T diff1 = norm(trans(DX) * DX);
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a3 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c3 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}


				//YANG function
				//**********************************************************************************************************************************
				if (test_yang)
				{
					//X(0, 0) = -37.8808594;
					//X(1, 0) = 106.6882324;
					//X(2, 0) = 27.7148438;

					//BFGSH
					std::cout << " Yang ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, YangFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 *  max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, YangFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand4 += efficiency;
					false_out_rand4 += false_out;
					total_efficiency_rand4 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a4 += min_cost;
					iter_a4 += iterations;
					jac_a4 += jac;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(1, 0) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a4 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

					if (detected)
					{
						eff_4++;
						min_cost_c4 += min_cost;
						iter_c4 += iterations;
						jac_c4 += jac;
						dx_rand_c4 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
					x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, YangFunction, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, 0.0001 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					if (n_par == 5) X2(4, 0) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(1, 0) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));
					dx_rand_a41 += norm(DX2 * trans(DX2));
					iter_a41 += iterations;
					min_cost_a41 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(1, 0) < DXT(1, 0)) && (DX2(2, 0) < DXT(2, 0)) && (DX2(3, 0) < DXT(3, 0));
					if (detected)
					{
					eff_41++;
					min_cost_c41 += min_cost;
					dx_rand_c41 += norm(DX2 * trans(DX2));
					iter_c41 += iterations;
					}

					T diff1 = norm(trans(DX) * DX);
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a4 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c4 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}


				// TUKEY function
				//**********************************************************************************************************************************
				if (test_tukey)
				{
					//BFGSH
					std::cout << " Tukey ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, TukeyFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
						x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, TukeyFunction, k, IX, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, mult3 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
						min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
						analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, XX, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand5 += efficiency;
					false_out_rand5 += false_out;
					total_efficiency_rand5 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a5 += min_cost;
					iter_a5 += iterations;
					jac_a5 += jac;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(4, 0) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(1, 0) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a5 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(1, 0) < DXT(1, 0)) && (DX(2, 0) < DXT(2, 0)) && (DX(3, 0) < DXT(3, 0));

					if (detected)
					{
						eff_5++;
						min_cost_c5 += min_cost;
						iter_c5 += iterations;
						jac_c5 += jac;
						dx_rand_c5 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					/*
					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV2 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test_out, pl_reference_out, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions, jac), FAnalyzeProjV3 <T>(nl_test_out, pl_reference_out, nl_projected, meridians, parallels, faces_test, proj,
					x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresRot2Method)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ4 <T>(nl_test_out, pl_reference_out, proj, aspect, R_def, q1, q2, analysis_parameters.print_exceptions, jac), FAnalyzeProjV4 <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					R_def, q1, q2, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, iterations, TukeyFunction, k, IX, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, 0.0001 * max_diff, output);
					else if (analysis_parameters.analysis_method == NonLinearLeastSquaresShiftsMethod)
					min_cost = NonLinearLeastSquares::BFGSH(FAnalyzeProjJ <T>(nl_test_out, pl_reference_out, proj, aspect, analysis_parameters.print_exceptions, jac), FAnalyzeProjV <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj,
					analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W2, X2, Y, V, XMIN, XMAX, iterations, alpha_mult, nu, eps, max_iter, max_diff, output);


					if (n_par == 5) X2(4, 0) = 0;
					if (strcmp(proj_text, "eqdc") == 0)
					X2(1, 0) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					T diff2 = norm(DX2 * trans(DX2));

					dx_rand_a51 += norm(DX2 * trans(DX2));
					iter_a51 += iterations;
					min_cost_a51 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(1, 0) < DXT(1, 0)) && (DX2(2, 0) < DXT(2, 0)) && (DX2(3, 0) < DXT(3, 0));
					if (detected)
					{
					eff_51++;
					min_cost_c51 += min_cost;
					dx_rand_c51 += norm(DX2 * trans(DX2));
					iter_c51 += iterations;
					}

					T diff1 = norm(trans(DX) * DX);
					*/

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a5 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c5 << "\t DNWSC = " << dx_rand_c11 << '\n';
				}
			}

			//output = &output_file_bfgs;
			//*output << noise << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_2 << '\t' << j2 << '\t' << eff_2 * 100.0 / n_tests << '\n';

			output = &output_file_de_rand1;
			*output << perc << '\t' << efficiency_rand1 / n_tests << '\t' << false_out_rand1 / n_tests << '\t' << total_efficiency_rand1 / n_tests << '\t' << dx_rand_a1 << '\t' << dx_rand_c1 << '\t' << iter_a1 << '\t' << iter_c1 << '\t' << jac_a1 << '\t' << jac_c1 << '\t' << eff_1 * 100.0 / n_tests << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\n';
			*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << iter_a2 << '\t' << iter_c2 << '\t' << jac_a2 << '\t' << jac_c2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\n';
			*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << iter_a3 << '\t' << iter_c3 << '\t' << jac_a3 << '\t' << jac_c3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\n';
			*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << iter_a4 << '\t' << iter_c4 << '\t' << jac_a4 << '\t' << jac_c4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\n';
			*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << iter_a5 << '\t' << iter_c5 << '\t' << jac_a5 << '\t' << jac_c5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\n';
			*output << perc << '\t' << "0" << '\t' << "0" << '\t' << "0" << '\t' << dx_rand_a11 << '\t' << dx_rand_c11 << '\t' << iter_a11 << '\t' << iter_c11 << '\t' << jac_a11 << '\t' << jac_c11 << '\t' << eff_11 * 100.0 / n_tests << '\t' << min_cost_a11 << '\t' << min_cost_c11 << '\n';
			//*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << dx_rand_a21 << '\t' << dx_rand_c21 << '\t' << iter_a2 << '\t' << iter_c2 << '\t' << iter_a21 << '\t' << iter_c21 << '\t' << jac_a2 << '\t' << jac_c2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << eff_21 * 100.0 / n_tests << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << min_cost_a21 << '\t' << min_cost_c21 << '\n';
			//*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << dx_rand_a31 << '\t' << dx_rand_c31 << '\t' << iter_a3 << '\t' << iter_c3 << '\t' << iter_a31 << '\t' << iter_c31 << '\t' << jac_a3 << '\t' << jac_c3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << eff_31 * 100.0 / n_tests << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\t' << min_cost_a31 << '\t' << min_cost_c31 << '\n';
			//*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << dx_rand_a41 << '\t' << dx_rand_c41 << '\t' << iter_a4 << '\t' << iter_c4 << '\t' << iter_a41 << '\t' << iter_c41 << '\t' << jac_a4 << '\t' << jac_c4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << eff_41 * 100.0 / n_tests << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\t' << min_cost_a41 << '\t' << min_cost_c41 << '\n';
			//*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << dx_rand_a51 << '\t' << dx_rand_c51 << '\t' << iter_a5 << '\t' << iter_c5 << '\t' << iter_a51 << '\t' << iter_c51 << '\t' << jac_a5 << '\t' << jac_c5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << eff_51 * 100.0 / n_tests << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\t' << min_cost_a51 << '\t' << min_cost_c51 << '\n';


			//Close file
			output_file_de_rand1.close();
		}
	}
}



template <typename T>
void CartAnalysis::batchTestNelderMeadOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
	const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, const T max_diff, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	//Outlier test for hybrid Nelder Mead method
	//Perform batch tests
	unsigned short n_tests = 100;

	const unsigned short n_items = nl_test.size();

	//Initialize random number generatorres
	srand((unsigned)time(0));

	//Create file names
	char output_file_text_de_rand1[256];

	strcpy(output_file_text_de_rand1, "nelder.log");

	//New output file variables
	static std::ofstream  output_file_de_rand1;

	unsigned int n_par = 6;
	if (analysis_parameters.analysis_method == SimplexRotMethod)
		n_par = 7;
	else if (analysis_parameters.analysis_method == SimplexRot2Method)
		n_par = 5;
	else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
		n_par = 8;

	//Divide P
	for (unsigned int i = 0; i < n_items; i++)
	{
		nl_test[i]->setX(nl_test[i]->getX() / 100000);
		nl_test[i]->setY(nl_test[i]->getY() / 100000);
	}

	//Max coordinate errors
	const T sx = 5.0 / 3 / sqrt(2.0) / 1000;
	const T sy = sx;
	const T slat = 4.5 / 3.0 / sqrt(2.0);
	const T slon = slat;

	// Get intervals
	const T latp_min = proj->getLatPInterval().min_val;
	const T latp_max = proj->getLatPInterval().max_val;
	const T lonp_min = proj->getLonPInterval().min_val;
	const T lonp_max = proj->getLonPInterval().max_val;
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;
	T R = 0.1;

	std::cout << latp_min << "  " << latp_max;

	T mult1 = 1, mult2 = 0.5, mult3 = 0.0001;
	unsigned int mode = 2;

	//Perform for all noises
	for (T k = 2.0; k <= 3; k += 0.5)
	{
		//Open output log file
		output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app);

		std::cout << " >>>>>>>>>>>>> k = " << k << " <<<<<<<<<<<<<<<<<<" << '\n';

		//Write to file
		output = &output_file_de_rand1;
		*output << "\n K = " << k << '\n';

		//Close files
		output_file_de_rand1.close();

		//Process for all noise levels
		unsigned int percs[] = { 10, 20, 30, 40, 45, 50 };
		for (unsigned int i_per = 0; i_per < 6; i_per++)
		{

			//Open output log files
			output_file_de_rand1.open(output_file_text_de_rand1, std::ofstream::app);

			//Create variables for measurements
			T  min_cost_a1 = 0, min_cost_a2 = 0, min_cost_a3 = 0, min_cost_a4 = 0, min_cost_a5 = 0,
				min_cost_a11 = 0, min_cost_a21 = 0, min_cost_a31 = 0, min_cost_a41 = 0, min_cost_a51 = 0,
				min_cost_c1 = 0, min_cost_c2 = 0, min_cost_c3 = 0, min_cost_c4 = 0, min_cost_c5 = 0,
				min_cost_c11 = 0, min_cost_c21 = 0, min_cost_c31 = 0, min_cost_c41 = 0, min_cost_c51 = 0,
				total_efficiency_rand1 = 0, total_efficiency_rand2 = 0, total_efficiency_rand3 = 0, total_efficiency_rand4 = 0, total_efficiency_rand5 = 0,
				dx_rand_a1 = 0, dx_rand_a2 = 0, dx_rand_a3 = 0, dx_rand_a4 = 0, dx_rand_a5 = 0,
				dx_rand_c1 = 0, dx_rand_c2 = 0, dx_rand_c3 = 0, dx_rand_c4 = 0, dx_rand_c5 = 0,
				dx_rand_a11 = 0, dx_rand_a21 = 0, dx_rand_a31 = 0, dx_rand_a41 = 0, dx_rand_a51 = 0,
				dx_rand_c11 = 0, dx_rand_c21 = 0, dx_rand_c31 = 0, dx_rand_c41 = 0, dx_rand_c51 = 0,
				diff1 = 0, diff2 = 0, jac_a11 = 0, jac_c11 = 0;
			unsigned int eff_1 = 0, eff_2 = 0, eff_3 = 0, eff_4 = 0, eff_5 = 0,
				eff_11 = 0, eff_21 = 0, eff_31 = 0, eff_41 = 0, eff_51 = 0,
				iter_a1 = 0, iter_a2 = 0, iter_a3 = 0, iter_a4 = 0, iter_a5 = 0,
				iter_a11 = 0, iter_a21 = 0, iter_a31 = 0, iter_a41 = 0, iter_a51 = 0,
				jac_a1 = 0, jac_a2 = 0, jac_a3 = 0, jac_a4 = 0, jac_a5 = 0,
				iter_c1 = 0, iter_c2 = 0, iter_c3 = 0, iter_c4 = 0, iter_c5 = 0,
				iter_c11 = 0, iter_c21 = 0, iter_c31 = 0, iter_c41 = 0, iter_c51 = 0,
				jac_c1 = 0, jac_c2 = 0, jac_c3 = 0, jac_c4 = 0, jac_c5 = 0, jac11 = 0,
				efficiency_rand1 = 0, efficiency_rand2 = 0, efficiency_rand3 = 0, efficiency_rand4 = 0, efficiency_rand5 = 0,
				false_out_rand1 = 0, false_out_rand2 = 0, false_out_rand3 = 0, false_out_rand4 = 0, false_out_rand5 = 0;

			const T perc = percs[i_per];
			//const T perc = 30;

			std::cout << " *********** perc = " << perc << "************* \n";

			//Repeat for all tests
			for (unsigned int i = 0; i < n_tests; i++)
			{
				Container <Node3DCartesian <T> *> nl_test_out = nl_test;
				Container <Point3DGeographic <T> *> pl_reference_out = pl_reference;

				std::cout << "Test " << i << "/" << n_tests << '\n';

				//Create random permutation
				unsigned int n_rand = (perc * n_items) / 100;
				//n_rand = 6;

				Matrix <unsigned int> I = RandomPermutation::randperm(n_items, n_rand);

				/*
				Matrix <unsigned int> I(1, 6);

				I(0, 0) = 0;
				I(0, 1) = 2;
				I(0, 2) = 4;
				I(0, 3) = 6;
				I(0, 4) = 8;
				I(0, 5) = 10;
				*/
				//Create initial weight matrix
				Matrix <T> WI = ones(2 * n_items, 2 * n_items, 1.0);

				for (unsigned int j = 0; j < n_rand; j++)
				{
					WI(I(0, j), I(0, j)) = 0.0;
					WI(I(0, j) + n_items, I(0, j) + n_items) = 0.0;
				}
				/*
				nl_test_out[I(0, 0)]->setX(nl_test[I(0, 0)]->getX() + 3 * sx);
				nl_test_out[I(0, 0)]->setY(nl_test[I(0, 0)]->getY() - 3 * sy);
				nl_test_out[I(0, 1)]->setX(nl_test[I(0, 1)]->getX() + 4 * sx);
				nl_test_out[I(0, 1)]->setY(nl_test[I(0, 1)]->getY() - 4 * sy);
				nl_test_out[I(0, 2)]->setX(nl_test[I(0, 2)]->getX() + 4 * sx);
				nl_test_out[I(0, 2)]->setY(nl_test[I(0, 2)]->getY() - 4 * sy);
				nl_test_out[I(0, 3)]->setX(nl_test[I(0, 3)]->getX() + 5 * sx);
				nl_test_out[I(0, 3)]->setY(nl_test[I(0, 3)]->getY() - 5 * sy);
				nl_test_out[I(0, 4)]->setX(nl_test[I(0, 4)]->getX() - 3 * sx);
				nl_test_out[I(0, 4)]->setY(nl_test[I(0, 4)]->getY() + 3 * sy);
				nl_test_out[I(0, 5)]->setX(nl_test[I(0, 5)]->getX() - 4 * sx);
				nl_test_out[I(0, 5)]->setY(nl_test[I(0, 5)]->getY() + 4 * sy);

				*/

				//Contaminated randomly selected points on P: outliers
				for (unsigned int j = 0; j < n_rand; j++)
				{
					//Get coordinates
					T x = nl_test[I(0, j)]->getX();
					T y = nl_test[I(0, j)]->getY();

					//Get random numbers
					const T r1 = ((T)rand() / (RAND_MAX));
					const T r2 = ((T)rand() / (RAND_MAX));

					//Get sign
					const T sig_x = copysign(1.0, r1 - 0.5);
					const T sig_y = copysign(1.0, r2 - 0.5);

					const T r3 = ((T)rand() / (RAND_MAX));
					const T r4 = ((T)rand() / (RAND_MAX));

					//Get random numbers
					const T dx = 3 + r3 * 6;
					const T dy = 3 + r4 * 6;

					//Errors
					const T deltax = sig_x * dx * sx;
					const T deltay = sig_y * dy * sy;

					//Add errors
					x = x + mult1 * deltax;
					y = y + mult1 * deltay;

					//Create new points
					Node3DCartesian <T> *p = new Node3DCartesian <T>(x, y);

					//Update coordinates
					nl_test_out[I(0, j)]->setX(x);
					nl_test_out[I(0, j)]->setY(y);
				}


				//Contaminated P: random errors, all items
				if (mode == 2)
				{
					for (unsigned int j = 0; j < n_items; j++)
					{
						//Get cooridnates
						T x1 = nl_test_out[j]->getX();
						T y1 = nl_test_out[j]->getY();
						T lat1 = pl_reference_out[j]->getLat();
						T lon1 = pl_reference_out[j]->getLon();

						//Get random numbers
						const T r1 = ((T)rand() / (RAND_MAX));
						const T r2 = ((T)rand() / (RAND_MAX));
						const T r3 = ((T)rand() / (RAND_MAX));
						const T r4 = ((T)rand() / (RAND_MAX));

						//Get signs
						const T sig_x = copysign(1.0, r1 - 0.5);
						const T sig_y = copysign(1.0, r2 - 0.5);
						const T sig_lat = copysign(1.0, r3 - 0.5);
						const T sig_lon = copysign(1.0, r4 - 0.5);

						//Get random numbers
						const T dx = ((T)rand() / (RAND_MAX));
						const T dy = ((T)rand() / (RAND_MAX));
						const T dlat = ((T)rand() / (RAND_MAX));
						const T dlon = ((T)rand() / (RAND_MAX));

						//Errors
						const T deltax = sig_x * dx * sx;
						const T deltay = sig_y * dy * sy;
						const T deltalat = sig_lat * dlat * slat; //
						const T deltalon = sig_lon * dlon * slon; //

						//slon multiplier
						const T m = 1.0 / cos(lat1 * M_PI / 180);

						//Add errors
						x1 = x1 + mult2 * deltax;
						y1 = y1 + mult2 * deltay;
						lat1 = lat1 + mult2 * deltalat;
						lon1 = lon1 + m * mult2 * deltalon;

						//Update first point
						nl_test_out[j]->setX(x1);
						nl_test_out[j]->setY(y1);

						//Update second point
						pl_reference_out[j]->setLat(lat1);
						pl_reference_out[j]->setLon(lon1);
					}
				}

				//Lowest and highest values, range
				T noise = 50;
				T low_R = (1 - noise / 100.0) * 6378, high_R = (1 + noise / 100.00) * 6378;
				T low_dx = -(1 + noise / 100.0) * 10000000, high_dx = (1 + noise / 100.00) * 1000000;
				T range_R = high_R - low_R;
				T range_dx = high_dx - low_dx;

				//Random values
				T rand_R = (low_R + range_R*rand() / (RAND_MAX + 1.0));
				T rand_latp = latp_min + (latp_max - latp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lonp = lonp_min + (lonp_max - lonp_min) * rand() / (RAND_MAX + 1.0);
				T rand_lat0 = lat0_min + (lat0_max - lat0_min) * rand() / (RAND_MAX + 1.0);
				T rand_dx = low_dx + range_dx*rand() / (RAND_MAX + 1.0);
				T rand_dy = low_dx + range_dx*rand() / (RAND_MAX + 1.0);

				Matrix <T> X(1, n_par), XMIN(1, n_par), XMAX(1, n_par), Y(2 * n_items, 1), W(2 * n_items, 2 * n_items, 0.0, 1), V(2 * n_items, 1),
					XEQDC(1, n_par), XLAEA(1, n_par), XMERC(1, n_par), XOK(1, n_par, 1), DXT(1, n_par);

				//Store randomly generated values
				X(0, 0) = rand_R;
				X(0, 1) = rand_latp;
				X(0, 2) = rand_lonp;
				X(0, 3) = rand_lat0;

				//Set intervals
				XMIN(0, 0) = 0; XMAX(0, 0) = low_R + range_R;
				XMIN(0, 1) = latp_min; XMAX(0, 1) = latp_max + 10; // +10 for eqdc
				XMIN(0, 2) = lonp_min; XMAX(0, 2) = lonp_max;
				XMIN(0, 3) = lat0_min; XMAX(0, 3) = lat0_max;
				XMIN(0, 4) = 0.0; XMAX(0, 4) = 0.0;

				//Correct values
				XEQDC(0, 0) = 0.06380;
				XEQDC(0, 1) = 90;
				XEQDC(0, 2) = 0;
				XEQDC(0, 3) = 50;
				XEQDC(0, 4) = 0.0;

				XLAEA(0, 0) = 0.06380;
				XLAEA(0, 1) = 52;
				XLAEA(0, 2) = 10;
				XLAEA(0, 3) = 0.0;
				XLAEA(0, 4) = 0.0;

				XMERC(0, 0) = 0.06380;
				XMERC(0, 1) = 0;
				XMERC(0, 2) = 90;
				XMERC(0, 3) = 47.0;
				XMERC(0, 4) = 0.0;

				//Threshold
				DXT(0, 0) = 6380;
				DXT(0, 1) = 180;
				DXT(0, 2) = 360;
				DXT(0, 3) = 180;

				//Set intervals
				if (analysis_parameters.analysis_method == SimplexMethod)
				{
					XMIN(0, 5) = 0.0; XMAX(0, 5) = 0;
				}

				else if (analysis_parameters.analysis_method == SimplexRotMethod)
				{
					X(0, 5) = 1;
					X(0, 6) = 80;

					XMIN(0, 5) = 0.0; XMAX(0, 5) = 10000000;
					XMIN(0, 6) = -MAX_LON; XMAX(0, 6) = MAX_LON;

				}

				else if (analysis_parameters.analysis_method == SimplexRot2Method)
				{
					X(0, 0) = rand_latp;
					X(0, 1) = rand_lonp;
					X(0, 2) = rand_lat0;
					X(0, 3) = 0;
					X(0, 4) = 1;

					XMIN(0, 0) = latp_min; XMAX(0, 0) = latp_max;
					XMIN(0, 1) = lonp_min; XMAX(0, 1) = lonp_max;
					XMIN(0, 2) = lat0_min; XMAX(0, 2) = lat0_max;
					XMIN(0, 3) = 0.0; XMAX(0, 3) = 0.0;
					XMIN(0, 4) = 0.0; XMAX(0, 4) = 0.0;

					//Correct values
					XEQDC(0, 0) = 90;
					XEQDC(0, 1) = 0;
					XEQDC(0, 2) = 50;
					XEQDC(0, 3) = 0.0;

					XLAEA(0, 0) = 52;
					XLAEA(0, 1) = 10;
					XLAEA(0, 2) = 0.0;
					XLAEA(0, 3) = 0.0;

					XMERC(0, 0) = 0;
					XMERC(0, 1) = 90;
					XMERC(0, 2) = 47.0;
					XMERC(0, 3) = 0.0;

					DXT(0, 0) = 180;
					DXT(0, 1) = 360;
					DXT(0, 2) = 180;
				}

				else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
				{
					//X(5, 0) = rand_dx;
					//X(6, 0) = rand_dy;
					X(0, 7) = 1;

					XMIN(0, 5) = -1.0e09; XMAX(0, 5) = 1.0e9;
					XMIN(0, 6) = -1.0e09; XMAX(0, 6) = 1.0e9;
					XMIN(0, 7) = 0.0; XMAX(0, 7) = 10000000;
				}


				//Update threshold
				DXT = DXT * 1.0 / 200;

				//Get projection ID
				const char * proj_text = proj->getName();
				if (strcmp(proj_text, "eqdc") == 0)
				{
					XOK = XEQDC;
				}

				else if (strcmp(proj_text, "laea") == 0)
				{
					XOK = XLAEA;
				}

				else if (strcmp(proj_text, "merc") == 0)
				{
					XOK = XMERC;
				}


				//Print actual values
				std::cout << "Noise = " << noise << '\t' << "Test " << i << "/" << n_tests << '\n';
				std::cout << "R= " << rand_R << "   latp= " << rand_latp << "   lonp= " << rand_lonp << "   lato0= " << rand_lat0 << '\n';

				//1.0
				//7e6
				//1.0e2
				//2.0e8
				//1.0e7;
				const T res2 = 2.0e8;


				unsigned int it_res = 0;

				bool test_huber = true, test_andrew = true, test_danish = true, test_yang = true, test_tukey = true;


				//HUBER function
				//**********************************************************************************************************************************
				if (test_huber)
				{
					//BFGSH
					std::cout << " Huber ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);
					/*
					XX(0, 0) = -45.9613037;
					XX(1, 0) =57.9748535;
					XX(2, 0) = 21.2915039;
					XX(3, 0) = 0;
					XX(4, 0) = 1;
					*/

					//XX(0, 0) = 90;
					//XX(1, 0) = 0;
					//XX(2, 0) = 10;
					//XX.print();

					//nl_test_out.print();

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, 3 * max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand1 += efficiency;
					false_out_rand1 += false_out;
					total_efficiency_rand1 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a1 += min_cost;
					iter_a1 += iterations;
					jac_a1 += it_res;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					std::cout << "Differ = \n";

					dx_rand_a1 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_1++;
						min_cost_c1 += min_cost;
						iter_c1 += iterations;
						jac_c1 += it_res;
						dx_rand_c1 += norm(DX * trans(DX));
					}

					else
					{
						//XX.print();
					}

					XX.print();

					//************* Non-weighted version ********************
					Matrix <T> W2 = eye(2 * n_items, 2 * n_items, 1.0), X2 = X;
					q1 = 1, q2 = 1, R_def = 1; jac = 0;
					nl_projected.clear();

					analysis_parameters.remove_outliers = false;

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, W2, X2, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W2, X2, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, HuberFunction, k, IX, output), XMIN, XMAX, W2, X2, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W2, X2, Y, V, iterations, eps, max_iter, output);

					if (n_par == 5) X2(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						X2(0, 1) = 0;

					Matrix <T> DX2 = abs(abs(X2) - XOK);
					diff2 = norm(DX2 * trans(DX2));
					dx_rand_a11 += norm(DX2 * trans(DX2));
					iter_a11 += iterations;
					jac_a11 += it_res;
					min_cost_a11 += min_cost;

					detected = (DX2(0, 0) < DXT(0, 0)) && (DX2(0, 1) < DXT(0, 1)) && (DX2(0, 2) < DXT(0, 2)) && (DX2(0, 3) < DXT(0, 3));
					if (detected)
					{
						eff_11++;
						min_cost_c11 += min_cost;
						iter_c11 += iterations;
						jac_c11 += it_res;
						dx_rand_c11 += norm(DX2 * trans(DX2));
					}


					diff1 = norm(trans(DX) * DX);

					X2.print();

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a1 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c1 << "\t DNWSC = " << dx_rand_c11 << "\t EFFW = " << eff_1 * 100.0 / (i + 1) << "\t EFFNW = " << eff_11 * 100.0 / (i + 1) << '\n';

				}

				//ANDREW function
				//**********************************************************************************************************************************
				if (test_andrew)
				{
					//BFGSH
					std::cout << " Andrew ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, AndrewFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand2 += efficiency;
					false_out_rand2 += false_out;
					total_efficiency_rand2 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a2 += min_cost;
					iter_a2 += iterations;
					jac_a2 += it_res;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a2 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_2++;
						min_cost_c2 += min_cost;
						iter_c2 += iterations;
						jac_c2 += it_res;
						dx_rand_c2 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a2 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c2 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}

				//DANISH function
				//**********************************************************************************************************************************
				if (test_danish)
				{
					//BFGSH
					std::cout << " Danish ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, DanishFunction2, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand3 += efficiency;
					false_out_rand3 += false_out;
					total_efficiency_rand3 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a3 += min_cost;
					iter_a3 += iterations;
					jac_a3 += it_res;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a3 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_3++;
						min_cost_c3 += min_cost;
						iter_c3 += iterations;
						jac_c3 += it_res;
						dx_rand_c3 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);


					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a3 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c3 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}


				//YANG function
				//**********************************************************************************************************************************
				if (test_yang)
				{
					//BFGSH
					std::cout << " Yang ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, YangFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, YangFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);

					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand4 += efficiency;
					false_out_rand4 += false_out;
					total_efficiency_rand4 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a4 += min_cost;
					iter_a4 += iterations;
					jac_a4 += it_res;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a4 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_4++;
						min_cost_c4 += min_cost;
						iter_c4 += iterations;
						jac_c4 += it_res;
						dx_rand_c4 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);


					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a4 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c4 << "\t DNWSC = " << dx_rand_c11 << '\n';

				}


				// TUKEY function
				//**********************************************************************************************************************************
				if (test_tukey)
				{
					//BFGSH
					std::cout << " Tukey ";
					T min_cost = 0, q1 = 1, q2 = 1, R_def = 1;
					unsigned int jac = 0;
					output = &output_file_de_rand1;
					nl_projected.clear();
					analysis_parameters.remove_outliers = true;

					Matrix <T> XX = X;
					W = eye(2 * n_items, 2 * n_items, 1.0);
					Matrix <unsigned int> IX = ones(2 * n_items, 1, 1.0);

					if (analysis_parameters.analysis_method == SimplexMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV2S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRotMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV3S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, XMIN, XMAX, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexRot2Method)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjV4S <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, R_def, q1, q2, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, TukeyFunction, k, IX, output), XMIN, XMAX, W, XX, Y, V, iterations, mult3 * eps, max_iter, output);
					else if (analysis_parameters.analysis_method == SimplexShiftsMethod)
						min_cost = SimplexMethod::NelderMead(FAnalyzeProjVS <T>(nl_test_out, pl_reference_out, meridians, parallels, faces_test, proj, analysis_parameters, ObliqueAspect, best_sample,
						total_created_and_analyzed_samples_projection, it_res, output), XMIN, XMAX, W, XX, Y, V, iterations, eps, max_iter, output);


					//Efficiency
					T efficiency = 0, false_out = 0, total_efficiency;
					for (unsigned int j = 0; j < 2 * n_items; j++)
					{
						if ((WI(j, j) == 0) && (IX(j, 0) < 0.001)) efficiency++;
						if ((WI(j, j) == 1) && (IX(j, 0) < 0.001)) false_out++;
					}

					efficiency = 100.0 * efficiency / (2 * n_rand);
					false_out = 100.0 * false_out / (2 * n_items);
					total_efficiency = efficiency * (1 - false_out / 100.0);

					//All cases
					efficiency_rand5 += efficiency;
					false_out_rand5 += false_out;
					total_efficiency_rand5 += total_efficiency;
					std::cout << "Efficiency = " << efficiency << "   False outliers = " << false_out << "    Total efficiency = " << total_efficiency << '\n';

					min_cost_a5 += min_cost;
					iter_a5 += iterations;
					jac_a5 += it_res;
					std::cout << "(" << iterations << ") ";
					std::cout << min_cost;

					// Succcess case
					if (n_par == 5) XX(0, 4) = 0;
					if ((strcmp(proj_text, "eqdc") == 0) && (n_par == 5))
						XX(0, 1) = 0;

					Matrix <T> DX = abs(abs(XX) - XOK);
					dx_rand_a5 += norm(DX * trans(DX));
					bool detected = (DX(0, 0) < DXT(0, 0)) && (DX(0, 1) < DXT(0, 1)) && (DX(0, 2) < DXT(0, 2)) && (DX(0, 3) < DXT(0, 3));

					if (detected)
					{
						eff_5++;
						min_cost_c5 += min_cost;
						iter_c5 += iterations;
						jac_c5 += it_res;
						dx_rand_c5 += norm(DX * trans(DX));
					}

					else
					{
						XX.print();
					}

					diff1 = norm(trans(DX) * DX);

					std::cout << "DW = " << diff1 << "\t DNW = " << diff2 << "\t DWS = " << dx_rand_a5 << "\t DNWS = " << dx_rand_a11 << "\t DWSC = " << dx_rand_c5 << "\t DNWSC = " << dx_rand_c11 << '\n';
				}
			}

			//output = &output_file_bfgs;
			//*output << noise << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << iter_2 << '\t' << j2 << '\t' << eff_2 * 100.0 / n_tests << '\n';

			output = &output_file_de_rand1;
			*output << perc << '\t' << efficiency_rand1 / n_tests << '\t' << false_out_rand1 / n_tests << '\t' << total_efficiency_rand1 / n_tests << '\t' << dx_rand_a1 << '\t' << dx_rand_c1 << '\t' << iter_a1 << '\t' << iter_c1 << '\t' << jac_a1 << '\t' << jac_c1 << '\t' << eff_1 * 100.0 / n_tests << '\t' << min_cost_a1 << '\t' << min_cost_c1 << '\n';
			*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << iter_a2 << '\t' << iter_c2 << '\t' << jac_a2 << '\t' << jac_c2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\n';
			*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << iter_a3 << '\t' << iter_c3 << '\t' << jac_a3 << '\t' << jac_c3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\n';
			*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << iter_a4 << '\t' << iter_c4 << '\t' << jac_a4 << '\t' << jac_c4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\n';
			*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << iter_a5 << '\t' << iter_c5 << '\t' << jac_a5 << '\t' << jac_c5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\n';
			*output << perc << '\t' << "0" << '\t' << "0" << '\t' << "0" << '\t' << dx_rand_a11 << '\t' << dx_rand_c11 << '\t' << iter_a11 << '\t' << iter_c11 << '\t' << jac_a11 << '\t' << jac_c11 << '\t' << eff_11 * 100.0 / n_tests << '\t' << min_cost_a11 << '\t' << min_cost_c11 << '\n';
			//*output << perc << '\t' << efficiency_rand2 / n_tests << '\t' << false_out_rand2 / n_tests << '\t' << total_efficiency_rand2 / n_tests << '\t' << dx_rand_a2 << '\t' << dx_rand_c2 << '\t' << dx_rand_a21 << '\t' << dx_rand_c21 << '\t' << iter_a2 << '\t' << iter_c2 << '\t' << iter_a21 << '\t' << iter_c21 << '\t' << jac_a2 << '\t' << jac_c2 << '\t' << eff_2 * 100.0 / n_tests << '\t' << eff_21 * 100.0 / n_tests << '\t' << min_cost_a2 << '\t' << min_cost_c2 << '\t' << min_cost_a21 << '\t' << min_cost_c21 << '\n';
			//*output << perc << '\t' << efficiency_rand3 / n_tests << '\t' << false_out_rand3 / n_tests << '\t' << total_efficiency_rand3 / n_tests << '\t' << dx_rand_a3 << '\t' << dx_rand_c3 << '\t' << dx_rand_a31 << '\t' << dx_rand_c31 << '\t' << iter_a3 << '\t' << iter_c3 << '\t' << iter_a31 << '\t' << iter_c31 << '\t' << jac_a3 << '\t' << jac_c3 << '\t' << eff_3 * 100.0 / n_tests << '\t' << eff_31 * 100.0 / n_tests << '\t' << min_cost_a3 << '\t' << min_cost_c3 << '\t' << min_cost_a31 << '\t' << min_cost_c31 << '\n';
			//*output << perc << '\t' << efficiency_rand4 / n_tests << '\t' << false_out_rand4 / n_tests << '\t' << total_efficiency_rand4 / n_tests << '\t' << dx_rand_a4 << '\t' << dx_rand_c4 << '\t' << dx_rand_a41 << '\t' << dx_rand_c41 << '\t' << iter_a4 << '\t' << iter_c4 << '\t' << iter_a41 << '\t' << iter_c41 << '\t' << jac_a4 << '\t' << jac_c4 << '\t' << eff_4 * 100.0 / n_tests << '\t' << eff_41 * 100.0 / n_tests << '\t' << min_cost_a4 << '\t' << min_cost_c4 << '\t' << min_cost_a41 << '\t' << min_cost_c41 << '\n';
			//*output << perc << '\t' << efficiency_rand5 / n_tests << '\t' << false_out_rand5 / n_tests << '\t' << total_efficiency_rand5 / n_tests << '\t' << dx_rand_a5 << '\t' << dx_rand_c5 << '\t' << dx_rand_a51 << '\t' << dx_rand_c51 << '\t' << iter_a5 << '\t' << iter_c5 << '\t' << iter_a51 << '\t' << iter_c51 << '\t' << jac_a5 << '\t' << jac_c5 << '\t' << eff_5 * 100.0 / n_tests << '\t' << eff_51 * 100.0 / n_tests << '\t' << min_cost_a5 << '\t' << min_cost_c5 << '\t' << min_cost_a51 << '\t' << min_cost_c51 << '\n';


			//Close file
			output_file_de_rand1.close();
		}
	}
}


template <typename T>
void CartAnalysis::batchTestAmountOfPoints(Container <Sample <T> > &sl, Container <Projection <T> *> &proj_list, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
	unsigned int & total_created_or_thrown_samples, std::ostream * output)
{
	//Batch test: dependencance on the amount of analyzed features, n = 5, 10, 15, 20, 25
	const unsigned int n = nl_test.size();

	Container <Node3DCartesian <T> *> nl_test_temp(nl_test), nl_test_temp2;
	Container <Point3DGeographic <T> *> nl_reference_temp(pl_reference), nl_reference_temp2;

	//Get nearest value
	unsigned int n_points = n, iterations = 0;

	while (n_points > 5)
	{
		//Find min_max box
		const Node3DCartesian <T> *p_x_max = *std::max_element(nl_test_temp.begin(), nl_test_temp.end(), sortPointsByX());
		const Node3DCartesian <T> *p_x_min = *std::min_element(nl_test_temp.begin(), nl_test_temp.end(), sortPointsByX());
		const Node3DCartesian <T> *p_y_max = *std::max_element(nl_test_temp.begin(), nl_test_temp.end(), sortPointsByY());
		const Node3DCartesian <T> *p_y_min = *std::min_element(nl_test_temp.begin(), nl_test_temp.end(), sortPointsByY());

		//Create min-max box
		const Node3DCartesian <T> *m1 = new Node3DCartesian <T>(p_x_min->getX(), p_y_min->getY());
		const Node3DCartesian <T> *m2 = new Node3DCartesian <T>(p_x_max->getX(), p_y_min->getY());
		const Node3DCartesian <T> *m3 = new Node3DCartesian <T>(p_x_max->getX(), p_y_max->getY());
		const Node3DCartesian <T> *m4 = new Node3DCartesian <T>(p_x_min->getX(), p_y_max->getY());
		
		//Create lists of closest points
		std::map<double, unsigned int> v1, v2, v3, v4;
		std::map<double, unsigned int> ::iterator iv1, iv2, iv3, iv4;

		for (unsigned int j = 0; j < n_points; j++ )
		{
			//Get point
			const Node3DCartesian <T> *node = nl_test_temp[j];

			//Measure its distance to min max box points
			const T d1 = EuclDistance::getEuclDistance(node, m1);
			const T d2 = EuclDistance::getEuclDistance(node, m2);
			const T d3 = EuclDistance::getEuclDistance(node, m3);
			const T d4 = EuclDistance::getEuclDistance(node, m4);

			//Find closest vertex for a point
			//Results stored in minimum1 and index1
			T minimum1 = d1, minimum2 = d3;
			std::map<double, unsigned int> *pv1 = &v1, *pv2 = &v3;
			
			//Closest to the second vertex
			if (d2 < d1)
			{
				minimum1 = d2; pv1 = &v2;
			}

			//Closest to the fourth vertex
			if (d4 < d3)
			{
				minimum2 = d4; pv2 = &v4;
			}

			//Compare both
			if (minimum2 < minimum1)
			{
				minimum1 = minimum2;
				pv1 = pv2;
			}

			//Add to the list
			(*pv1)[minimum1] = j;
		}
		
		//Set initial n = [30, 25, 20, 15, 10...]
		if (iterations == 0)
		{
			for (unsigned int i = 30; i > 5; i -= 5)
			{
				if (n / i >= 1)
				{
					n_points = i;
					break;
				}
			}
		}

		//Decrement amount of analyzed features
		else n_points -= 5;

		std::cout << "\n\n >>> Test, analyzed " << n_points << " points: \n\n";
		*output << "\n\n >>> Test, analyzed " << n_points << " points: \n\n";

		//Amount of points per vertex
		const unsigned int npt = n_points / 4;

		//Size of maps
		const unsigned int nv1 = v1.size(), nv2 = v2.size(), nv3 = v3.size(), nv4 = v4.size();

		//Correct amount of points closest to each mbr vertex 
		unsigned int np1 = 0, np2 = 0, np3 = 0, np4 = 0;

		//How many avaiable vertices are there
		unsigned int n1 = std::min(nv1, npt), n2 = std::min(nv2, npt), n3 = std::min(nv3, npt), n4 = std::min(nv4, npt);
		const unsigned int rem = n_points % 4 + n_points - n1 - n2 - n3 - n4;
		
		if (rem != 0)
			n1 += rem;
		
		//Read from the first map
		for (iv1 = v1.begin(); iv1 != v1.end() && np1 < n1; iv1++) 
		{
			//Get index of the point
			const unsigned int index = iv1->second;

			//Add copy of the point to the list
			nl_test_temp2.push_back(nl_test_temp[index]->clone());
			nl_reference_temp2.push_back(nl_reference_temp[index]->clone());

			np1++;
		}

		//First map is shorter
		unsigned int diff1 = (np1 < n1 ? n1 - np1 : 0);
		n2 += diff1;

		//Read from the second map
		for (iv2 = v2.begin(); iv2 != v2.end() && np2 < n2; iv2++)
		{
			//Get index of the point
			const unsigned int index = iv2->second;

			//Add copy of the point to the list
			nl_test_temp2.push_back(nl_test_temp[index]->clone());
			nl_reference_temp2.push_back(nl_reference_temp[index]->clone());

			np2++;
		}

		//Second map is shorter
		unsigned int diff2 = (np2 < n2 ? n2 - np2 : 0);
		n3 += diff2;

		//Read from the third map
		for (iv3 = v3.begin(); iv3 != v3.end() && np3 < n3; iv3++)
		{
			//Get index of the point
			const unsigned int index = iv3->second;
			
			//Add copy of the point to the list
			nl_test_temp2.push_back(nl_test_temp[index]->clone());
			nl_reference_temp2.push_back(nl_reference_temp[index]->clone());

			np3++;
		}

		//nl_reference_temp2.print(output);

		//Third map is shorter
		unsigned int diff3 = (np3 < n3 ? n3 - np3 : 0);
		n4 += diff3;

		//Read from the fourth map
		for (iv4 = v4.begin(); iv4 != v4.end() && np4 < n4; iv4++)
		{
			//Get index of the point
			const int index = iv4->second;

			//Add copy of the point to the list
			nl_test_temp2.push_back(nl_test_temp[index]->clone());
			nl_reference_temp2.push_back(nl_reference_temp[index]->clone());

			np4++;
		}

		//Convert amount of analyzed features to char
		char points_text[5];
		sprintf(points_text, "%d", n_points );

		//Create file names
		char output_file1[MAX_TEXT_LENGTH], output_file2[MAX_TEXT_LENGTH], *output_file3 = "results_amount1.log", *output_file4 = "results_amount2.log";
		strcpy(output_file1, "test"); strcpy(output_file2, "reference");
		strcat(output_file1, points_text); strcat(output_file2, points_text);
		strcat(output_file1, ".txt"); strcat(output_file2, ".txt");

		//Open files
		std::filebuf fb1, fb2, fb3, fb4;
		fb1.open(output_file1, std::ios::out); fb2.open(output_file2, std::ios::out);
		fb3.open(output_file3, std::ios::app); fb4.open(output_file4, std::ios::app);
		std::ostream os1(&fb1); std::ostream os2(&fb2);
		std::ostream os3(&fb3); std::ostream os4(&fb4);
		
		//Write to files
		nl_test_temp2.print(&os1);
		nl_reference_temp2.print(&os2);
		
		//Close files
		fb1.close();
		fb2.close();

		//Perform analysis
		sl.clear();
		computeAnalysisForAllSamplesNLS(sl, proj_list, nl_test_temp2, nl_reference_temp2, meridians, parallels, faces_test, analysis_parameters,
			total_created_or_thrown_samples, output);

		//CartAnalysis::computeAnalysisForAllSamplesDE(sl, proj_list, nl_test_temp2, nl_reference_temp2, meridians, parallels, faces_test, analysis_parameters,
		//	total_created_or_thrown_samples, output);

		//Sort computed results
		output->flush();
		std::cout << ">> Sorting all samples...";
		CartAnalysis::sortSamplesByComputedRatios ( sl, analysis_parameters.analysis_type );
		std::cout << " Completed."  << std::endl << std::endl;

		//Print all results into output
		std::cout << "Print results";
	
		CartAnalysis::printResults2(sl, nl_test_temp2, nl_reference_temp2, analysis_parameters, output);

		//List of best candidates
		CartAnalysis::printResults3(sl, nl_test_temp2, nl_reference_temp2, analysis_parameters, 0, &os3);
		
		//List of second best candidates
		CartAnalysis::printResults3(sl, nl_test_temp2, nl_reference_temp2, analysis_parameters, 1, &os4);

		//Close files
		fb3.close();
		fb4.close();

		//Assign new list to the old one
		nl_test_temp = nl_test_temp2;
		nl_reference_temp = nl_reference_temp2;

		//Clear lists
		nl_test_temp2.clear();
		nl_reference_temp2.clear();

		//Delete min-max box
		delete m1;
		delete m2;
		delete m3;
		delete m4;

		//Increment iterations
		iterations++;
	}
}



template <typename T>


void CartAnalysis::batchTestG(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, const TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
	unsigned short &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, const unsigned int n, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, Container <Sample <T> > &sl, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output)
{
	unsigned short iterations_tot = 0, effic = 0, n_tests = 300;
	T cost_tot = 0, cost_good = 0;
	srand((unsigned)time(0));

	const TMinMax <T> lon_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLon()))->getLon());
	const TMinMax <T> lat_interval((*std::min_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat(),
		(*std::max_element(pl_reference.begin(), pl_reference.end(), sortPointsByLat()))->getLat());

	TMinMax <T> latp_interval_heur = proj->getLatPIntervalH(lat_interval);
	TMinMax <T> lonp_interval_heur = proj->getLonPIntervalH(lon_interval);
	TMinMax <T> lat0_interval = proj->getLat0Interval();

	//Remember old values
	Point3DGeographic <T> cart_pole = proj->getCartPole();
	const T lat0 = proj->getLat0();
	const T lon0 = proj->getLon0();
	const T dx = proj->getDx();
	const T dy = proj->getDy();
	const T c = proj->getC();

	//Start time
	time_t start, end;
	time(&start);
	const int m = nl_test.size();

	for (unsigned int i = 0; i < n_tests; i++)
	{
		Matrix <T> X(n, 1), Y(2 * m, 1);
		Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1);

		//Intervals initialization
		const T R_min = 3000, R_max = 10000;
		const T latp_min = -70, latp_max = 70, lonp_min = -150, lonp_max = 150, lat0_max = 80, lat0_min = 10;
		const T rR = R_max - R_min + 1;
		const T rlatp = latp_max - latp_min + 1;
		const T rlonp = lonp_max - lonp_min + 1;
		const T rlat0 = lat0_max - lat0_min + 1;

		//Random numbers
		const T R_r = R_min + int(rR*rand() / (RAND_MAX + 1.0));
		const T latp_r = latp_min + int(rlatp * rand() / (RAND_MAX + 1.0));
		const T lonp_r = lonp_min + int(rlonp * rand() / (RAND_MAX + 1.0));
		const T lat0_r = lat0_min + int(rlat0 * rand() / (RAND_MAX + 1.0));

		//Store randomly generated values
		X(0, 0) = 6380;
		X(1, 0) = (aspect == NormalAspect ? 90 : latp_r);
		X(2, 0) = (aspect == NormalAspect ? 0 : lonp_r);
		X(3, 0) = lat0_r;
		X(5, 0) = 1.0;

		/*
		X(0, 0) = 6539;
		X(1, 0) = 59;
		X(2, 0) = -3;
		X(3, 0) = 77;
		X(5, 0) = 1.0;
		*/
		/*
		X(0, 0) = 7254.0000000;
		X(1, 0) = 29.0000000;
		X(2, 0) = -29.0000000;
		X(3, 0) = 23.0000000;
		X(5, 0) = 1.0;
		*/
		//Print actual values
		std::cout << "Test " << i << "/" << n_tests << '\n';
		std::cout << "R= " << R_r << "   latp= " << latp_r << "   lonp= " << lonp_r << "   lato0= " << lat0_r << '\n';

		//Clear list of points
		nl_projected.clear();

		//Perform analysis
		T min_cost = 0;
		min_cost = (analysis_parameters.analysis_method == NonLinearLeastSquaresMethod ?
			NonLinearLeastSquares::BFGSH(FAnalyzeProjJ2 <T>(nl_test, pl_reference, proj, aspect, analysis_parameters.print_exceptions), FAnalyzeProjV2 <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
			analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, iterations, alpha_mult, nu, eps, max_iter, max_diff, output) :
			NonLinearLeastSquares::BFGSH(FAnalyzeProjJ3 <T>(nl_test, pl_reference, nl_projected, proj, aspect, x_mass_reference, y_mass_reference, analysis_parameters.print_exceptions), FAnalyzeProjV3 <T>(nl_test, pl_reference, nl_projected, meridians, parallels, faces_test, proj,
			x_mass_reference, y_mass_reference, analysis_parameters, aspect, best_sample, total_created_and_analyzed_samples_projection, output), FAnalyzeProjC <double>(), W, X, Y, V, iterations, alpha_mult, nu, eps, max_iter, max_diff, output));


		if (n >= 6)
		{
			//Set properties of the sample
			best_sample.setR(X(0, 0));
			best_sample.setLatP(X(1, 0));
			best_sample.setLonP(X(2, 0));
			best_sample.setLat0(X(3, 0));
			best_sample.setLon0(X(4, 0));
		}

		else
		{
			best_sample.setLatP(X(0, 0));
			best_sample.setLonP(X(1, 0));
			best_sample.setLat0(X(2, 0));
			best_sample.setLon0(X(3, 0));
		}


		if (analysis_parameters.analysis_method == NonLinearLeastSquaresRotMethod)
			best_sample.setAlpha(X(6, 0));


		//Add to the list
		sl.push_back(best_sample);

		//Restore projection properties after analysis
		proj->setCartPole(cart_pole);
		proj->setLat0(lat0);
		proj->setLon0(lon0);
		proj->setDx(dx);
		proj->setDy(dy);
		proj->setC(c);

		//Test result
		iterations_tot += iterations;
		cost_tot += norm(trans(V) * W * V);

		//Correct result
		//if ( best_sample.getHomotheticTransformationRatio() < 5.0e6)
		//1.0e-4
		//3e0
		//5e6
		//3.0e-3
		//3.0e-2
		std::cout << "\nC= " << min_cost << "   ";
		if (min_cost < 3.0e0)
		{
			effic++;
			cost_good += norm(trans(V) * W * V);
		}

		//Bad result
		else
		{
			*output << "Bad convergence, latp: " << latp_r << ",  lonp:" << lonp_r << ",  lat0:" << lat0_r << '\n';
			std::cout << "Bad convergence, latp: " << latp_r << ",  lonp:" << lonp_r << ",  lat0:" << lat0_r << '\n';

			X.print(output);
			X.print();
		}

		std::cout << i << " ";
	}

	//Time difference
	time(&end);
	float time_diff = difftime(end, start);

	//Print results
	*output << "***** RESULTS ***** \n" << '\n';
	*output << "Efficiency: " << effic << '\n';
	*output << "Iterations: " << iterations_tot << '\n';
	*output << "Cost good: " << cost_good << '\n';
	*output << "Cost total: " << cost_tot << '\n';
	*output << "Time [s]: " << time_diff << '\n';
}



#endif