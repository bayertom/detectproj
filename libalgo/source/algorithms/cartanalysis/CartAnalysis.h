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

#ifndef CartAnalysis_H
#define CartAnalysis_H

#include <list>

#include "libalgo/source/structures/point/Point3DGeographic.h"

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/transformation/Transformation2D.h"

#include "libalgo/source/algorithms/geneticalgorithms/DifferentialEvolution.h"

//Forward declarations
template <typename T >
class Node3DCartesian;

template <typename T >
class Node3DCartesianProjected;

template <typename T >
class Point3DGeographic;

template <typename T>
class Projection;

template <typename T>
class Projection;

template <typename T >
class Sample;

template <typename T >
struct TRansacResults;

template <typename T >
struct TMeridiansList;

template <typename T >
struct TParallelsList;

template <typename T>
class sortSamplesByAllRatios;

template <typename T>
class FAnalyzeProjJ;

template <typename T>
class FAnalyzeProjJ2;

template <typename T>
class FAnalyzeProjJ3;

template <typename T>
class FAnalyzeProjJ4;

template <typename T>
class FAnalyzeProjV;

template <typename T>
class FAnalyzeProjVS;

template <typename T>
class FAnalyzeProjVS2;

template <typename T>
class FAnalyzeProjVDE;

template <typename T>
class FAnalyzeProjV2;

template <typename T>
class FAnalyzeProjV2S;

template <typename T>
class FAnalyzeProjV2S2;

template <typename T>
class FAnalyzeProjV2DE;

template <typename T>
class FAnalyzeProjV3;

template <typename T>
class FAnalyzeProjV3S;

template <typename T>
class FAnalyzeProjV3S2;

template <typename T>
class FAnalyzeProjV3DE;

template <typename T>
class FAnalyzeProjV4;

template <typename T>
class FAnalyzeProjV4S;

template <typename T>
class FAnalyzeProjV4S2;

template <typename T>
class FAnalyzeProjV4DE;

template <typename T>
class FAnalyzeProjC;

template <typename T>
class removeProjectionPolePositions;

class sortProjectionPolePositionsByCompCriterium;
class sortProjectionPolePositionsByLat;


//Coordinates of the cartographic pole + complex criterium value
template <typename T>
struct TProjectionPolePosition
{
        Point3DGeographic <T> cart_pole;
        T lat0;
        T complex_crit;

        TProjectionPolePosition () :  cart_pole ( MAX_LAT, 0.0 ), lat0 ( ), complex_crit ( 0.0 ) {}
        TProjectionPolePosition ( const T latp_, const T lonp_, const T lat0_, const T complex_crit_ ) :  cart_pole ( latp_, lonp_ ), lat0 ( lat0_ ),  complex_crit ( complex_crit_ ) {}
};


//Parameters of the analyzed projection extracted from command line
template <typename T>
struct TAnalyzedProjParameters
{
        char proj_name[64];
        T lat0, latp, lonp;

        //Set bad values in default constructor to distinguish, if a parameter was set from command line or not
        TAnalyzedProjParameters () : lat0 ( 1000.0 ), latp ( 1000.0 ), lonp ( 1000.0 ) {}
};


//List of analyzed projections extracted from command line
template <typename T>
struct TAnalyzedProjParametersList
{
        typedef std::list < TAnalyzedProjParameters <T> > Type;
};


//Result of absolute analysis will be matched using circle or Tissot Indicatrix
typedef enum
{
        MatchCircle = 0,
        MatchTissotIndicatrix,
} TMatchPointsType;


//Collect matched points to the list
typedef enum
{
        CollectOn = 0,
        CollectOff,
} TMatchedPointsCollectType;


//Aspect of a projection used to generate all appropriate pole positions
typedef enum
{
        NormalAspect = 0,			//Projection in the normal aspect
        TransverseAspect,			//Projection in the transverse aspect
        ObliqueAspect,				//Projection in the oblique aspect
} TProjectionAspect;


//Method of the analysis
typedef enum
{
        SimplexMethod = 0,			//Simplex method without rotation (global optmization)
	DifferentialEvolutionMethod,		//Differential evolution without rotation (global optmization)
	NonLinearLeastSquaresMethod,		//Non linear least squares solution  (local optmization)
	SimplexRotMethod,			//Simplex method involving map rotation (global optmization)
	DifferentialEvolutionRotMethod,		//Differential evolution involving map rotation (global optmization)
	NonLinearLeastSquaresRotMethod,		//Non linear least squares solution involving map rotation (local optmization)
	SimplexRot2Method,			//Simplex method involving map rotation (global optmization)
	DifferentialEvolutionRot2Method,	//Differential evolution involving map rotation (global optmization)
	NonLinearLeastSquaresRot2Method,	//Non linear least squares solution involving map rotation (local optmization), scaled
	SimplexShiftsMethod,			//Simplex method involving map rotation (global optmization)
	DifferentialEvolutionShiftsMethod,	//Differential evolution involving map rotation (global optmization)
	NonLinearLeastSquaresShiftsMethod

} TAnalysisMethod;


// Parameters of the projection usedin the incremental algorithm
template <typename T>
struct TSampleProjection
{
	T weight;
	Projection <T> * projection;
	TProjectionAspect aspect;
	Matrix <T> X, XMIN, XMAX;

	TSampleProjection( const unsigned int n) : weight(0.0), projection(NULL), aspect (NormalAspect), X(n, 1), XMIN(n, 1), XMAX(n, 1) {}
};


//Parametrs of the analysis get from the command line
template <typename T>
struct TAnalysisParameters
{
        //Internal structure: Information which cartometric analysis will be performed
        struct TAnalysisType
        {
                //Performed analysis
                bool a_cnd, a_homt, a_helt, a_gn_tf, a_vd_tf, a_vd_id;

                //Constructors
                TAnalysisType () :  a_cnd ( false ), a_homt ( false ), a_helt ( false ), a_gn_tf ( false ), a_vd_tf ( false ), a_vd_id ( false ) {}

				TAnalysisType(const bool status) : a_cnd(status), a_homt(status), a_helt(status), a_gn_tf(status), a_vd_tf(status), a_vd_id (status) {}

                TAnalysisType ( const bool a_cnd_, const bool a_homt_, const bool a_helt_, const bool a_gn_tf_, const bool a_vd_tf_ , const bool a_vd_id_) :
						a_cnd(a_cnd_), a_homt(a_homt_), a_helt(a_helt_), a_gn_tf(a_gn_tf_), a_vd_tf(a_vd_tf_), a_vd_id(a_vd_id_) {}
        };

        //Analysis type
        TAnalysisType analysis_type;

        TAnalysisMethod analysis_method;

        //Compare results to circle or Tissot Indicatrix
        TMatchPointsType match_method;

        //Switches
        bool perform_heuristic, analyze_normal_aspect, analyze_transverse_aspect, analyze_oblique_aspect, remove_outliers, correct_rotation, print_exceptions;

        //Total printed results
        unsigned short exported_graticule, printed_results, analysis_repeat;

        //Parameters of the analysis: steps, increments, intervals and sensitivity
        T lon0, lat_step, lon_step, latp_step, lonp_step, lat0_step, analyzed_proj_lat0, heuristic_sensitivity_ratio, heuristic_sensitivity_increment, max_error;

        //Analyzed projections
        Container <Projection <T> *> analyzed_projections;

        //Test and reference files
        char * test_file, *reference_file, *projections_file;

        TAnalysisParameters () : analysis_type ( false ), analysis_method ( NonLinearLeastSquaresMethod ), match_method ( MatchTissotIndicatrix ), perform_heuristic ( false ), analyze_normal_aspect ( false ), analyze_transverse_aspect ( false ), analyze_oblique_aspect ( false ),
                remove_outliers ( false ), correct_rotation ( false ), print_exceptions ( false ), exported_graticule ( 3 ), printed_results ( 50 ), analysis_repeat ( 0 ), lon0 ( 0.0 ), lat_step ( 10.0 ), lon_step ( 10.0 ), latp_step ( 10.0 ), lonp_step ( 10.0 ),
                lat0_step ( 10.0 ), heuristic_sensitivity_ratio ( 3.0 ), heuristic_sensitivity_increment ( 2.0 ), max_error ( 1.0 ), analyzed_projections ( ),
                test_file ( NULL ), reference_file ( NULL ), projections_file ( "projections.txt" ) { };

        TAnalysisParameters ( const bool status ) : analysis_type ( status ), analysis_method ( NonLinearLeastSquaresMethod ), match_method ( MatchTissotIndicatrix ), perform_heuristic ( status ), analyze_normal_aspect ( status ), analyze_transverse_aspect ( status ), analyze_oblique_aspect ( status ),
                remove_outliers ( status ), correct_rotation ( status ), print_exceptions ( status ), exported_graticule ( 3 ), printed_results ( 50 ), analysis_repeat ( 0 ), lon0 ( 0.0 ), lat_step ( 10.0 ), lon_step ( 10.0 ), latp_step ( 10.0 ), lonp_step ( 10.0 ),
                lat0_step ( 10.0 ), heuristic_sensitivity_ratio ( 3.0 ), heuristic_sensitivity_increment ( 2.0 ), max_error ( 1.0 ), analyzed_projections ( ),
                test_file ( NULL ), reference_file ( NULL ), projections_file ( "projections.txt" ) { }
};



//Cartometric analysis
class CartAnalysis
{
        public:
                template <typename T>
                static void computeAnalysisForAllSamplesGS ( Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
                                typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
                                unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout );

                template <typename T>
                static void computeAnalysisForAllSamplesSim ( Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
                                typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
                                unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout );


                template <typename T>
                static void computeAnalysisForAllSamplesDE ( Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
                                typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
                                unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout );

                template <typename T>
                static void computeAnalysisForAllSamplesNLS ( Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
                                typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
                                unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout );

		//template <typename T>
		//static void analyzeProjectionIncr(Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
		//	typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
		//	unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout);

                template <typename T>
                static bool checkSample ( const typename TMeridiansList <T> ::Type & meridians, const typename TParallelsList <T> ::Type & parallels, const Container <Node3DCartesian <T> *> &nl_test,
                                          const Container <Node3DCartesianProjected <T> *> & nl_projected, const T treshold = 0.5 );

                template <typename T>
                static void findLatPLonPIntervals ( const Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TMinMax <T> &latp_interval_heur, TMinMax <T> &lonp_interval_heur );

                template <typename T>
                static void createOptimalLatPLonPPositions ( const Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, const TMinMax <T> &latp_interval, const TMinMax <T> &lonp_interval, const TAnalysisParameters <T> & analysis_parameters,
                                const TProjectionAspect & proj_aspect, typename TItemsList <TProjectionPolePosition<T> >::Type & proj_pole_positions_list, std::ostream * output = &std::cout );

                template <typename T>
                static T computeAnalysisForOneSample ( Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
                                                       const Container <Face <T> *> &faces_test, Projection <T> *proj, const TAnalysisParameters <T> & analysis_parameters, Sample <T> &sample_res, bool singular_points_found, unsigned int & created_samples, std::ostream * output );

                template <typename T>
                static void printResults ( const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
                                           const TAnalysisParameters <T> & analysis_parameters, std::ostream * output = &std::cout );

		template <typename T>
		static void printResults2(const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
			const TAnalysisParameters <T> & analysis_parameters, std::ostream * output = &std::cout);

		template <typename T>
		static void printResults3(const Container <Sample <T> > &sl, const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference,
			const TAnalysisParameters <T> & analysis_parameters, unsigned int index, std::ostream * output = &std::cout);

                template <typename T>
                static void sortSamplesByComputedRatios ( Container <Sample <T> > &sl, const typename TAnalysisParameters <T>::TAnalysisType & analysis_type );

                template <typename T>
                static T getMatchRatioTissotIndicatrix ( const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference, const Matrix <T> &X, Matrix <T> &W, const Projection <T> *proj, TIndexList & matched_points_ind,
                                const TMatchedPointsCollectType collect_matched_points, const T map_scale, const T lambda = 0.5 );

                template <typename Point1, typename Point2>
                static typename Point1::Type getMatchRatioCircle ( const Container <Point1 *> &global_points, const Container <Point2 *> &transformed_points, TIndexList & matched_points,
                                const TMatchedPointsCollectType collect_matched_points, const typename Point1::Type radius, const typename Point1::Type lambda = 0.5 );

                template <typename Point1, typename Point2>
                static typename Point1::Type getMatchRatioTissotIndicatrix ( const Container <Point1 *> &global_points, const Container <Point2 *> &transformed_points, TIndexList & matched_points,
                                const TMatchedPointsCollectType collect_matched_points, const typename Point1::Type radius, const typename Point1::Type lambda = 0.5 );

                template <typename T>
                static T solutionDiversity ( Container <Sample <T> > &sl, Container <Projection <T> *> &pl, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, const unsigned int n_items );

		template <typename T>
		static void batchTestAmountOfPoints(Container <Sample <T> > &sl, Container <Projection <T> *> &proj_list, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,
			typename TMeridiansList <T> ::Type meridians, typename TParallelsList <T> ::Type parallels, const Container <Face <T> *> &faces_test, TAnalysisParameters <T> & analysis_parameters,
			unsigned int & total_created_or_thrown_samples, std::ostream * output = &std::cout);

        private:

                template <typename T>
                static void redLon ( const Container <Point3DGeographic <T> *> & pl_source, const T lon0, Container <Point3DGeographic <T> *> & pl_destination );

                template <typename T>
                static void redLon ( const Container <Point3DGeographic <T> *> & pl_source, const T lon0 );

				template <typename T>
				static T getInitialRadius(const Projection <T> *proj, const T scale, const Container <Point3DGeographic <T> *> & pl_reference );


                template <typename T>
                static void removeSingularPoints ( const Container <Node3DCartesian <T> *> & nl_source, const Container <Point3DGeographic <T> *> & pl_source, const Projection <T>  *proj, Container <Node3DCartesian <T> *> &nl_destination, Container <Point3DGeographic <T> *> &pl_destination,
                                                   typename TDevIndexPairs <T>::Type & non_singular_point_pairs );

                template <typename T>
                static void removeSingularPoints ( const Container <Node3DCartesian <T> *> & nl_source, const Container <Point3DGeographic <T> *> & pl_source, const Projection <T>  *proj, typename TDevIndexPairs <T>::Type & non_singular_point_pairs );


                template <typename T>
                static void correctMeridiansAndParrallels ( typename TMeridiansList <T> ::Type & meridians, typename TParallelsList <T> ::Type & parallels, typename TDevIndexPairs <T>::Type & point_pairs );


                template <typename T>
                static void correctPointsMeridiansAndParrallels ( const Container <Node3DCartesian <T> *> &nl_test_corr, const Container <Point3DGeographic <T> *> &pl_reference_corr, const typename TMeridiansList <T> ::Type & meridians, const typename TParallelsList <T> ::Type & parallels, const unsigned int n,
                                Container <Node3DCartesian <T> *> **p_nl_test, Container <Point3DGeographic <T> *> **p_pl_reference, typename TMeridiansList <T> ::Type ** p_meridians_new, typename TParallelsList <T> ::Type ** p_parallels_new, typename TMeridiansList <T> ::Type & meridians_new,
                                typename TParallelsList <T> ::Type & parallels_corr, typename TDevIndexPairs<T>::Type & point_pairs_new, bool &uncorrect_points_found );

                template <typename T>
                static void analyzeSampleCrossNearestNeighbourDistance ( Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const float mult_ratio );

                template <typename T>
                static void analyzeSampleHomotheticTransformationDeviation ( Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const TMatchPointsType match_type, const float mult_ratio );

                template <typename T>
                static void analyzeSampleHelmertTransformationDeviation ( Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const TMatchPointsType match_type, const float mult_ratio );

                template <typename T>
                static void analyzeSampleGeographicNetworkTurningFunctionRatio ( Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test,  const Container <Node3DCartesianProjected  <T> *> &nl_projected,
                                const typename TMeridiansList <T> ::Type & meridians, const typename TParallelsList <T> ::Type & parallels, const float mult_ratio );

                template <typename T, TDestructable destructable>
                static void analyzeSampleUsingVoronoiDiagramTurningFunctionRatio2 ( Sample <T> &s, const Container <Node3DCartesian <T> *> &nl_test, const Container <Node3DCartesianProjected  <T> *> &nl_projected, const Container <Face <T> *, destructable> &faces_test,
                                const TAnalysisParameters <T> & analysis_parameters, const float mult_ratio );

                template <typename T>
                static void setPositionForSortedSamples ( Container <Sample <T> > &sl, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type );

		template <typename T>
		static void batchTestShifts(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
			unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> &best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);

		template <typename T>
		static void batchTestNLSP( Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,  Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect,  typename TMeridiansList <T> ::Type &meridians,  typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
			const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);


		template <typename T>
		static void batchTestSimplex(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference,  Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init, 
			const T lonp_init, const T lat0_init, unsigned short &iterations, const T eps, unsigned short max_iter, const T max_diff, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);

		template <typename T>
		static void batchTestDiffEvolutionSchema(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init,
			const T lonp_init, const T lat0_init, unsigned int &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);

		template <typename T>
		static void batchTestDiffEvolutionAdaptiveControl(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init,
			const T lonp_init, const T lat0_init, unsigned int &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, const TMutationStrategy strategy, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);
		
		template <typename T>
		static void batchTestDiffEvolutionOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test, const T R0, const T latp_init,
			const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, const TMutationStrategy strategy, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output);

		template <typename T>
		static void batchTestNLSPOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
			const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);

		template <typename T>
		static void batchTestNelderMeadOutliers(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
			const T R0, const T latp_init, const T lonp_init, const T lat0_init, unsigned int  &iterations, const T eps, unsigned short max_iter, const T max_diff, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);

		template <typename T>
		static void batchTestG(Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, Container <Node3DCartesianProjected <T> *> &nl_projected, Projection <T> *proj, const TProjectionAspect aspect, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels, const Container <Face <T> *> &faces_test,
			unsigned short  &iterations, const T alpha_mult, const T nu, const T eps, unsigned short max_iter, const T max_diff, T &x_mass_reference, T &y_mass_reference, const unsigned int n, Sample <T> best_sample, TAnalysisParameters <T> & analysis_parameters, Container <Sample <T> > &sl, unsigned int &total_created_and_analyzed_samples_projection, std::ostream * output = &std::cout);



};

#include "CartAnalysis.hpp"

#endif
