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

#ifndef Sample_H
#define Sample_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

//Forward declarations
struct TAnalysisType;

template <typename T>
class Node3DCartesian;

template <typename T>
class Projection;


//Cartographic sample in specified projection for ccrtometric analysis
template <typename T>
class Sample
{
        private:
                //Map projection (overriden in batch processing)
                Projection <T> * proj;								//Map pojection related to the sample

                //Map projection properties
                T R;										//Sphere radius
                T latp;										//Latitude of the cartographic pole, overriden compared to the projection definition.
                T lonp;										//Longitude of the cartographic pole, overriden compared to the projectin definition.
                T lat0;										//Latitude of the undistorted parallel, overriden compared to the projection definition.
		T lat1;										//Latitude of the undistorted parallel, overriden compared to the projection definition.
		T lat2;										//Latitude of the undistorted parallel, overriden compared to the projection definition.
                T lon0;										//New central parallel
                T dx;										//Additive constant
                T dy;										//Additive constant
                T c;										//Any other constant of the projection needs to be minimized
		T alpha;									//Rotation estimated by M7 method

                //Estimated cartographic properties
                T scale_homt, scale_helt;							//Estimated scale of analyzed map
                T rotation;									//Estimated rotation of analyzed map / old rotation for corrected samples
		
                bool analyzed_proj_sample;							//Represents this sample an analyzed sample with known projection set by user?
                bool rotated_sample;								//Was the sample rotated?
                bool singular_points_found;							//Was any singular point found in reference file
                bool outliers_found;								//Was any outlier found in files

                //Results of analysis
                T cross_nearest_neighbour_distance_ratio;					//Cross nearest neighbour distance difference between two datasets
                T homothetic_transformation_ratio;						//Homothetic transformation ratio
                T helmert_transformation_ratio;							//Helmert transformation ratio
                T gn_turning_function_ratio;							//Turning function ratio for meridians and parallels in geographic network
                T voronoi_cell_turning_function_ratio;						//Voronoi cell turning function ratio
		T voronoi_cell_inner_distance_ratio;						//Voronoi cell inner distance + shape context ratio

		unsigned int iterations;							//Amount of iterations
		unsigned int residual_eval;							//Amount of residual evaluations

                //Percentage match
                unsigned int homothetic_transformation_perc_match;				//Percantage match ratio for homothetic transformation
                unsigned int helmert_transformation_perc_match;					//Percantage match ratio for helmert transformation

                //List of k-best points
                TIndexList non_singular_points_indices;						//List of acceptable points indices after removing singular points
                TIndexList k_best_points_indices;						//List of acceptable points indices after removing outliers

                //Matched points
                TIndexList homothetic_transformation_matched_points_indices;			//Indices of matched points: homothehtic transformation
                TIndexList helmert_transformation_matched_points_indices;			//Indices of matched points: helmert transformation

                //Positions of a sample by computed ratios
                int cross_nearest_neighbour_distance_ratio_position;				//Position of the sample sorted by the cross nearest distance
                int homothetic_transformation_ratio_position;					//Position of the sample sorted by the homothetic transformation
                int helmert_transformation_ratio_position;					//Position of the sample sorted by the helmert transformation
                int gn_turning_function_ratio_position;						//Position of the sample sorted by the geographic network tangent function
                int voronoi_cell_turning_function_ratio_position;				//Position of the sample sorted by the voronoi cell tangent function
		int voronoi_cell_inner_distance_ratio_position;
				
        public:
		Sample() : proj(NULL), R(1.0), latp(MAX_LAT), lonp(0.0), lat0(0.0), lat1(0.0), lat2(0.0), lon0(0.0), dx(0.0), dy(0.0), c(0.0), alpha(0.0), scale_homt(1.0), scale_helt(1.0), rotation(0.0), analyzed_proj_sample(false), rotated_sample(false), singular_points_found(false), outliers_found(false),
			cross_nearest_neighbour_distance_ratio(0.0), homothetic_transformation_ratio(0.0), helmert_transformation_ratio(0.0), gn_turning_function_ratio(0.0), voronoi_cell_turning_function_ratio(0.0), voronoi_cell_inner_distance_ratio(0.0), iterations (0), residual_eval(0), homothetic_transformation_perc_match(0), helmert_transformation_perc_match(0),
                        non_singular_points_indices ( 0 ), k_best_points_indices ( 0 ), homothetic_transformation_matched_points_indices ( 0 ), helmert_transformation_matched_points_indices ( 0 ), cross_nearest_neighbour_distance_ratio_position ( 1 ), homothetic_transformation_ratio_position ( 1 ),
                        helmert_transformation_ratio_position ( 1 ), gn_turning_function_ratio_position ( 1 ), voronoi_cell_turning_function_ratio_position ( 1 ) {}

		Sample(Projection <T> * proj_, const T R_, const T latp_, const T lonp_, const T lat0_, const T lat1_, const T lat2_, const T lon0_, const T dx_, const T dy_, const T c_, const T alpha_, const bool singular_points_found_, const bool outliers_found_, const TIndexList & non_singular_points_indices_, const TIndexList k_best_points_indices_) : proj(proj_),
			R(R_), latp(latp_), lonp(lonp_), lat0(lat0_), lat1(lat1_), lat2(lat2_), lon0(lon0_), dx(dx_), dy(dy_), c(c_), alpha(alpha_), scale_homt(1.0), scale_helt(1.0), rotation(0.0), analyzed_proj_sample(false), rotated_sample(false), singular_points_found(singular_points_found_), outliers_found(outliers_found_),
			cross_nearest_neighbour_distance_ratio(0.0), homothetic_transformation_ratio(0.0), helmert_transformation_ratio(0.0), gn_turning_function_ratio(0.0), voronoi_cell_turning_function_ratio(0.0), voronoi_cell_inner_distance_ratio(0.0), iterations(0), residual_eval(0), homothetic_transformation_perc_match(0), helmert_transformation_perc_match(0),
                        non_singular_points_indices ( non_singular_points_indices_ ), k_best_points_indices ( k_best_points_indices_ ), homothetic_transformation_matched_points_indices ( 0 ), helmert_transformation_matched_points_indices ( 0 ), cross_nearest_neighbour_distance_ratio_position ( 1 ), homothetic_transformation_ratio_position ( 1 ),
                        helmert_transformation_ratio_position ( 1 ), gn_turning_function_ratio_position ( 1 ), voronoi_cell_turning_function_ratio_position ( 1 ) {}

                virtual ~Sample();

        public:
                Projection <T> * getProj() const {return proj;}
                T getR() const {return R;}
                T getLatP() const {return latp;}
                T getLonP() const {return lonp;}
                T getLat0() const {return lat0;}
		T getLat1() const { return lat1; }
		T getLat2() const { return lat2; }
                T getLon0() const {return lon0;}
                T getDx() const {return dx;}
                T getDy() const {return dy;}
                T getC() const {return c;}
                T getScaleHomT() const {return scale_homt;}
                T getScaleHelT() const {return scale_helt;}
                T getRotation() const {return rotation;}
		T getAlpha() const {return alpha;}
                bool getAnalyzedProjectionSample() const {return analyzed_proj_sample;}
                bool getRotatedSample() const {return rotated_sample;}
                bool getSingularPointsFound() const {return singular_points_found;}
                bool getOutliersFound() const {return outliers_found;}
                T getCrossNearestNeighbourDistanceRatio() const {return cross_nearest_neighbour_distance_ratio;}
                T getHomotheticTransformationRatio() const {return homothetic_transformation_ratio;}
                T getHelmertTransformationRatio() const {return helmert_transformation_ratio;}
                T getGNTurningFunctionRatio() const {return gn_turning_function_ratio;}
                T getVoronoiCellTurningFunctionRatio() const {return voronoi_cell_turning_function_ratio;}
		T getVoronoiCellInnerDistanceRatio() const { return voronoi_cell_inner_distance_ratio; }
		unsigned int getIterations() const { return iterations; }
		unsigned int getResidualEval() const { return residual_eval; }

                unsigned int getHomotheticTransformationPercMatch() const { return homothetic_transformation_perc_match; }
                unsigned int getHelmertTransformationPercMatch() const { return helmert_transformation_perc_match; }
                TIndexList const & getNonSingularPointsIndices () const {return non_singular_points_indices;}
                TIndexList const & getKBestPointsIndices () const {return k_best_points_indices;}
                TIndexList const & getHomotheticTransformationMatchedPointsIndices () const {return homothetic_transformation_matched_points_indices;}
                TIndexList const & getHelmertTransformationMatchedPointsIndices () const {return helmert_transformation_matched_points_indices;}
                T getSampleCost ( const typename TAnalysisParameters <T>::TAnalysisType & analysis_type ) const;

        public:
                int getCrossNearestNeighbourDistanceRatioPosition() const {return cross_nearest_neighbour_distance_ratio_position;}
                int getHomotheticTransformationRatioPosition() const {return homothetic_transformation_ratio_position;}
                int getHelmertTransformationRatioPosition() const {return helmert_transformation_ratio_position;}
                int getGNTurningFunctionRatioPosition() const {return gn_turning_function_ratio_position;}
                int getVoronoiCellTurningFunctionRatioPosition() const {return voronoi_cell_turning_function_ratio_position;}
				int getVoronoiCellInnerDistanceRatioPosition() const { return voronoi_cell_inner_distance_ratio_position; }

        public:

                void setProj ( Projection <T> * proj_ ) { proj = proj_;}
                void setR ( const T R_ )  {R = R_;}
                void setLatP ( const T latp_ )  {latp = latp_;}
                void setLonP ( const T lonp_ )  {lonp = lonp_;}
                void setLat0 ( const T lat0_ )  {lat0 = lat0_;}
		void setLat1 (const T lat1_)  { lat1 = lat1_; }
		void setLat2( const T lat2_)  { lat2 = lat2_; }
                void setLon0 ( const T lon0_ )  {lon0 = lon0_;}
                void setDx ( const T dx_ )  {dx = dx_;}
                void setDy ( const T dy_ )  {dy = dy_;}
                void setC ( const T c_ )  {c = c_;}
                void setScaleHomT ( const T scale_homt_ ) {scale_homt = scale_homt_;}
                void setScaleHelT ( const T scale_helt_ ) {scale_helt = scale_helt_;}
                void setRotation ( const T rotation_ ) {rotation = rotation_;}
		void setAlpha ( const T alpha_ ) { alpha = alpha_;}
                void setAnalyzedProjectionSample ( const bool analyzed_proj_sample_ ) {analyzed_proj_sample = analyzed_proj_sample_;}
                void setRotatedSample ( const bool rotated_sample_ ) {rotated_sample = rotated_sample_;}
                void setSingularPointsFound ( const bool singular_points_found_ ) {singular_points_found = singular_points_found_;}
                void setOutliersFound ( const bool outliers_found_ ) {outliers_found = outliers_found_;}
                void setCrossNearestNeighbourDistanceRatio ( const T cross_nearest_neighbour_distance_ratio_ ) {cross_nearest_neighbour_distance_ratio = cross_nearest_neighbour_distance_ratio_;}
                void setHomotheticTransformationRatio ( const T homothetic_transformation_ratio_ ) {homothetic_transformation_ratio = homothetic_transformation_ratio_;}
                void setHelmertTransformationRatio ( const T helmert_transformation_ratio_ ) {helmert_transformation_ratio = helmert_transformation_ratio_;}
                void setGNTurningFunctionRatio ( const T turning_function_ratio_ ) {gn_turning_function_ratio = turning_function_ratio_;}
                void setVoronoiCellTurningFunctionRatio ( const T voronoi_cell_turning_function_ratio_ ) {voronoi_cell_turning_function_ratio = voronoi_cell_turning_function_ratio_;}
		void setVoronoiCellInnerDistanceRatio(const T voronoi_cell_inner_distance_ratio_) { voronoi_cell_inner_distance_ratio = voronoi_cell_inner_distance_ratio_; }
		void setIterations(const unsigned int iterations_) { iterations = iterations_; }
		void setResidualEval(const unsigned int residual_eval_) { residual_eval = residual_eval_; }

		void setHomotheticTransformationPercMatch ( const unsigned int homothetic_transformation_perc_match_ ) {homothetic_transformation_perc_match = homothetic_transformation_perc_match_;}
                void setHelmertTransformationPercMatch ( const unsigned int helmert_transformation_perc_match_ ) {helmert_transformation_perc_match = helmert_transformation_perc_match_;}
                void setNonSingularPointsIndices ( const TIndexList & non_singular_points_indices_ ) {non_singular_points_indices = non_singular_points_indices_;}
                void setKBestPointsIndices ( const TIndexList & k_best_points_indices_ ) {k_best_points_indices = k_best_points_indices_;}
                void setHomotheticTransformationMatchedPointsIndices ( const TIndexList & homothetic_transformation_matched_points_indices_ ) {homothetic_transformation_matched_points_indices = homothetic_transformation_matched_points_indices_;}
                void setHelmertTransformationMatchedPointsIndices ( const TIndexList & helmert_transformation_matched_points_indices_ ) {helmert_transformation_matched_points_indices = helmert_transformation_matched_points_indices_;}


        public:
                void setCrossNearestNeighbourDistanceRatioPosition ( const int cross_nearest_neighbour_distance_ratio_position_ ) {cross_nearest_neighbour_distance_ratio_position = cross_nearest_neighbour_distance_ratio_position_;}
                void setHomotheticTransformationRatioPosition ( const int homothetic_transformation_ratio_position_ ) {homothetic_transformation_ratio_position = homothetic_transformation_ratio_position_;}
                void setHelmertTransformationRatioPosition ( const int helmert_transformation_ratio_position_ ) {helmert_transformation_ratio_position = helmert_transformation_ratio_position_;}
                void setGNTurningFunctionRatioPosition ( const int gn_turning_function_ratio_position_ ) {gn_turning_function_ratio_position = gn_turning_function_ratio_position_;}
                void setVoronoiCellTurningFunctionRatioPosition ( const int voronoi_cell_turning_function_ratio_position_ ) {voronoi_cell_turning_function_ratio_position = voronoi_cell_turning_function_ratio_position_;}
				void setVoronoiCellInnerDistanceRatioPosition(const int voronoi_cell_inner_distance_ratio_position_) { voronoi_cell_inner_distance_ratio_position = voronoi_cell_inner_distance_ratio_position_; }



        public:
                virtual void print ( std::ostream * output = &std::cout ) const;
                friend void operator << ( std::ostream & output, const Sample <T> &s ) { s.print ( &output ) ;}
                void printSampleRatios ( const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, const unsigned int n, std::ostream * output ) const;
		void printSampleRatios2 (const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, const unsigned int n, std::ostream * output) const;
                void printSamplePositions ( const int position, const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, std::ostream * output ) const;
                void printSampleMatchedPoints ( const Container < Node3DCartesian <T> *> &nl_test, const Container <Point3DGeographic <T> *> &nl_reference, const int position,
                                                const typename TAnalysisParameters<T>::TAnalysisType & analysis_type, std::ostream * output ) const;
};

#include "Sample.hpp"

#endif
