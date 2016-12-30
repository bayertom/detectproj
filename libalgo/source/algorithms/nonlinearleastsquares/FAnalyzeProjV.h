// Description: Functor, compute matrix V of residuals for cartometric analysis

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


#ifndef FAnalyzeProjV_H
#define FAnalyzeProjV_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV
{
        private:

                Container <Node3DCartesian <T> *> &nl_test;
                Container <Point3DGeographic <T> *> &pl_reference;
                typename TMeridiansList <T> ::Type &meridians;
                typename TParallelsList <T> ::Type &parallels;
                const Container <Face <T> *> &faces_test;
                Projection <T> *proj;
                const TAnalysisParameters <T> &analysis_parameters;
                const TProjectionAspect aspect;
                Sample <T> &sample_res;
                unsigned int & created_samples;
                std::ostream * output;

        public:

                FAnalyzeProjV ( Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
                                const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, std::ostream * output_ )
                        : nl_test ( nl_test_ ), pl_reference ( pl_reference_ ), meridians ( meridians_ ), parallels ( parallels_ ), faces_test ( faces_test_ ),  proj ( proj_ ), analysis_parameters ( analysis_parameters_ ), aspect ( aspect_ ), sample_res ( sample_res_ ),
                          created_samples ( created_samples_ ), output ( output_ ) {}

                void operator () ( Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true )
                {

                        //Compute parameters of the V matrix: residuals
                        const unsigned int m = nl_test.size();

                        //Get lat0 min and lat0 max
                        const T lat0_min = proj->getLat0Interval().min_val;
                        const T lat0_max = proj->getLat0Interval().max_val;

                        //Normal aspect: lat0, lon0
                        if ( aspect == NormalAspect )
                        {
                                //Correct R, lat0, lon0
                                if ( X ( 0, 0 ) < 0.0 ) X ( 0, 0 ) = fabs ( X ( 0, 0 ) );

				//Set lat0 inside the interval
				if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

				//Set lon0
				if (X(4, 0) < MIN_LON)  X(4, 0) = MAX_LON - fmod(X(4, 0), MIN_LON);
				else if (X(4, 0) > MAX_LON)  X(4, 0) = MIN_LON - fmod(X(4, 0), MAX_LON);
				//if (fabs(X(4, 0)) > MAX_LON)  X(4, 0) = fmod(X(4, 0), 90);

                                //Set to interval
                                //if ( X ( 3, 0 ) < lat0_min ) X ( 3, 0 ) = lat0_min;
                                //if ( X ( 3, 0 ) > lat0_max ) X ( 3, 0 ) = lat0_max;
                        }

                        //Transverse aspect: lonp, lat0
                        else  if ( aspect == TransverseAspect )

                        {
                                //Correct R, lonp, lat0
                                if ( X ( 0, 0 ) < 0.0 ) X ( 0, 0 ) = fabs ( X ( 0, 0 ) );

                                //Subtract period
				if (X(2, 0) < MIN_LON)  X(2, 0) = MAX_LON + fmod(X(2, 0), MIN_LON);
				else if (X(2, 0) > MAX_LON)  X(2, 0) = MIN_LON + fmod(X(2, 0), MAX_LON);

				//Set lat0 inside the interval
				if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

				//Set lon0
				X(4, 0) = 0;
			}

                        //Oblique aspect: latp, lonp, lat0
                        else if ( aspect == ObliqueAspect )
                        {
				//Correct R, latp, lonp, lat0
				if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

				//Subtract period
				if (X(1, 0) < MIN_LAT)  X(1, 0) = MIN_LAT - fmod(X(1, 0), MIN_LAT);
				else if (X(1, 0) > MAX_LAT)  X(1, 0) = MAX_LAT - fmod(X(1, 0), MAX_LAT);

				if (X(2, 0) < MIN_LON)  X(2, 0) = MAX_LON + fmod(X(2, 0), MIN_LON);
				else if (X(2, 0) > MAX_LON)  X(2, 0) = MIN_LON + fmod(X(2, 0), MAX_LON);

				//Set lat0 inside the interval
				if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

				//Set lonp to zero, if latp = 90
				if (fabs(X(1, 0) - MAX_LAT) < 3.0)
				{
					//X(1, 0) = 90.0
					//X(2, 0) = 0.0;
				}

				//Set lon0
				X(4, 0) = 0;
                        }

                        //Set properties to the projection: ommit estimated radius, additional constants dx, dy
                        // They will be estimated again using the transformation
                        Point3DGeographic <T> cart_pole ( X ( 1, 0 ), X ( 2, 0 ) );
                        proj->setR ( X( 0, 0) );
                        proj->setCartPole ( cart_pole );
                        proj->setLat0 ( X ( 3, 0 ) );
                        proj->setLon0 ( X ( 4, 0 ) );
			proj->setDx ( X(5, 0) );
			proj->setDy ( X(6, 0) );
                        proj->setC ( X ( 7, 0 ) );

                        //Compute analysis for one sample
                        if ( compute_analysis )
                        {
                                try
                                {
                                        //Compute analysis
                                        try
                                        {
                                                CartAnalysis::computeAnalysisForOneSample ( nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, sample_res, false, created_samples, output );
                                        }

                                        //Throw exception
                                        catch ( Exception & error )
                                        {
                                                if ( analysis_parameters.print_exceptions )
                                                {
                                                        //Print error and info about projection properties
                                                        error.printException ( output );
                                                        *output << "proj = " << proj->getName() << "  latp = " << proj->getCartPole().getLat() << "  lonp = " << proj->getCartPole().getLon() << "  lat0 = " << proj->getLat0() << '\n';
                                                }
                                        }

                                        //Get index list of the sample
                                        TIndexList non_singular_points_indices = sample_res.getNonSingularPointsIndices();
                                        TIndexList k_best_points_indices = sample_res.getKBestPointsIndices();

                                        //Change weights in W matrix: weights of singular points or outliers are 0, otherwise they are 1
                                        unsigned int index_k_best_points = 0, n_k_best = k_best_points_indices.size(), n_points = pl_reference.size();
                                        int index_point = ( n_k_best > 0 ? non_singular_points_indices [ k_best_points_indices [index_k_best_points++] ] : - 1 );

                                        for ( int i = 0; ( i < n_points ) && ( n_k_best > 0 ); i++ )
                                        {
                                                //Set weight of point to 1 (it is not an outlier nor singular)
                                                if ( i == index_point )
                                                {
                                                        W ( index_point, index_point ) = 1.0; W ( index_point + n_points, index_point + n_points ) = 1.0;

                                                        if ( index_k_best_points < n_k_best ) index_point = non_singular_points_indices [ k_best_points_indices [index_k_best_points++] ];
                                                }

                                                //Set weight of point to zero (it is an outlier or singular)
                                                else
                                                {
                                                        W ( i, i ) = 0.0; W ( i + n_points, i + n_points ) = 0.0;
                                                }
                                        }
                                }

                                //Throw error
                                catch ( Exception & error )
                                {
                                        if ( analysis_parameters.print_exceptions ) error.printException();
                                }
                        }

                        //Compute new coordinates
                        for ( unsigned int i = 0; i < m; i++ )
                        {
                                T x_new = 0.0, y_new = 0.0;

				//Get type of the direction
				TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

				//Reduce lon
				const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(4, 0));

				try
				{
					//Convert geographic point to oblique aspect
					T lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0));
					T lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0), trans_lon_dir);

					for (unsigned int j = 0; j < 3; j++)
					{
						try
						{
							//Compute new coordinates: add shifts

							Y(i, 0) = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(7, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
							Y(i + m, 0) = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(6, 0), X(7, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);

							//Y(i, 0) = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(7, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false) + X(5, 0);
							//Y(i + m, 0) = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(7, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false) + X(6, 0);
						}

						//2 attempt to avoid the singularity
						catch (Exception &error)
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
						}
					}
				}

                                //Throw exception: bad conversion, a singular point
                                catch ( Exception & error )
                                {
					//Disable point from analysis: set weight to zero
					W(i, i) = 0; W(i + m, i + m) = 0;
                                }

                                //Compute coordinate differences (residuals): estimated - input
                                V ( i, 0 ) = ( Y ( i, 0 ) - nl_test [i]->getX() );
                                V ( i + m, 0 ) = ( Y ( i + m, 0 ) - nl_test [i]->getY() );
                        }
                }
};


#endif
