// Description: Functor, compute matrix V of residuals for cartometric analysis
// Method: Simplex method, M7 (5 determined parameters, without rotation, without radius)

// Copyright (c) 2010 - 2014
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


#ifndef FAnalyzeProjV4S_H
#define FAnalyzeProjV4S_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV4S
{
	private:

		Container <Node3DCartesian <T> *> &nl_test;
		Container <Point3DGeographic <T> *> &pl_reference;
		typename TMeridiansList <T> ::Type &meridians;
		typename TParallelsList <T> ::Type &parallels;
		const Container <Face <T> *> &faces_test;
		Projection <T> *proj;
		T &R;
		T &q1;
		T &q2;
		const TAnalysisParameters <T> &analysis_parameters;
		const TProjectionAspect aspect;
		Sample <T> &sample_res;
		unsigned int & created_samples;
		unsigned int &res_evaluation;
		TMEstimatorsWeightFunction me_function;
		T k;
		Matrix <unsigned int> &I;
		std::ostream * output;


	public:

		FAnalyzeProjV4S(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_, const Container <Face <T> *> &faces_test_, 
			Projection <T> *proj_, T &R_est_, T &q1_, T &q2_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &res_evaluation_, const TMEstimatorsWeightFunction &me_function_,
			const T k_, Matrix <unsigned int> &I_, std::ostream * output_) : nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R(R_est_), q1(q1_), q2(q2_), analysis_parameters(analysis_parameters_), 
			aspect(aspect_), sample_res(sample_res_), created_samples(created_samples_), res_evaluation(res_evaluation_), me_function(me_function_), k(k_), I(I_), output(output_) {}


		void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
		{

			//Compute parameters of the V matrix: residuals
			const unsigned int m = nl_test.size();
			const unsigned int m1 = X.rows();
			const unsigned int m2 = V.rows();

			//Get lat0 min and lat0 max
			const T lat0_min = proj->getLat0Interval().min_val;
			const T lat0_max = proj->getLat0Interval().max_val;

			//Process all simplex vertices
			for (unsigned int i = 0; i < m1; i++)
			{

				//Normal aspect: lat0, lon0
				if (aspect == NormalAspect)
				{
					//Subtract period
					if (fabs(X(i, 2)) > MAX_LAT) X(i, 2) = fmod(X(i, 2), 90);

					//Set lon0
					if (X(i, 3) < MIN_LON)  X(i, 3) = MAX_LON + fmod(X(i, 3), MIN_LON);
					else if (X(i, 3) > MAX_LON)  X(i, 3) = MIN_LON + fmod(X(i, 3), MAX_LON);
				}

				//Transverse aspect: lonp, lat0
				else  if (aspect == TransverseAspect)
				{
					//Subtract period
					if (X(i, 1) < MIN_LON)  X(i, 1) = MIN_LON - fmod(X(i, 1), MIN_LON);
					else if (X(i, 1) > MAX_LON)  X(i, 1) = MAX_LON - fmod(X(i, 1), MAX_LON);

					//Set lat0 inside the interval
					if (X(i, 2) < lat0_min || X(i, 2) > lat0_max) X(i, 2) = 0.5 * (lat0_min + lat0_max);

					//Set lon0
					X(i, 3) = 0;
				}

				//Oblique aspect: latp, lonp, lat0
				else if (aspect == ObliqueAspect)
				{
					//Subtract period
					if (X(i, 0) < MIN_LAT)  X(i, 0) = MIN_LAT - fmod(X(i, 0), MIN_LAT);
					else if (X(i, 0) > MAX_LAT)  X(i, 0) = MAX_LAT - fmod(X(i, 0), MAX_LAT);

					if (X(i, 1) < MIN_LON)  X(i, 1) = MAX_LON + fmod(X(i, 1), MIN_LON);
					else if (X(i, 1) > MAX_LON)  X(i, 1) = MIN_LON + fmod(X(i, 1), MAX_LON);

					//Set lat0 inside the interval
					if (X(i, 2) < lat0_min || X(i, 2) > lat0_max) X(i, 2) = 0.5 * (lat0_min + lat0_max);

					//Set lonp to zero, if latp = 90
					if (fabs(X(i, 0) - MAX_LAT) < 5)
					{
						//X(i, 0) = 90.0;
						//X(i, 1) = 0.0;
					}

					//Set lon0
					//X(i, 3) = 0;
					if (X(i, 3) < MIN_LON)  X(i, 3) = MAX_LON + fmod(X(i, 3), MIN_LON);
					else if (X(i, 3) > MAX_LON)  X(i, 3) = MIN_LON + fmod(X(i, 3), MAX_LON);
				}

				//Set properties to the projection: ommit estimated radius, additional constants dx, dy
				// They will be estimated again using the transformation
				Point3DGeographic <T> cart_pole(X(i, 0), X(i, 1));
				proj->setR(R);
				proj->setCartPole(cart_pole);
				proj->setLat0(X(i, 2));
				proj->setLat1(X(i, 2));
				proj->setLat2(X(i, 4));
				proj->setLon0(X(i, 3));
				proj->setDx(0.0);
				proj->setDy(0.0);
				proj->setC(X(i, 4));

				W = eye(2 * m, 2 * m, 1.0);

				//Compute coordinate differences (residuals): items of V matrix
				Container <Node3DCartesianProjected <T> *> nl_projected_temp;

				for (unsigned int j = 0; j < m; j++)
				{
					//Get type of the direction
					TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

					//Reduce lon
					const T lon_red = CartTransformation::redLon0(pl_reference[j]->getLon(), X(i, 3));

					T lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0;
					
					try
					{ 
						//Convert geographic point to oblique aspect
						lat_trans = CartTransformation::latToLatTrans(pl_reference[j]->getLat(), lon_red, X(i, 0), X(i, 1));
						lon_trans = CartTransformation::lonToLonTrans(pl_reference[j]->getLat(), lon_red, X(i, 0), X(i, 1), trans_lon_dir);

						for (unsigned int k = 0; k < 3; k++)
						{
							try
							{
								//Compute x, y coordinates
								x = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(i, 4), X(i, 2), X(i, 2), X(i, 4), false);
								y = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(i, 4), X(i, 2), X(i, 2), X(i, 4), false);

								//x = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), X(i, 4), X(i, 2), proj->getLat1(), proj->getLat2(), false);
								//y = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), X(i, 4), X(i, 2), proj->getLat1(), proj->getLat2(), false);
							}

							//2 attempt to avoid the singularity
							catch (Exception &error)
							{
								//Move in latitude direction
								if (k == 0)
								{
									if (lat_trans == MAX_LAT) 
										lat_trans -= GRATICULE_ANGLE_SHIFT;
									else 
										lat_trans += GRATICULE_ANGLE_SHIFT;
								}

								//Move in longitude direction
								else if (k == 1)
								{
									if (lon_trans == MAX_LON) 
										lon_trans -= GRATICULE_ANGLE_SHIFT;
									else 
										lon_trans += GRATICULE_ANGLE_SHIFT;
								}

								//Neither first nor the second shhifts do not bring improvement
								else if (k == 2)
								{
									throw;
								}
							}
						}
					}

					catch (Exception &error)
					{
						//Disable point from analysis: set weight to zero
						W(j, j) = 0; W(j + m, j + m) = 0;
					}

					//Create new cartographic point
					Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

					//Add point to the list
					nl_projected_temp.push_back(n_projected);
				}

				/*
				// Weighted Helmert transformation
				Matrix <T> C(2, 2), beta(4, 1), P(m, 2), Q(m, 2), A(2 * m, 4);

				for (unsigned int j = 0; j < m; j++)
				{
					P(j, 0) = nl_test[j]->getX(); P(j, 1) = nl_test[j]->getY();
					Q(j, 0) = nl_projected_temp[j]->getX(); Q(j, 1) = nl_projected_temp[j]->getY();
				}

				HelmertTransformation2D::getTransformKey2(P, Q, W, A, beta, Y, C);

				//Get centers of gravity
				T x_mass_reference = C(1, 0);
				T y_mass_reference = C(1, 1);
				T x_mass_test = beta(2, 0);		//Determined from the transformation
				T y_mass_test = beta(3, 0);		//Determined from the trandformation

				//Get coefficients
				q1 = beta(0, 0);
				q2 = beta(1, 0);
				T alpha = atan2(q2, q1) * 180.0 / M_PI;
				*/
						
				//Computer centers of mass for both systems P, P'
				T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

				for (unsigned int j = 0; j < m; j++)
				{
					//Use only non singular points
					//if (W(j, j) != 0.0)
					{
						x_mass_test += nl_test[j]->getX();
						y_mass_test += nl_test[j]->getY();

						x_mass_reference += nl_projected_temp[j]->getX();
						y_mass_reference += nl_projected_temp[j]->getY();
					}
				}

				//Centers of mass
				x_mass_test = x_mass_test / m;
				y_mass_test = y_mass_test / m;
				x_mass_reference = x_mass_reference / m;
				y_mass_reference = y_mass_reference / m;

				//Compute scale using the least squares adjustment: h = inv (A'WA)A'WL, weighted Helmert transformation
				T sum_xy_1 = 0, sum_xy_2 = 0, sum_xx_yy = 0;
				for (unsigned int j = 0; j < m; j++)
				{
					sum_xy_1 = sum_xy_1 + (nl_test[j]->getX() - x_mass_test) * (nl_projected_temp[j]->getX() - x_mass_reference) +
						(nl_test[j]->getY() - y_mass_test) * (nl_projected_temp[j]->getY() - y_mass_reference);
					sum_xy_2 = sum_xy_2 + (nl_test[j]->getY() - y_mass_test) * (nl_projected_temp[j]->getX() - x_mass_reference) -
						(nl_test[j]->getX() - x_mass_test) * (nl_projected_temp[j]->getY() - y_mass_reference);
					sum_xx_yy = sum_xx_yy + (nl_projected_temp[j]->getX() - x_mass_reference) * (nl_projected_temp[j]->getX() - x_mass_reference) +
						(nl_projected_temp[j]->getY() - y_mass_reference) * (nl_projected_temp[j]->getY() - y_mass_reference);
				}

				//Transformation ratios
				q1 = sum_xy_1 / sum_xx_yy;
				q2 = sum_xy_2 / sum_xx_yy;

				//Rotation
				const T alpha = atan2(q2, q1) * 180.0 / M_PI;
				//std::cout << "rot = " << alpha;
				//std::cout << "   scale = " << sqrt(q1*q1 + q2*q2);
			
				//Outliers
				if (analysis_parameters.remove_outliers)
				{

					//Remove outliers
					Matrix <T> PR(m, 2), QR(m, 2), Eps(2 * m, 1);
					for (unsigned int j = 0; j < m; j++)
					{
						PR(j, 0) = (nl_test[j]->getX() - x_mass_test);
						PR(j, 1) = (nl_test[j]->getY() - y_mass_test);

						QR(j, 0) = (q1 * (nl_projected_temp[j]->getX() - x_mass_reference) - q2 * (nl_projected_temp[j]->getY() - y_mass_reference));
						QR(j, 1) = (q2 * (nl_projected_temp[j]->getX() - x_mass_reference) + q1 * (nl_projected_temp[j]->getY() - y_mass_reference));
					}

					//Remove  outliers
					T eps_init = 0, eps = 0;
					unsigned int iterations = 0;
					Outliers::findOutliersME(PR, QR, k, 1.0e-10, SimilarityScheme, me_function, 30, W, I, Eps, eps_init, eps, iterations);
				}

				//Compute coordinate differences (residuals): estimated - input
				Matrix <T> v(2 * m, 1);
				for (unsigned int j = 0; j < m; j++)
				{
					//Use only non singular points
					//if (W(j, j) != 0.0)
					{
						const T vx = (q1 * (nl_projected_temp[j]->getX() - x_mass_reference) - q2 * (nl_projected_temp[j]->getY() - y_mass_reference) - (nl_test[j]->getX() - x_mass_test));
						const T vy = (q2 * (nl_projected_temp[j]->getX() - x_mass_reference) + q1 * (nl_projected_temp[j]->getY() - y_mass_reference) - (nl_test[j]->getY() - y_mass_test));

						///Input is the best vector: first vertex of the simplex
						if ((m1 == 1) && (m2 > 1))
						{
							V(j, 0) = vx;
							V(j + m, 0) = vy;
						}

						//Input is a simplex
						else
						{
							v(j, 0) = vx;
							v(j + m, 0) = vy;
						}
					}
				}

				//Compute squares of residuals: only if input is a full simplex
				if (((m1 > 1) || (m2 == 1)))
				{
					V(i, 0) = MatrixOperations::norm(MatrixOperations::trans(v) * W * v);
				}

				//Set DX, DY, compute it from the first vertex of the simplex, which stores the current best solution
				if (i == 0)
				{
					//Compute DX, DY
					const T dx = x_mass_test - x_mass_reference * q1 + y_mass_reference * q2;
					const T dy = y_mass_test - x_mass_reference * q2 - y_mass_reference * q1;

					sample_res.setDx(dx);
					sample_res.setDy(dy);

					//Set rotation
					sample_res.setRotation(alpha);

					//Perform scaling
					R *= sqrt(q1 * q1 + q2 * q2);
					q1 = q1 * R;
					q2 = q2 * R;
				}

				res_evaluation++;
			}
		}

};


#endif
