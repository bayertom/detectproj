// Description: Functor, compute matrix V of residuals for cartometric analysis
// Method: NLSP, M7 (5 determined parameters, without rotation, without radius)

// Copyright (c) 2010 - 2015
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


#ifndef FAnalyzeProjV4_H
#define FAnalyzeProjV4_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"

#include "libalgo/source/algorithms/geneticalgorithms/FAnalyzeProjV4DEL.h"

//Forwards declarations
template <typename T>
class Projection;

template <typename T>
class FAnalyzeProjV4DEL;


// Functor, compute matrix V of residuals for cartometric analysis
template <typename T>
class FAnalyzeProjV4
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
	bool &enable_additional_lon0_analysis;
	std::ostream * output;

public:

	FAnalyzeProjV4(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, T &R_est_, T &q1_, T &q2_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_,
		unsigned int & created_samples_, unsigned int &res_evaluation_, const TMEstimatorsWeightFunction &me_function_, const T k_, Matrix <unsigned int> &I_, bool & enable_additional_lon0_analysis_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R(R_est_), q1(q1_), q2(q2_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), res_evaluation(res_evaluation_), me_function(me_function_), k(k_), I(I_), enable_additional_lon0_analysis(enable_additional_lon0_analysis_), output(output_) {}


	void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
	{

		//Compute parameters of the V matrix: residuals
		const unsigned int m = nl_test.size();
		const unsigned int n = X.rows();

		//Get lat0 min and lat0 max
		const T lat0_min = proj->getLat0Interval().min_val;
		const T lat0_max = proj->getLat0Interval().max_val;

		//Normal aspect: lat0, lon0
		if (aspect == NormalAspect)
		{
			//Set lat0 inside the interval
			//if (X(2, 0) < lat0_min || X(2, 0) > lat0_max) X(2, 0) = 0.5 * (lat0_min + lat0_max);

			//Set lat0 to interval
			if (fabs(X(2, 0)) > MAX_LAT)  X(2, 0) = fmod(X(2, 0), 90);

			//Set lon0
			if (X(3, 0) < MIN_LON)  X(3, 0) = MAX_LON + fmod(X(3, 0), MIN_LON);
			else if (X(3, 0) > MAX_LON)  X(3, 0) = MIN_LON + fmod(X(3, 0), MAX_LON);

			//Set lon0 to interval
			//if (fabs(X(3, 0)) > MAX_LON) 	X(3, 0) = fmod(X(3, 0), 90);
		}

		//Transverse aspect: lonp, lat0
		else  if (aspect == TransverseAspect)

		{
			//Subtract period
			if (X(1, 0) < MIN_LON)  X(1, 0) = MAX_LON + fmod(X(1, 0), MIN_LON);
			else if (X(1, 0) > MAX_LON)  X(1, 0) = MIN_LON + fmod(X(1, 0), MAX_LON);

			//Set lat0 inside the interval
			if (X(2, 0) < lat0_min || X(2, 0) > lat0_max) X(2, 0) = 0.5 * (lat0_min + lat0_max);

			//Set lon0
			X(3, 0) = 0;
		}

		//Oblique aspect: latp, lonp, lat0
		else if (aspect == ObliqueAspect)
		{
			
			//Subtract period
			if (X(0, 0) < MIN_LAT)  X(0, 0) = MIN_LAT - fmod(X(0, 0), MIN_LAT);
			else if (X(0, 0) > MAX_LAT)  X(0, 0) = MAX_LAT - fmod(X(0, 0), MAX_LAT);

			if (X(1, 0) < MIN_LON)  X(1, 0) = MAX_LON + fmod(X(1, 0), MIN_LON);
			else if (X(1, 0) > MAX_LON)  X(1, 0) = MIN_LON + fmod(X(1, 0), MAX_LON);
			
			//Set lat0 inside the interval
			if (X(2, 0) < lat0_min) X(2, 0) = lat0_min + fabs(fmod(X(2, 0), MIN_LAT));
			if (X(2, 0) > lat0_max) X(2, 0) = lat0_max - fabs(fmod(X(2, 0), MAX_LAT));
			
			//Set lonp to zero, if latp = 90
			if (fabs(X(0, 0) - MAX_LAT) < 5)
			{
				//X(0, 0) = 90.0;
				//X(1, 0) = 0.0;
			}

			//Set lon0
			//X(3, 0) = 0;
			if (X(3, 0) < MIN_LON)  X(3, 0) = MAX_LON + fmod(X(3, 0), MIN_LON);
			else if (X(3, 0) > MAX_LON)  X(3, 0) = MIN_LON + fmod(X(3, 0), MAX_LON);
		}

		// Set properties of the projection
		Point3DGeographic <T> cart_pole(X(0, 0), X(1, 0));
		proj->setR(R);
		proj->setCartPole(cart_pole);
		proj->setLat0(X(2, 0));
		proj->setLat1(X(2, 0));
		proj->setLat2(X(4, 0));
		proj->setLon0(X(3, 0));
		proj->setDx(0.0);
		proj->setDy(0.0);
		proj->setC(X(4, 0));

		//Additional determination of lon0 using the differential evolution
		T min_cost = 0, min_cost_lon0 = 0;
		Matrix <T> XL(1, 1);

		if ( (enable_additional_lon0_analysis) && (aspect == NormalAspect) && (compute_analysis) )
		{
			//Create matrices
			const unsigned int population = 4 * n, max_gen = 5;
			unsigned int res_evaluations_lon0 = 0,  iterations_lon0 = 0;
			const T CR = 0.8, eps = 0.1;
			T res_aver = 0, res_max = 0;

			Matrix <unsigned int> IX(2 * m, 1);
			Matrix <T>   XLMIN(1, 1), XLMAX(1, 1), XAVER(1, 1), YL(2 * m, 1), WL(2 * m, 2 * m, 0.0, 1), VL(2 * m, 1), F(1, 1);
			
			//Initialize bounds
			F(0, 0) = 0.5;
			XLMIN(0, 0) = -MAX_LON; XLMAX(0, 0) = MAX_LON;
			XL(0, 0) = X(3, 0);

			//Compute DE solution
			min_cost_lon0 = DifferentialEvolution::diffEvolution(FAnalyzeProjV4DEL <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res,
				created_samples, res_evaluations_lon0, me_function, k, IX, output), population, eps, max_gen, F, CR, DERandBest1Strategy, MFDE, WL, XL, YL, VL, XLMIN, XLMAX, XAVER, res_aver, res_max, iterations_lon0, true);
		}

		//Compute residuals for the NLSP solution
		evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);
		min_cost = norm(trans(V)*W*V);

		//Compare NLSP and DE solutions, minimum is flat
		if ((enable_additional_lon0_analysis) && (aspect == NormalAspect) && (fabs(min_cost - min_cost_lon0) > 1.0e-5) && (compute_analysis))
		{
			X(3, 0) = XL(0, 0);
		}

		//Assign lon0 to proj
		proj->setLon0(X(3, 0));
	}


	void evaluateResiduals(const Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
		const Container <Face <T> *> &faces_test, Projection <T> *proj, T &R, T &q1, T &q2, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples, unsigned int &res_evaluation, const TMEstimatorsWeightFunction &me_function, const T k, Matrix <unsigned int> &I, std::ostream * output)
	{
		//Evaluate residuals
		const unsigned int m = nl_test.size();

		W = eye(2 * m, 2 * m, 1.0);

		//Compute coordinate differences (residuals): items of V matrix
		Container <Node3DCartesianProjected <T> *> nl_projected_temp;

		for (unsigned int i = 0; i < m; i++)
		{
			T lat = pl_reference[i]->getLat();

			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Reduce lon
			T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(3, 0));

			T lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0;

			try
			{
				//Convert geographic point to oblique aspect
				//lat_trans = CartTransformation::latToLatTrans(lat, lon_red, X(0, 0), X(1, 0));
				//lon_trans = CartTransformation::lonToLonTrans(lat, lon_red, lat_trans, X(0, 0), X(1, 0), trans_lon_dir);

				//More safe lat/lon conversion
				for (unsigned int j = 0; j < 3; j++)
				{
					try
					{
						//Convert geographic point to oblique aspect
						lat_trans = CartTransformation::latToLatTrans(lat , lon_red, X(0, 0), X(1, 0));
						lon_trans = CartTransformation::lonToLonTrans(lat, lon_red, X(0, 0), X(1, 0), trans_lon_dir);
					}

					//2 attempt to avoid the singularity
					catch (Exception &error)
					{
						//Move in latitude direction
						if (j == 0)
						{
							if (lat == MAX_LAT || lat == X(0, 0))
								lat -= GRATICULE_ANGLE_SHIFT;
							else 
								lat += GRATICULE_ANGLE_SHIFT;
						}

						//Move in longitude direction
						else if (j == 1)
						{
							if (lon_red == MAX_LON || lon_red == X(1, 0))
								lon_red -= GRATICULE_ANGLE_SHIFT;
							else 
								lon_red += GRATICULE_ANGLE_SHIFT;
						}

						//Neither first nor the second shhifts do not bring improvement
						else if (j == 2)
						{
							throw;
						}
					}
				}
				
				//More safe coordinate computation
				for (unsigned int j = 0; j < 3; j++)
				{
					try
					{
						//Compute x, y coordinates
						x = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(4, 0), X(2, 0), X(2, 0), X(4, 0), false);
						y = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(4, 0), X(2, 0), X(2, 0), X(4, 0), false);

						//x = CartTransformation::latLonToX(proj->getXEquat(), proj->getFThetaEquat(), proj->getTheta0Equat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(4, 0), X(2, 0), proj->getLat1(), proj->getLat2(), false);
						//y = CartTransformation::latLonToY(proj->getYEquat(), proj->getFThetaEquat(), proj->getTheta0Equat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), 0.0, X(4, 0), X(2, 0), proj->getLat1(), proj->getLat2(), false);


						//x = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), X(4, 0), X(2, 0), proj->getLat1(), proj->getLat2(), false);
						//y = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, R, proj->getA(), proj->getB(), X(4, 0), X(2, 0), proj->getLat1(), proj->getLat2(), false);
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

			catch (Exception &error)
			{
				//Disable point from analysis: set weight to zero
				W(i, i) = 0; W(i + m, i + m) = 0;
				std::cout << "Error";
				error.printException();
			}

			//Create new cartographic point
			Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

			//Add point to the list
			nl_projected_temp.push_back(n_projected);
		}

		//Computer centers of mass for both systems P, P'
		T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

		for (unsigned int i = 0; i < m; i++)
		{
			//Use only non singular points
			//if (W(i, i) != 0.0)
			{
				x_mass_test += nl_test[i]->getX();
				y_mass_test += nl_test[i]->getY();

				x_mass_reference += nl_projected_temp[i]->getX();
				y_mass_reference += nl_projected_temp[i]->getY();
			}
		}

		//Centers of mass
		x_mass_test = x_mass_test / m;
		y_mass_test = y_mass_test / m;
		x_mass_reference = x_mass_reference / m;
		y_mass_reference = y_mass_reference / m;

		//Compute scale using the least squares adjustment: h = inv (A'WA)A'WL, weighted Helmert transformation
		T sum_xy_1 = 0, sum_xy_2 = 0, sum_xx_yy = 0;
		for (unsigned int i = 0; i < m; i++)
		{
			//if (W(i, i) != 0.0)
			{
				sum_xy_1 = sum_xy_1 + (nl_test[i]->getX() - x_mass_test) * (nl_projected_temp[i]->getX() - x_mass_reference) +
					(nl_test[i]->getY() - y_mass_test) * (nl_projected_temp[i]->getY() - y_mass_reference);
				sum_xy_2 = sum_xy_2 + (nl_test[i]->getY() - y_mass_test) * (nl_projected_temp[i]->getX() - x_mass_reference) -
					(nl_test[i]->getX() - x_mass_test) * (nl_projected_temp[i]->getY() - y_mass_reference);
				sum_xx_yy = sum_xx_yy + (nl_projected_temp[i]->getX() - x_mass_reference) * (nl_projected_temp[i]->getX() - x_mass_reference) +
					(nl_projected_temp[i]->getY() - y_mass_reference) * (nl_projected_temp[i]->getY() - y_mass_reference);
			}
		}

		//Transformation ratios
		q1 = sum_xy_1 / sum_xx_yy;
		q2 = sum_xy_2 / sum_xx_yy;

		//Rotation
		const T alpha = atan2(q2, q1) * 180.0 / M_PI;

		//Outliers
		if (analysis_parameters.remove_outliers)
		{

			//Remove outliers
			Matrix <T> PR(m, 2), QR(m, 2), Eps(2 * m, 1);
			for (unsigned int i = 0; i < m; i++)
			{
				PR(i, 0) = (nl_test[i]->getX() - x_mass_test);
				PR(i, 1) = (nl_test[i]->getY() - y_mass_test);

				QR(i, 0) = (q1 * (nl_projected_temp[i]->getX() - x_mass_reference) - q2 * (nl_projected_temp[i]->getY() - y_mass_reference));
				QR(i, 1) = (q2 * (nl_projected_temp[i]->getX() - x_mass_reference) + q1 * (nl_projected_temp[i]->getY() - y_mass_reference));
			}

			//Remove  outliers
			T eps_init = 0, eps = 0;
			unsigned int iterations = 0;
			Outliers::findOutliersME(PR, QR, k, 1.0e-10, SimilarityScheme, me_function, 30, W, I, Eps, eps_init, eps, iterations);
		}

		//Compute coordinate differences (residuals): estimated - input
		for (unsigned int i = 0; i < m; i++)
		{
			//Use only non singular points
			//if (W(i, i) != 0.0)
			{
				V(i, 0) = (q1 * (nl_projected_temp[i]->getX() - x_mass_reference) - q2 * (nl_projected_temp[i]->getY() - y_mass_reference) - (nl_test[i]->getX() - x_mass_test));
				V(i + m, 0) = (q2 * (nl_projected_temp[i]->getX() - x_mass_reference) + q1 * (nl_projected_temp[i]->getY() - y_mass_reference) - (nl_test[i]->getY() - y_mass_test));
			}
		}

		//Compute DX, DY
		const T dx = x_mass_test - x_mass_reference * q1 + y_mass_reference * q2;
		const T dy = y_mass_test - x_mass_reference * q2 - y_mass_reference * q1;

		sample_res.setDx(dx);
		sample_res.setDy(dy);

		//V.print();
		//std::cout << q1 << "  " << q2;

		//Set rotation
		sample_res.setRotation(alpha);

		//Perform scaling
		R *= sqrt(q1 * q1 + q2 * q2);
		
		//R *= sqrt(q1 * q1 + q2 * q2);

		//q1 /= 10;
		//q2 /= 10;
		//R = 6348;
		//std::cout << "R= " << R << "  \n";
		//*output << "R= " << R << "  \n";
		//*output << "alpha= " << alpha << "  \n";

		res_evaluation++;
	}
};

#endif
