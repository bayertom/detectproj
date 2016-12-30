// Description: Functor, compute matrix V of residuals for cartometric analysis
// Method: NLSP, M6 (without rotation)

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


#ifndef FAnalyzeProjV2_HPP
#define FAnalyzeProjV2_HPP
/*
#include "FAnalyzeProjV2B.h"

template <typename T>
void FAnalyzeProjV2<T>:: operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis)
{

	//Compute parameters of the V matrix: residuals
	//Get lat0 min and lat0 max
	const T lat0_min = proj->getLat0Interval().min_val;
	const T lat0_max = proj->getLat0Interval().max_val;

	//Normal aspect: lat0, lon0
	if (aspect == NormalAspect)
	{
		//Correct R, lat0, lon0
		if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

		//Set lat0 inside the interval
		//if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

		//Set lon0
		if (X(4, 0) < MIN_LON) X(4, 0) = X(4, 0) + 360;
		else if (X(4, 0) > MAX_LON) X(4, 0) = X(4, 0) - 360;

		//Set to interval
		if (fabs(X(3, 0)) > MAX_LAT)  X(3, 0) = fmod(X(3, 0), 90);
		//if (X(3, 0) < lat0_min) X(3, 0) = lat0_min;

		//if (X(3, 0) > lat0_max) X(3, 0) = lat0_max;
	}

	//Transverse aspect: lonp, lat0
	else  if (aspect == TransverseAspect)

	{
		//Correct R, lonp, lat0
		if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

		//Subtract period
		//if (fabs(X(2, 0)) > MAX_LON) X(2, 0) = fmod(X(2, 0), 180);

		if (X(2, 0) < MIN_LON)  X(2, 0) = MIN_LON - fmod(X(2, 0), MIN_LON);
		else if (X(2, 0) > MAX_LON)  X(2, 0) = MAX_LON - fmod(X(2, 0), MAX_LON);

		//Set lat0 inside the interval
		if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

		//Set lon0
		X(4, 0) = 0;

	}

	//Oblique aspect: latp, lonp, lat0
	else if (aspect == ObliqueAspect)
	{
		//Correct R, latp, lonp, lat0
		if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

		//Subtract period
		//if (fabs(X(1, 0)) > MAX_LAT)  X(1, 0) = fmod(X(1, 0), 90);

		//if (fabs(X(2, 0)) > MAX_LON)  X(2, 0) = fmod(X(2, 0), 180);

		//if (fabs(X(3, 0)) > MAX_LAT)  X(3, 0) = fmod(X(3, 0), 90);

		if (X(1, 0) < MIN_LAT)  X(1, 0) = MIN_LAT - fmod(X(1, 0), MIN_LAT);
		else if (X(1, 0) > MAX_LAT)  X(1, 0) = MAX_LAT - fmod(X(1, 0), MAX_LAT);

		if (X(2, 0) < MIN_LON)  X(2, 0) = MIN_LON - fmod(X(2, 0), MIN_LON);
		else if (X(2, 0) > MAX_LON)  X(2, 0) = MAX_LON - fmod(X(2, 0), MAX_LON);

		//Set lat0 inside the interval
		if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

		//Set lonp to zero, if latp = 90
		//if (fabs(X(1, 0) - MAX_LAT) < 1.0)  X(2, 0) = 0.0;

		//Set lon0
		X(4, 0) = 0;

		//Set lonp to zero, if latp = 90
		if (fabs(X(1, 0) - MAX_LAT) < 5.0)
		{
			//X(1, 0) = 90.0;
			//X(2, 0) = 0.0;
		}
	}

	//Compute residuals
	evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test,proj, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);

	//Evaluate bisection, if the coordinate function is the linear function of lon (cylindtrical projections)
	//if (enable_bisection)
	{
		unsigned short iterations = 0;
		const unsigned short max_iterations = 20, m = X.rows(), n = V.rows();

		T  xmin, fmin, lon0;
		const T eps = 0.1, max_diff = 0.01;

		//Set limits for bisection: copy of the currennt solution + [m, n]
		/*
		Matrix <T> A(m + 2, 1);
		A.submat(X, 0, 0);
		A(6, 0) = m; A(7, 0) = n;
		Matrix <T> B = A;
		A(4, 0) = -180; B(4, 0) = 180;
		*/
		Matrix <T> A = X, B = X;
		A(4, 0) = -180; B(4, 0) = 180;
		X.print();

		//Run bisection
		//Bisection::bisection(evaluateResidualsW <T>, A, B, eps, max_diff, xmin, fmin, iterations, max_iterations);

		Bisection::bisection(FAnalyzeProjV2B <T>(nl_test, pl_reference, meridians, parallels, faces_test, proj,
			analysis_parameters, aspect, sample_res, created_samples, res_evaluation, HuberFunction, k, I, output), A, B, eps, max_diff, xmin, fmin, iterations, max_iterations);

		//Assign determined bisection
		X(4, 0) = xmin;
		//X(4, 0) = -120;

	}
}



template <typename T>
inline void evaluateResiduals(const Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples, unsigned short &res_evaluation, const TMEstimatorsWeightFunction &me_function, const T k, Matrix <unsigned int> &I, std::ostream * output)
{
	//Evaluate residuals
	const unsigned int m = nl_test.size();

	//Set properties to the projection: ommit estimated radius, additional constants dx, dy
	// They will be estimated again using the transformation
	Point3DGeographic <T> cart_pole(X(1, 0), X(2, 0));
	proj->setCartPole(cart_pole);
	proj->setR(X(0, 0));
	proj->setLat0(X(3, 0));
	proj->setLon0(X(4, 0));
	proj->setDx(0.0);
	proj->setDy(0.0);
	proj->setC(0.0);

	W = eye(2 * m, 2 * m, 1.0);

	//Compute coordinate differences (residuals): items of V matrix
	Container <Node3DCartesianProjected <T> *> nl_projected_temp;

	for (unsigned int i = 0; i < m; i++)
	{
		//Get type of the direction
		TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

		//Reduce lon
		const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(4, 0));

		T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

		try
		{
			//Convert geographic point to oblique position: use a normal direction of converted longitude
			lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0));
			lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, lat_trans, X(1, 0), X(2, 0), trans_lon_dir);

			//Compute x, y coordinates
			x = ArithmeticParser::parseEq(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
			y = ArithmeticParser::parseEq(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
		}

		catch (Error &error)
		{
			//Disable point from analysis: set weight to zero
			W(i, i) = 0; W(i + m, i + m) = 0;
		}

		//Create new cartographic point
		Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

		//Add point to the list
		nl_projected_temp.push_back(n_projected);
	}

	//nl_projected_temp.print(output);

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

	//nl_projected_temp.print();

	//Outliers
	if (analysis_parameters.remove_outliers)
	{
		//Remove outliers
		Matrix <T> PR(m, 2), QR(m, 2), Eps(2 * m, 1);
		for (unsigned int i = 0; i < m; i++)
		{
			PR(i, 0) = (nl_test[i]->getX() - x_mass_test);
			PR(i, 1) = (nl_test[i]->getY() - y_mass_test);

			QR(i, 0) = (nl_projected_temp[i]->getX() - x_mass_reference);
			QR(i, 1) = (nl_projected_temp[i]->getY() - y_mass_reference);
		}

		//Remove  outliers
		T eps_init = 0, eps = 0;
		unsigned int iterations = 0;
		Outliers::findOutliersME(PR, QR, k, 1.0e-10, ScaleShiftsScheme, me_function, 30, W, I, Eps, eps_init, eps, iterations);
	}


	//Compute coordinate differences (residuals): estimated - input
	for (unsigned int i = 0; i < m; i++)
	{
		//Use only non singular points
		//if (W(i, i) != 0.0)
		{
			V(i, 0) = ((nl_projected_temp[i]->getX() - x_mass_reference) - (nl_test[i]->getX() - x_mass_test));
			V(i + m, 0) = ((nl_projected_temp[i]->getY() - y_mass_reference) - (nl_test[i]->getY() - y_mass_test));
		}
	}

	//Set DX, DY
	sample_res.setDx(x_mass_test - x_mass_reference);
	sample_res.setDy(y_mass_test - y_mass_reference);

	res_evaluation++;
}
*/

#endif
