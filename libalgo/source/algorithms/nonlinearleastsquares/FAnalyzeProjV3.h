// Description: Functor, compute matrix V of residuals for cartometric analysis
// Method: NLSP, M7 (including rotation)

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


#ifndef FAnalyzeProjV3_H
#define FAnalyzeProjV3_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV3
{
private:

	Container <Node3DCartesian <T> *> &nl_test;
	Container <Point3DGeographic <T> *> &pl_reference;
	Container <Node3DCartesianProjected <T> *> &nl_projected;
	typename TMeridiansList <T> ::Type &meridians;
	typename TParallelsList <T> ::Type &parallels;
	const Container <Face <T> *> &faces_test;
	Projection <T> *proj;
	T &x_mass_reference;
	T &y_mass_reference;
	const TAnalysisParameters <T> &analysis_parameters;
	const TProjectionAspect aspect;
	Sample <T> &sample_res;
	unsigned int & created_samples;
	std::ostream * output;

public:

	FAnalyzeProjV3(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, Container <Node3DCartesianProjected <T> *> &nl_projected_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, T &x_mass_reference_, T &y_mass_reference_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), nl_projected(nl_projected_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), x_mass_reference(x_mass_reference_), y_mass_reference(y_mass_reference_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), output(output_) {}

	void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
	{

		//Compute parameters of the V matrix: residuals
		const unsigned int m = nl_test.size();

		//Get lat0 min and lat0 max
		const T lat0_min = proj->getLat0Interval().min_val;
		const T lat0_max = proj->getLat0Interval().max_val;

		//Normal aspect: lat0, lon0
		if (aspect == NormalAspect)
		{
			//Correct R, lat0, lon0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Set lat0 inside the interval
			if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

			//Set lon0
			if (X(4, 0) < MIN_LON)  X(4, 0) = MAX_LON + fmod(X(4, 0), MIN_LON);
			else if (X(4, 0) > MAX_LON)  X(4, 0) = MIN_LON + fmod(X(4, 0), MAX_LON);

			//Set c inside the interval
			//if (fabs(X(5, 0)) < 0) X(5, 0) = -X(5, 0);

			//if (fabs(X(5, 0)) > MAX_C) X(5, 0) = 2 * MAX_C - X(5, 0);

			//Subtract period
			if (fabs(X(6, 0)) > MAX_LON) X(6, 0) = fmod(X(6, 0), 180);
		}

		//Transverse aspect: lonp, lat0
		else  if (aspect == TransverseAspect)

		{
			//Correct R, lonp, lat0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Subtract period
			if (X(2, 0) < MIN_LON)  X(2, 0) = MAX_LON + fmod(X(2, 0), MIN_LON);
			else if (X(2, 0) > MAX_LON)  X(2, 0) = MIN_LON + fmod(X(2, 0), MAX_LON);

			//Set lat0 inside the interval
			if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

			//Set lon0
			X(4, 0) = 0;

			//Set c inside the interval
			//if (fabs(X(5, 0)) < 0) X(5, 0) = -X(5, 0);

			//if (fabs(X(5, 0)) > MAX_C) X(5, 0) = 2 * MAX_C - X(5, 0);

			//Subtract period
			if (fabs(X(6, 0)) > MAX_LON) X(6, 0) = fmod(X(6, 0), 180);
		}

		//Oblique aspect: latp, lonp, lat0
		else if (aspect == ObliqueAspect)
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
			if (fabs(fabs(X(1, 0)) - MAX_LAT) < 1.0)
			{
				//X(1, 0) = 90.0;
				X(2, 0) = 0.0;
			}

			//Set lon0
			X(4, 0) = 0;

			//Set c inside the interval
			//if (fabs(X(5, 0)) < 0) X(5, 0) = -X(5, 0);

			//if (fabs(X(5, 0)) > MAX_C) X(5, 0) = 2 * MAX_C - X(5, 0);

			//Subtract period
			if (fabs(X(6, 0)) > MAX_LON) X(6, 0) = fmod(X(6, 0), 180);
				
		}

		//Set properties to the projection: ommit estimated radius, additional constants dx, dy
		// They will be estimated again using the transformation
		Point3DGeographic <T> cart_pole(X(1, 0), X(2, 0));
		proj->setR(X(0, 0));
		proj->setCartPole(cart_pole);
		proj->setLat0(X(3, 0));
		proj->setLon0(X(4, 0));
		proj->setDx(0.0);
		proj->setDy(0.0);
		proj->setC(X(5, 0));

		//Get aplha
		const T alpha = X(6, 0);

		//Clear points from the previous iteration
		nl_projected.clear();

		//Compute coordinate differences (residuals): items of V matrix
		for (unsigned int i = 0; i < m; i++)
		{
			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Reduce lon
			const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(4, 0));

			T lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0;

			try
			{
				//Convert geographic point to oblique aspect
				lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0));
				lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0), trans_lon_dir);

				for (unsigned int j = 0; j < 3; j++)
				{
					try
					{
						//Compute x, y coordinates
						x = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), 0.0, X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
						y = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), 0.0, X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);

						//x = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
						//y = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
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
			}

			//Create new cartographic point
			Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

			//Add point to the list
			nl_projected.push_back(n_projected);
		}

		//Computer centers of mass for both systems P, P'
		unsigned int n_points = 0;
		T x_mass_test = 0.0, y_mass_test = 0.0;
		x_mass_reference = 0.0, y_mass_reference = 0.0;

		for (unsigned int i = 0; i < m; i++)
		{
			//Use only non singular points
			if (W(i, i) != 0.0)
			{
				x_mass_test += nl_test[i]->getX();
				y_mass_test += nl_test[i]->getY();

				x_mass_reference += nl_projected[i]->getX();
				y_mass_reference += nl_projected[i]->getY();

				n_points++;
			}
		}

		//Centers of mass
		x_mass_test = x_mass_test / n_points;
		y_mass_test = y_mass_test / n_points;
		x_mass_reference = x_mass_reference / n_points;
		y_mass_reference = y_mass_reference / n_points;

		//Compute coordinate differences (residuals): estimated - input
		for (unsigned int i = 0; i < m; i++)
		{
			//Use only non singular points
			if (W(i, i) != 0.0)
			{
				V(i, 0) = ((nl_projected[i]->getX() - x_mass_reference) * cos(alpha * M_PI / 180) - (nl_projected[i]->getY() - y_mass_reference) * sin(alpha * M_PI / 180) - (nl_test[i]->getX() - x_mass_test));
				V(i + m, 0) = ((nl_projected[i]->getX() - x_mass_reference) * sin(alpha * M_PI / 180) + (nl_projected[i]->getY() - y_mass_reference) * cos(alpha * M_PI / 180) - (nl_test[i]->getY() - y_mass_test));
			}
		}

		//Compute DX, DY
		T dx = x_mass_test - x_mass_reference * cos(alpha * M_PI / 180) + y_mass_reference * sin(alpha * M_PI / 180);
		T dy = y_mass_test - x_mass_reference * sin(alpha * M_PI / 180) - y_mass_reference * cos(alpha * M_PI / 180);

		//Store shifts and rotation
		sample_res.setDx(dx);
		sample_res.setDy(dy);
		sample_res.setRotation(alpha);

	}
};


#endif
