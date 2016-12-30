// Description: Functor, compute matrix V of squares of residuals for cartometric analysis
// Method: Differential evolution, determined shifts

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


#ifndef FAnalyzeProjVDE_H
#define FAnalyzeProjVDE_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjVDE
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

	unsigned int &iter;

public:

	FAnalyzeProjVDE(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &iter_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), iter(iter_), output(output_) {}

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
			if (X(0, 3) < lat0_min || X(0, 3) > lat0_max) X(0, 3) = 0.5 * (lat0_min + lat0_max);

			//Set lon0
			if (X(0, 4) < MIN_LON)  X(0, 4) = MAX_LON + fmod(X(0, 4), MIN_LON);
			else if (X(0, 4) > MAX_LON)  X(0, 4) = MIN_LON + fmod(X(0, 4), MAX_LON);
		}

		//Transverse aspect: lonp, lat0
		else  if (aspect == TransverseAspect)

		{
			//Correct R, lonp, lat0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Subtract period
			if(X(0, 2) < MIN_LON)  X(0, 2) = MAX_LON + fmod(X(0, 2), MIN_LON);
			else if (X(0, 2) > MAX_LON)  X(0, 2) = MIN_LON + fmod(X(0, 2), MAX_LON);

			//Set lat0 inside the interval
			if (X(0, 3) < lat0_min || X(0, 3) > lat0_max) X(0, 3) = 0.5 * (lat0_min + lat0_max);

			//Set lon0
			X(0, 4) = 0;
		}

		//Oblique aspect: latp, lonp, lat0
		else if (aspect == ObliqueAspect)
		{
			//Correct R, latp, lonp, lat0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Subtract period
			if (X(0, 1) < MIN_LAT)  X(0, 1) = MIN_LAT - fmod(X(0, 1), MIN_LAT);
			else if (X(0, 1) > MAX_LAT)  X(0, 1) = MAX_LAT - fmod(X(0, 1), MAX_LAT);

			if (X(0, 2) < MIN_LON)  X(0, 2) = MAX_LON + fmod(X(0, 2), MIN_LON);
			else if (X(0, 2) > MAX_LON)  X(0, 2) = MIN_LON + fmod(X(0, 2), MAX_LON);

			//Set lat0 inside the interval
			if (X(0, 3) < lat0_min || X(0, 3) > lat0_max) X(0, 3) = 0.5 * (lat0_min + lat0_max);

			//Set lonp to zero, if latp = 90
			if (fabs(X(0, 1) - MAX_LAT) < 1.0)
			{
				//X(0, 1) = 90.0;
				//X(0, 2) = 0.0;
			}

			//Set lon0
			X(0, 4) = 0;
		}

		//Set properties to the projection: ommit estimated radius, additional constants dx, dy
		// They will be estimated again using the transformation
		Point3DGeographic <T> cart_pole(X(0, 1), X(0, 2));
		proj->setR(X(0, 0));
		proj->setCartPole(cart_pole);
		proj->setLat0(X(0, 3));
		proj->setLat1(X(0, 3));
		proj->setLat2(X(0, 7));
		proj->setLon0(X(0, 4));
		proj->setDx(X(0, 5));
		proj->setDy(X(0, 6));
		proj->setC(X(0, 7));

		//Compute new coordinates
		for (unsigned int i = 0; i < m; i++)
		{
			T x_new, y_new;

			//Get type of the direction
			TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

			//Reduce lon
			const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(0, 4));

			try
			{
				//Convert geographic point to oblique aspect
				const T lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2));
				const T lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2), trans_lon_dir);

				for (unsigned int j = 0; j < 3; j++)
				{
					try
					{
						//Compute new coordinates: add shifts
						Y(i, 0) = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 5), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false);
						Y(i + m, 0) = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 6), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false);

						//Y(i, 0) = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false) + X(0, 5);
						//Y(i + m, 0) = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false) + X(0, 6);
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
			catch (Exception & error)
			{
				x_new = 1.0;
				y_new = 0.0;
			}

			//Compute coordinate differences (residuals): estimated - input
			V(i, 0) = (Y(i, 0) - nl_test[i]->getX());
			V(i + m, 0) = (Y(i + m, 0) - nl_test[i]->getY());
		}

		iter++;
	}

};


#endif