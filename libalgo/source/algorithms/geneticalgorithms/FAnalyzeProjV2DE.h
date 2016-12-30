// Description: Functor, compute matrix V of squares of residuals for cartometric analysis
// Method: Differential evolution, no shifts are determined

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


#ifndef FAnalyzeProjV2DE_H
#define FAnalyzeProjV2DE_H



#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"

//Forward declaration
template <typename T>
class Projection;


// Functor, compute matrix V of squares of residuals for cartometric analysis
template <typename T>
class FAnalyzeProjV2DE
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
		unsigned int &res_evaluations;
		TMEstimatorsWeightFunction me_function;
		T k;
		Matrix <unsigned int> &I;
                std::ostream * output;


        public:

                FAnalyzeProjV2DE ( Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
			const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &res_evaluations_, const TMEstimatorsWeightFunction &me_function_, const T k_, Matrix <unsigned int> &I_, std::ostream * output_)
                        : nl_test ( nl_test_ ), pl_reference ( pl_reference_ ), meridians ( meridians_ ), parallels ( parallels_ ), faces_test ( faces_test_ ),  proj ( proj_ ), analysis_parameters ( analysis_parameters_ ), aspect ( aspect_ ), sample_res ( sample_res_ ),
			created_samples(created_samples_), res_evaluations(res_evaluations_), me_function(me_function_), k(k_), I(I_), output(output_) {}

                void operator () ( Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true )
                {

                        //Compute squares of residuals V = v' * W *v
                        const unsigned int m = nl_test.size();
                        const unsigned int m1 = X.rows();

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
				if (X(0, 2) < MIN_LON)  X(0, 2) = MAX_LON - fmod(X(0, 2), MIN_LON);
				else if (X(0, 2) > MAX_LON)  X(0, 2) = MIN_LON - fmod(X(0, 2), MAX_LON);

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
			proj->setCartPole(cart_pole);
			proj->setLat0(X(0, 3));
			proj->setLat1(X(0, 3));
			proj->setLat2(X(0, 5));
			proj->setLon0(X(0, 4));
			proj->setDx(0.0);
			proj->setDy(0.0);
			proj->setC(X(0, 5));
			proj->setR(X(0, 0));

			//Compute residuals
			evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj,analysis_parameters, aspect, sample_res, created_samples, res_evaluations, me_function, k, I, output);
                }
                
};



template <typename T>
inline void evaluateResiduals(const Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples, unsigned int &res_evaluations, const TMEstimatorsWeightFunction &me_function, const T k, Matrix <unsigned int> &I, std::ostream * output)
{
	//Evaluate residuals
	const unsigned int m = nl_test.size();

	W = eye(2 * m, 2 * m, 1.0);

	//Compute coordinate differences (residuals): items of V matrix
	Container <Node3DCartesianProjected <T> *> nl_projected_temp;

	for (unsigned int i = 0; i < m; i++)
	{
		//Get type of the direction
		TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

		//Reduce lon
		const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(0, 4));

		//Convert geographic point to oblique aspect
		T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

		try
		{
			//Convert geographic point to oblique aspect
			lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2));
			lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2), trans_lon_dir);

			for (unsigned int j = 0; j < 3; j++)
			{
				try
				{
					//Compute x, y coordinates
					x = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), 0.0, X(0, 5), X(0, 3), X(0, 3), X(0, 5), false);
					y = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), 0.0, X(0, 5), X(0, 3), X(0, 3), X(0, 5), false);

					//x = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 5), X(0, 3), proj->getLat1(), proj->getLat2(), false);
					//y = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 5), X(0, 3), proj->getLat1(), proj->getLat2(), false);
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
		nl_projected_temp.push_back(n_projected);
	}

	//Computer centers of mass for both systems P, P'
	T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

	for (unsigned int i = 0; i < m; i++)
	{
		//Use only non singular points
		//if ( W ( i, i ) != 0.0 )
		{
			x_mass_test += nl_test[i]->getX();
			y_mass_test += nl_test[i]->getY();

			x_mass_reference += nl_projected_temp[i]->getX();
			y_mass_reference += nl_projected_temp[i]->getY();
		}
	}

	x_mass_test = x_mass_test / m;
	y_mass_test = y_mass_test / m;
	x_mass_reference = x_mass_reference / m;
	y_mass_reference = y_mass_reference / m;

	//Outliers
	if (analysis_parameters.remove_outliers)
	{
		//Remove outliers
		Matrix <T> PR(m, 2), QR(m, 2), Eps(2 * m, 1);
		for (unsigned int i = 0; i < m; i++)
		{
			PR(i, 0) = nl_test[i]->getX() - x_mass_test;
			PR(i, 1) = nl_test[i]->getY() - y_mass_test;

			QR(i, 0) = nl_projected_temp[i]->getX() - x_mass_reference;
			QR(i, 1) = nl_projected_temp[i]->getY() - y_mass_reference;
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
		//if ( W ( i, i ) != 0.0 )
		{
			V(i, 0) = ((nl_projected_temp[i]->getX() - x_mass_reference) - (nl_test[i]->getX() - x_mass_test));
			V(i + m, 0) = ((nl_projected_temp[i]->getY() - y_mass_reference) - (nl_test[i]->getY() - y_mass_test));
		}
	}

	//Set DX, DY
	sample_res.setDx(x_mass_test - x_mass_reference);
	sample_res.setDy(y_mass_test - y_mass_reference);

	res_evaluations++;

}


#endif
