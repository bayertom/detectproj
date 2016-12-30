// Description: Functor, compute matrix V of squares of residuals for cartometric analysis
// Method: Simplex Method with rotation

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


#ifndef FAnalyzeProjV3S_H
#define FAnalyzeProjV3S_H

template <typename T>
class Sample;

#include "libalgo/source/structures/projection/Sample.h"

//#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


//Forward declarations
template <typename T>
class Projection;


template <typename T>
struct TAnalysisParameters;


template <typename T>
class FAnalyzeProjV3S
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
	unsigned int &iter;
	Matrix <T> &XMIN, &XMAX;
	std::ostream * output;


public:

	FAnalyzeProjV3S(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &iter_, Matrix <T> &XMIN_, Matrix <T> &XMAX_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), iter(iter_), XMIN(XMIN_), XMAX(XMAX_), output(output_) {}

	void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis)
	{

		//Compute squares of residuals V = v' * W *v
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
				//Correct R, lat0, lon0
				if (X(i, 0) < 0.0) X(i, 0) = fabs(X(i, 0));

				//Subtract period
				if (X(i, 3) < lat0_min || X(i, 3) > lat0_max) X(i, 3) = 0.5 * (lat0_min + lat0_max);

				//Set lon0
				if (X(i, 4) < MIN_LON)  X(i, 4) = MAX_LON + fmod(X(i, 4), MIN_LON);
				else if (X(i, 4) > MAX_LON)  X(i, 4) = MIN_LON + fmod(X(i, 4), MAX_LON);

				//Subtract period
				if (fabs(X(i, 6)) > MAX_LAT) X(i, 6) = fmod(X(i, 6), 180);
			}

			//Transverse aspect: lonp, lat0
			else  if (aspect == TransverseAspect)
			{
				//Correct R, lonp, lat0
				if (X(i, 0) < 0.0) X(i, 0) = fabs(X(i, 0));

				//Subtract period
				if (X(i, 2) < MIN_LON)  X(i, 2) = MIN_LON - fmod(X(i, 2), MIN_LON);
				else if (X(i, 2) > MAX_LON)  X(i, 2) = MAX_LON - fmod(X(i, 2), MAX_LON);

				//Set to interval
				if (X(i, 3) < lat0_min || X(i, 3) > lat0_max) X(i, 3) = 0.5 * (lat0_min + lat0_max);

				//Set lon0
				X(i, 4) = 0;

				//Subtract period
				if (fabs(X(i, 6)) > MAX_LAT) X(i, 6) = fmod(X(i, 6), 180);
			}

			//Oblique aspect: latp, lonp, lat0
			else if (aspect == ObliqueAspect)
			{
				//Correct R, latp, lonp, lat0
				if (X(i, 0) < 0.0) X(i, 0) = fabs(X(i, 0));

				//Subtract period
				if(X(i, 1) < MIN_LAT)  X(i, 1) = MIN_LAT - fmod(X(i, 1), MIN_LAT);
				else if (X(i, 1) > MAX_LAT)  X(i, 1) = MAX_LAT - fmod(X(i, 1), MAX_LAT);

				if (X(i, 2) < MIN_LON)  X(i, 2) = MAX_LON + fmod(X(i, 2), MIN_LON);
				else if (X(i, 2) > MAX_LON)  X(i, 2) = MIN_LON + fmod(X(i, 2), MAX_LON);

				//Set lat0 inside the interval
				//if (X(i, 3) < lat0_min || X(i, 3) > lat0_max) X(i, 3) = 0.5 * (lat0_min + lat0_max);
				
				//Set lonp to zero, if latp = 90
				if (fabs(fabs(X(i, 1)) - MAX_LAT) <= 3.0)
				{
					//X(i, 1) = 90.0;
					//X(i, 2) = 0.0;
				}
				
				//Set lon0
				X(i, 4) = 0;

				//Subtract period
				if (fabs(X(i, 6)) > MAX_LON) X(i, 6) = fmod(X(i, 6), 180);
			}


			//Set properties to the projection: ommit estimated radius, additional constants dx, dy
			// They will be estimated again using the transformation
			Point3DGeographic <T> cart_pole(X(i, 1), X(i, 2));
			proj->setCartPole(cart_pole);
			proj->setLat0(X(i, 3));
			proj->setLon0(X(i, 4));
			proj->setDx(0.0);
			proj->setDy(0.0);
			proj->setC(X(i, 5));
			proj->setR(X(i, 0));

			//Get alpha
			const T alpha = X(i, 6);

			//Compute coordinate differences (residuals): items of V matrix
			Container <Node3DCartesianProjected <T> *> nl_projected_temp;

			for (unsigned int j = 0; j < m; j++)
			{
				//Get type of the direction
				TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

				//Reduce lon
				const T lon_red = CartTransformation::redLon0(pl_reference[j]->getLon(), X(i, 4));

				T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

				try
				{
					//Convert geographic point to oblique aspect
					lat_trans = CartTransformation::latToLatTrans(pl_reference[j]->getLat(), lon_red, X(i, 1), X(i, 2));
					lon_trans = CartTransformation::lonToLonTrans(pl_reference[j]->getLat(), lon_red, X(i, 1), X(i, 2), trans_lon_dir);

					for (unsigned int k = 0; k < 3; k++)
					{
						try
						{
							//Compute x, y coordinates
							x = CartTransformation::latLonToX(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(i, 0), proj->getA(), proj->getB(), 0.0, X(i, 5), X(i, 3), proj->getLat1(), proj->getLat2(), false);
							y = CartTransformation::latLonToY(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat_trans, lon_trans, X(i, 0), proj->getA(), proj->getB(), 0.0, X(i, 5), X(i, 3), proj->getLat1(), proj->getLat2(), false);
							
							//x = ArithmeticParser::parseEquation(proj->getXEquat(), lat_trans, lon_trans, X(i, 0), proj->getA(), proj->getB(), X(i, 5), X(i, 3), proj->getLat1(), proj->getLat2(), false);
							//y = ArithmeticParser::parseEquation(proj->getYEquat(), lat_trans, lon_trans, X(i, 0), proj->getA(), proj->getB(), X(i, 5), X(i, 3), proj->getLat1(), proj->getLat2(), false);
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

				//Bad conversion: singular point
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

			//Computer centers of mass
			unsigned int n_points = 0;
			T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

			for (unsigned int j = 0; j < m; j++)
			{
				//Use only non singular points
				if (W(j, j) != 0.0)
				{
					x_mass_test += nl_test[j]->getX();
					y_mass_test += nl_test[j]->getY();

					x_mass_reference += nl_projected_temp[j]->getX();
					y_mass_reference += nl_projected_temp[j]->getY();

					n_points++;
				}
			}

			x_mass_test = x_mass_test / n_points;
			y_mass_test = y_mass_test / n_points;
			x_mass_reference = x_mass_reference / n_points;
			y_mass_reference = y_mass_reference / n_points;

			//Compute coordinate differences (residuals): estimated - input
			Matrix <T> v(2 * m, 1);

			for (unsigned int j = 0; j < m; j++)
			{
				//Use only non singular points
				if (W(j, j) != 0.0)
				{
					const T vx = ((nl_projected_temp[j]->getX() - x_mass_reference) * cos(alpha * M_PI / 180) - (nl_projected_temp[j]->getY() - y_mass_reference) * sin(alpha * M_PI / 180) - (nl_test[j]->getX() - x_mass_test));
					const T vy = ((nl_projected_temp[j]->getX() - x_mass_reference) * sin(alpha * M_PI / 180) + (nl_projected_temp[j]->getY() - y_mass_reference) * cos(alpha * M_PI / 180) - (nl_test[j]->getY() - y_mass_test));

					//Input is the best vector
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

			//Compute squares of residuals: only if input is a simplex
			if (((m1 > 1) || (m2 == 1)))
			{
				V(i, 0) = MatrixOperations::norm(MatrixOperations::trans(v) * W * v);
			}

			//Set DX, DY, compute it from the first vertex of the simplex
			if (i == 0)
			{
				//Compute DX, DY
				T dx = x_mass_test - x_mass_reference * cos(alpha * M_PI / 180) + y_mass_reference * sin(alpha * M_PI / 180);
				T dy = y_mass_test - x_mass_reference * sin(alpha * M_PI / 180) - y_mass_reference * cos(alpha * M_PI / 180);

				//Set DX, DY
				sample_res.setDx(dx);
				sample_res.setDy(dy);
			}

			iter++;
		}
	}


};


#endif
