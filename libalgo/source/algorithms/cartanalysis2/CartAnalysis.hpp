// Description: Performs cartometric analysis (i.e. detection of the maps projection)

// Copyright (c) 2015 - 2016
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


#ifndef CartAnalysis_HPP
#define CartAnalysis_HPP


#include "libalgo/source/comparators2/sortPointsByLat.h"
#include "libalgo/source/comparators2/sortPointsByLon.h"

#include "libalgo/source/algorithms/nonlinearleastsquares2/NonLinearLeastSquares.h"
#include "libalgo/source/algorithms/nonlinearleastsquares2/FJM7.h"
#include "libalgo/source/algorithms/nonlinearleastsquares2/FJM8.h"
#include "libalgo/source/algorithms/nonlinearleastsquares2/FRM7.h"
#include "libalgo/source/algorithms/nonlinearleastsquares2/FRM8.h"
#include "libalgo/source/algorithms/geneticalgorithms/FRM7DE.h"
#include "libalgo/source/algorithms/geneticalgorithms/FRM8DE.h"
#include "libalgo/source/algorithms/simplexmethod/FRM7NM.h"
#include "libalgo/source/algorithms/simplexmethod/FRM8NM.h"
#include "libalgo/source/algorithms/geneticalgorithms/DifferentialEvolution.h"
#include "libalgo/source/algorithms/simplexmethod/SimplexMethod.h"
#include "libalgo/source/algorithms/transformation2/HelmertTransformation2D.h"

#include "libalgo/source/io/Format.h"


template <typename T>
void CartAnalysis::analyzeProjection(const TVector<Point3DCartesian<T> > &test_points, const TVector <Point3DGeographic<T> > &reference_points, 
	TListS <Projection<T> > &projections, TResults <T> &results, const TAnalysisMethod &method, std::ostream &output)
{
	//Analyze all projections
	const unsigned int n = test_points.size();
	
	//Geographic extent of the analyzed territory
	const T lat_min = (std::min_element(reference_points.begin(), reference_points.end(), sortPointsByLat()))->getLat();
	const T lat_max = (std::max_element(reference_points.begin(), reference_points.end(), sortPointsByLat()))->getLat();
	const T lon_min = (std::min_element(reference_points.begin(), reference_points.end(), sortPointsByLon()))->getLon();
	const T lon_max = (std::max_element(reference_points.begin(), reference_points.end(), sortPointsByLon()))->getLon();
	const T lat_aver = 0.5 * (lat_min + lat_max), lon_aver = 0.5 * (lon_min + lon_max);

	//Initialize cartographic parameters
	const T lat1 = (lat_aver == 0 ? 10 : lat_aver);
	const T lat2 = lat1 + 10;
	const T latp = std::min(90.0 - fabs(lat_aver), 80.0);
	const T lonp = (lon_aver > 90 ? lon_aver - 270 : (latp == 90 ? 0 : lon_aver + 90));
	const T lon0 =  lon_aver;

	//Process all projections
	for (const auto proj : projections)
	{
		//Set parameters of the projection
		proj->setCartPole(Point3DGeographic <T>(latp, lonp));
		proj->setLat1(lat1);
		proj->setLat2(lat2);
		proj->setLon0(lon0);
		proj->setC(1.0);
		std::cout << proj->getName();
		
		try
		{
			//Get initial scale
			TVector <Point3DCartesian<T> > reference_points_projected;

			//Apply projection proj(Q->P'): convert geographic points to the cartesian 
			for (const auto p : reference_points)
			{
				//(lat, lon) -> (lat_trans, lon_trans)
				const T lat_trans = CartTransformation::latToLatTrans(p.getLat(), p.getLon(), proj->getCartPole().getLat(), proj->getCartPole().getLon());
				const T lon_trans = CartTransformation::lonToLonTrans(p.getLat(), p.getLon(), proj->getCartPole().getLat(), proj->getCartPole().getLon(), proj->getLonDir());

				//Reduce longitude
				const T lon_transr = CartTransformation::redLon0(lon_trans, proj->getLon0());

				//(lat_trans, lon_trans) -> (X, Y)
				const T XR = proj->getX(lat_trans, lon_transr);
				const T YR = proj->getY(lat_trans, lon_transr);

				//Add point to the list
				reference_points_projected.push_back(Point3DCartesian<T>(XR, YR));
			}

			//Compute 2D Helmert transformation between P, P'
			TVector <T> weights(test_points.size(), 1.0);
			TTransformationKeyHelmert2D <T> key;
			HelmertTransformation2D::getTransformKey(test_points, reference_points_projected, weights, key);
			T R_0 = sqrt(key.c1 * key.c1 + key.c2 * key.c2) * proj->getR();

			//Remember old radius
			const T R = proj->getR();

			//Set new radius
			proj->setR(R_0);

			//Get initial matrices: X, A, B
			unsigned int  m = (method < 3 ? 7 : 6);
			Matrix <T> X = (method < 3 ? X0M7<T>(m, proj, method) : X0M8<T>(m, proj, method));
			Matrix <T> A = (method < 3 ? AM7<T>(m, proj, R_0) : AM8<T>(m, proj));
			Matrix <T> B = (method < 3 ? BM7<T>(m, proj, R_0) : BM8<T>(m, proj));

			//Transpose matrices (Differential evolution, Nelder-Mead)
			Matrix <T> XT = MatrixOperations::trans(X), AT = MatrixOperations::trans(A), BT = MatrixOperations::trans(B);
			Matrix <T> XAVER = XT;

			//Apply optimization method
			unsigned int iterations = 0, population = 2 * m * A.rows(), max_iter_nls = 100, max_gen = 75, max_iter_nm = 500;
			const T alpha = 0.0001, nu = 0.0001, max_error = 1.0e-10, eps = 1.0e-10, CR = 0.8;
			T q1 = key.c1, q2 = key.c2, dx = 0, dy = 0, min_cost = MAX_FLOAT, res_aver = 0, res_max = 0;
			
			Matrix <T> W(2 * n, 2 * n, 0.0, 1), Y(2 * n, 1), V(2 * n, 1), F(1, 1); F(0, 0) = 0.5;

			// Method M7, Non-linear least squares
			if (method == NLSM7)
			{
				min_cost = NonLinearLeastSquares::BFGSH(FJM7<T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir()), FRM7 <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), dx, dy), W, X, Y, V, A, B, iterations, alpha, nu, max_error,  max_iter_nls);
			}

			// Method M8, Non-linear least squares
			else if (method == NLSM8)
			{
				min_cost = NonLinearLeastSquares::BFGSH(FJM8<T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), R_0, q1, q2), FRM8 <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), R_0, q1, q2, dx, dy), W, X, Y, V, A, B, iterations, alpha, nu, max_error, max_iter_nls);
			}

			// Method M7, Differential evolution
			else if (method == DEM7)
			{
				min_cost = DifferentialEvolution::diffEvolution(FRM7DE <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), dx, dy), population, eps, max_gen, F, CR, DERand1Strategy, MFDE, W, XT, Y, V, AT, BT, XAVER, res_aver, res_max, iterations, true);
				X = MatrixOperations::trans(XT);
			}

			// Method M8, Differential evolution
			else if (method == DEM8)
			{
				min_cost = DifferentialEvolution::diffEvolution(FRM8DE <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), R_0, q1, q2, dx, dy), population, eps, max_gen, F, CR, DERand1Strategy, MFDE, W, XT, Y, V, AT, BT, XAVER, res_aver, res_max, iterations, true);
				X = MatrixOperations::trans(XT);
			}

			// Method M7, Nelder-Mead optimization
			else if (method == NMM7)
			{
				min_cost = SimplexMethod::NelderMead(FRM7NM <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), dx, dy), W, XT, Y, V, AT, BT, iterations, eps, max_iter_nm, true);
				X = MatrixOperations::trans(XT);
			}

			// Method M8, Nelder-Mead optimization
			else if (method == NMM8)
			{
				min_cost = SimplexMethod::NelderMead(FRM8NM  <T>(test_points, reference_points, proj->getX(), proj->getY(), proj->getLonDir(), R_0, q1, q2, dx, dy), W, XT, Y, V, AT, BT, iterations, eps, max_iter_nm, true);
				X = MatrixOperations::trans(XT);
			}

			//Set parameters to the projections
			method < 3 ? proj->setR(X(0, 0)) : proj->setR(R_0);
			method < 3 ? proj->setCartPole(Point3DGeographic<T>(X(1, 0), X(2, 0))) : proj->setCartPole(Point3DGeographic<T>(X(0, 0), X(1, 0)));
			method < 3 ? proj->setLat1(X(3, 0)) : proj->setLat1(X(2, 0));
			method < 3 ? proj->setLat2(X(4, 0)) : proj->setLat2(X(3, 0));
			method < 3 ? proj->setLon0(X(5, 0)) : proj->setLon0(X(4, 0));
			method < 3 ? proj->setC(X(6, 0)) : proj->setC(X(5, 0));
			proj->setDx(dx);
			proj->setDy(dy);

			//Set parameters of the map
			const T map_scale = R / proj->getR() * 1000 ;
			const T rotation =  (method < 3 ? 0 : atan2(q2, q1) * 180.0 / M_PI);

			//Add result to the list of results
			TResult<T> res = { proj, map_scale, rotation, iterations };
			results[min_cost] = res;
		}

		catch (Exception &e)
		{
			e.printException();
			std::cout << proj->getName();
		}
	}
}


template <typename T>
Matrix <T> CartAnalysis:: X0M7(const unsigned int m, const std::shared_ptr<Projection <T> > proj, const TAnalysisMethod &method)
{
	Matrix <T> X0(m, 1);

	//Set initial solution
	X0(0, 0) = proj->getR();
	X0(1, 0) = (method < 2 ? proj->getCartPole().getLat() : MAX_LAT); //Add initial solution to the simplex/population
	X0(2, 0) = (method < 2 ? proj->getCartPole().getLon() : 10.0); //Add initial solution to the simplex/population
	X0(3, 0) = proj->getLat1();
	X0(4, 0) = proj->getLat2();
	X0(5, 0) = proj->getLon0();
	X0(6, 0) = proj->getC();

	return X0;
}


template <typename T>
Matrix <T> CartAnalysis::X0M8(const unsigned int m, const std::shared_ptr<Projection <T> > proj, const TAnalysisMethod &method )
{
	Matrix <T> X0(m, 1);

	//Set initial solution
	X0(0, 0) = (method < 2 ? proj->getCartPole().getLat() : MAX_LAT); //Add initial solution to the simplex/population
	X0(1, 0) = (method < 2 ? proj->getCartPole().getLon() : 10.0); //Add initial solution to the simplex/population
	X0(2, 0) = proj->getLat1();
	X0(3, 0) = proj->getLat2();
	X0(4, 0) = proj->getLon0();
	X0(5, 0) = proj->getC();

	return X0;
}


template <typename T>
Matrix <T> CartAnalysis::AM7(const unsigned int m, const std::shared_ptr <Projection <T> > proj, const T & scale)
{
	//Return matrix of the lower bounds
	Matrix <T> A(m, 1);

	//Set initial solution
	A(0, 0) = 0.01 * scale;
	A(1, 0) = -90;
	A(2, 0) = -180;
	A(3, 0) = proj->getLat1Interval().min;
	A(4, 0) = proj->getLat1Interval().min;
	A(5, 0) = -180;
	A(6, 0) = 0;

	return A;
}


template <typename T>
Matrix <T> CartAnalysis::BM7(const unsigned int m, const std::shared_ptr <Projection <T> > proj, const T & scale)
{
	//Return matrix of the upper bounds
	Matrix <T> B(m, 1);

	//Set initial solution
	B(0, 0) = 100.0 * scale;
	B(1, 0) = 90;
	B(2, 0) = 180;
	B(3, 0) = proj->getLat1Interval().max;
	B(4, 0) = proj->getLat1Interval().max;
	B(5, 0) = 180;
	B(6, 0) = 100;

	return B;
}


template <typename T>
Matrix <T> CartAnalysis::AM8(const unsigned int m, const std::shared_ptr <Projection <T> > proj)
{
	Matrix <T> A(m, 1);

	//Set initial solution
	A(0, 0) = -90;
	A(1, 0) = -180;
	A(2, 0) = proj->getLat1Interval().min;
	A(3, 0) = proj->getLat1Interval().min;
	A(4, 0) = -180;
	A(5, 0) = 0;

	return A;
}


template <typename T>
Matrix <T> CartAnalysis::BM8(const unsigned int m, const std::shared_ptr <Projection <T> > proj)
{
	Matrix <T> B(m, 1);

	//Set initial solution
	B(0, 0) = 90;
	B(1, 0) = 180;
	B(2, 0) = proj->getLat1Interval().max;
	B(3, 0) = proj->getLat1Interval().max;
	B(4, 0) = 180;
	B(5, 0) = 100;

	return B;
}

template <typename T>
void CartAnalysis::printResults(TResults <T> &results, std::ostream &output)
{
	//Print results
	output << std::showpoint << std::fixed << std::right;

	//Create header
	output << std::endl
		<< std::setw(4) << "#"
		<< std::setw(13) << "Family"
		<< std::setw(8) << "Proj."
		<< std::setw(9) << "Resid."
		<< std::setw(12) << "R"
		<< std::setw(7) << "latP"
		<< std::setw(7) << "lonP"
		<< std::setw(7) << "lat1"
		<< std::setw(7) << "lat2"
		<< std::setw(7) << "lon0"
		<< std::setw(7) << "C"
		<< std::setw(15) << "Scale"
		<< std::setw(10) << "Rotation"
		<< std::setw(10) << "Iter."
		<< std::endl;

	//Print all candidates
	unsigned int index = 1;
	for (const auto res : results)
	{
		//Get projection
		const std::shared_ptr <Projection <T>> proj = res.second.proj;

		//Set parameters
		output << std::setw(4) << std::fixed << std::right;

		//Print results
		output << index
			<< std::setw(13) << proj->getFamily()
			<< std::setw(8) << proj->getName()
			<< std::setw(9) << std::setprecision(2) << Format::modScientific(res.first)
			<< std::setw(12) << std::setprecision(3) << proj->getR()
			<< std::setw(7) << std::setprecision(1) << proj->getCartPole().getLat()
			<< std::setw(7) << std::setprecision(1) << proj->getCartPole().getLon()
			<< std::setw(7) << std::setprecision(1) << proj->getLat1()
			<< std::setw(7) << std::setprecision(1) << proj->getLat2()
			<< std::setw(7) << std::setprecision(1) << proj->getLon0()
			<< std::setw(7) << std::setprecision(1) << proj->getC()
			<< std::setw(15) << std::setprecision(1) << res.second.map_scale
			<< std::setw(10) << std::setprecision(2) << res.second.map_rotation
			<< std::setw(10) << res.second.iterations
			<< std::endl;

		// output << proj->getName() << '\t' << res.first << '\t' << proj->getR() << '\t' << proj->getCartPole().getLat() << '\t' << proj->getCartPole().getLon() << '\t' << proj->getLat1() << '\t' << proj->getLat2() << '\t' << proj->getC() << '\t' << res.second.map_scale << '\t' << res.second.map_rotation << '\n';

		index++;
	}
}


#endif