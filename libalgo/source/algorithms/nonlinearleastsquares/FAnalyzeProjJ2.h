// Description: Functor, create Jacobi matrix J for cartometric analysis, method M6 (6 determined parameters)
// Elements are computed numerically using the Stirling method

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


#ifndef FAnalyzeProjJ2_H
#define FAnalyzeProjJ2_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

#include "libalgo/source/algorithms/numderivative/FProjEquationDerivative6Var.h"

//Forward declarations
template <typename T>
class Projection;

template <typename T>
class Node3DCartesian;

template <typename T>
class Point3DGeographic;

template <typename T>
class FAnalyzeProjJ2
{
        private:
                //List of points
                const Container <Node3DCartesian <T> *> &nl_test;
                const Container <Point3DGeographic <T> *> &pl_reference;

                //Map projection and analyzed aspect
                const Projection <T> *proj;
                const TProjectionAspect aspect;

		bool &enable_additional_lon0_analysis;
                const bool print_exceptions;
		
		unsigned int &jac_evaluation;


        public:

		FAnalyzeProjJ2(const Container <Node3DCartesian <T> *> &nl_test_, const Container <Point3DGeographic <T> *> &pl_reference_, Projection <T> *proj_, const TProjectionAspect aspect_, bool &enable_additional_lon0_analysis_, const bool print_exceptions_, unsigned int &jac_evaluation_)
			: nl_test(nl_test_), pl_reference(pl_reference_), proj(proj_), aspect(aspect_), enable_additional_lon0_analysis(enable_additional_lon0_analysis_), print_exceptions(print_exceptions_), jac_evaluation(jac_evaluation_) {}


                void operator () ( const Matrix <T> &X, Matrix <T> &J )
                {
                        //Compute parameters of the Jacobi Matrix A
                        //Jacobi matrix J = [ d_R, d_latp, d_lonp, d_lat0, d_lon0, d_dx, d_dy, c]
                        const unsigned int m = nl_test.size();
                        unsigned int correct_derivatives = m;

                        //Create matrix XT (1, 6) from X ( transposed )
                        Matrix <T> XT = MatrixOperations::trans ( X );

                        //Get type of the direction
                        TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

                        //Create temporary A matrix
                        Matrix <T> J_T = J;

                        //Process all points: compute matrix of partial derivatives
                        for ( unsigned int i = 0; i < m; i++ )
                        {
                                try
                                {
					//Get lat, lob of the current point
					T lat = pl_reference[i]->getLat();
					T lon = pl_reference[i]->getLon();
	
					//Correct lat
					if (fabs(lat - MAX_LAT) < ANGLE_ROUND_ERROR)
						lat -=  3 * NUM_DERIV_STEP;
					else if (fabs(lat - MIN_LAT) < ANGLE_ROUND_ERROR)
						lat += 3 * NUM_DERIV_STEP;

					//Correct lon
					const T lon_red = CartTransformation::redLon0(lon, proj->getLon0());
					if (fabs(lon_red - MAX_LON) < ANGLE_ROUND_ERROR)
						lon -= 3 * NUM_DERIV_STEP;
					else if (fabs(lon_red - MIN_LON) < ANGLE_ROUND_ERROR)
						lon += 3 * NUM_DERIV_STEP;


                                        //Normal aspect: R, lat0, lon0, c
                                        if ( aspect == NormalAspect )
                                        {
                                                //Upper part of the matrix:  latp=90, lonp=0, lat0, lon0
						J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat, lon, proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat, lon, proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat, lon, proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  lat, lon, proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);


                                                //Lower part of the matrix: lat0, lon0: R, latp=90, lonp=0, lat0, lon0
						J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
                                        }

                                        //Transverse aspect: lonp, lat0
                                        else  if ( aspect == TransverseAspect )
                                        {
                                                //Upper part of the matrix: R, latp=0, lonp, lat0, lon0=0
						J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);

                                                //Lower part of the matrix: R, latp=0, lonp, lat0, lon0=0
						J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
                                        }

                                        //Oblique aspect: latp, lonp, lat0
                                        else
                                        {
                                                //Upper part of the matrix: R, latp, lonp, lat0, lon0=0
						J_T(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(), pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
						J_T(i, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
						
                                                //Lower part of the matrix: R, latp, lonp, lat0, lon0=0
						J_T(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
						J_T(i + m, 5) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, print_exceptions);
                                        }
                                }

                                //Decrease amount of corrected points
                                catch ( Exception & error )
                                {
                                        correct_derivatives --;
                                }
                        }

                        //Compute column sums
                        T sum0X = 0.0, sum1X = 0.0, sum2X = 0.0, sum3X = 0.0, sum4X = 0.0, sum5X = 0.0,
                          sum0Y = 0.0, sum1Y = 0.0, sum2Y = 0.0, sum3Y = 0.0, sum4Y = 0.0, sum5Y = 0.0;

                        for ( unsigned int i = 0; i < m; i++ )
                        {
                                //Sums of X derivatives
                                sum0X += J_T ( i, 0 );
                                sum1X += J_T ( i, 1 );
                                sum2X += J_T ( i, 2 );
                                sum3X += J_T ( i, 3 );
                                sum4X += J_T ( i, 4 );
                                sum5X += J_T ( i, 5 );

                                //Sums of Y derivatives
                                sum0Y += J_T ( i + m, 0 );
                                sum1Y += J_T ( i + m, 1 );
                                sum2Y += J_T ( i + m, 2 );
                                sum3Y += J_T ( i + m, 3 );
                                sum4Y += J_T ( i + m, 4 );
                                sum5Y += J_T ( i + m, 5 );
                        }

                        //Compute Jacobi matrix
                        for ( unsigned int i = 0; i < m; i++ )
                        {
                                //X derivatives
                                J ( i, 0 ) = J_T ( i, 0 ) - sum0X / m;
                                J ( i, 1 ) = J_T ( i, 1 ) - sum1X / m;
                                J ( i, 2 ) = J_T ( i, 2 ) - sum2X / m;
                                J ( i, 3 ) = J_T ( i, 3 ) - sum3X / m;
                                J ( i, 4 ) = J_T ( i, 4 ) - sum4X / m;
                                J ( i, 5 ) = J_T ( i, 5 ) - sum5X / m;

                                //Y derivatives
                                J ( i + m, 0 ) = J_T ( i + m, 0 ) - sum0Y / m;
                                J ( i + m, 1 ) = J_T ( i + m, 1 ) - sum1Y / m;
                                J ( i + m, 2 ) = J_T ( i + m, 2 ) - sum2Y / m;
                                J ( i + m, 3 ) = J_T ( i + m, 3 ) - sum3Y / m;
                                J ( i + m, 4 ) = J_T ( i + m, 4 ) - sum4Y / m;
                                J ( i + m, 5 ) = J_T ( i + m, 5 ) - sum5Y / m;
                        }
			
			//Normal aspect, 4 th column are of zeros, lon0 can not be determined by NLSP, call bisection
			if (aspect == NormalAspect)
			{
				const T sum2cols = sum2(J(0, 2 * m - 1, 4, 4));
				//std::cout << sum2cols << " ";

				//if (sum2cols < 1.0e-8)
					enable_additional_lon0_analysis = true;
				//else
					//enable_additional_lon0_analysis = false;
			}

                        //Not enough points
                        if ( correct_derivatives < 3 )
                        {
				X.print();
				pl_reference.print();
                                throw BadDataException ( "BadDataException: not enough correct partial derivatives, maybe error in equation. ", "Can not compute Jacobi matrix." );
                        }

			jac_evaluation++;
                }

};

#endif
