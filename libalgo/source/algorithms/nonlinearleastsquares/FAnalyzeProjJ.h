// Description: Functor, create Jacobi matrix J for cartometric analysis, mrthod M1 (8 determined parameters)
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


#ifndef FAnalyzeProjJ_H
#define FAnalyzeProjJ_H

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
class FAnalyzeProjJ
{
        private:
                //List of points
                const Container <Node3DCartesian <T> *> &nl_test;
                const Container <Point3DGeographic <T> *> &pl_reference;

                //Map projection and analyzed aspect
                const Projection <T> *proj;
                const TProjectionAspect aspect;

                const bool print_exceptions;

		unsigned int &iter;

        public:

                FAnalyzeProjJ ( const Container <Node3DCartesian <T> *> &nl_test_, const Container <Point3DGeographic <T> *> &pl_reference_, Projection <T> *proj_, const TProjectionAspect aspect_, const bool print_exceptions_, unsigned int &iter_ )
			: nl_test(nl_test_), pl_reference(pl_reference_), proj(proj_), aspect(aspect_), print_exceptions(print_exceptions_), iter(iter_) {}


                void operator () ( const Matrix <T> &X, Matrix <T> &J )
                {
                        //Compute parameters of the Jacobi Matrix J
                        //Jacobian J = [ d_R, d_latp, d_lonp, d_lat0, d_lon0, d_dx, d_dy]
                        const unsigned int m = nl_test.size();
                        unsigned int correct_derivatives = m;

                        //Create matrix XT (1, 7) from X ( transposed )
                        Matrix <T> XT = trans ( X );

                        //Get type of the direction
                        TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

                        //Process all points: compute Jacobi matrix
                        for ( unsigned int i = 0; i < m; i++ )
                        {
                                try
                                {
                                        //Normal aspect: lat0, lon0
                                        if ( aspect == NormalAspect )
                                        {
						J(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>( proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J(i, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
                                                J ( i, 5 ) = 1.0;
                                                J ( i, 6 ) = 0.0;
						J(i, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);

                                                //Lower part of the matrix: lat0, lon0: R, latp=90, lonp=0, lat0, lon0, dx, dy, c
						J(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 4) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, print_exceptions);
                                                J ( i + m, 5 ) = 0.0;
                                                J ( i + m, 6 ) = 1.0;
						J(i + m, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);
                                        }

                                        //Transverse aspect: lonp, lat0
                                        else  if ( aspect == TransverseAspect )
                                        {
                                                //Upper part of the matrix: R, latp=0, lonp, lat0, lon0=0, dx, dy, c
						J(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
                                                J ( i, 5 ) = 1.0;
                                                J ( i, 6 ) = 0.0;
						J(i, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);

                                                //Lower part of the matrix: R, latp=0, lonp, lat0, lon0=0, dx, dy, c
						J(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
                                                J ( i + m, 5 ) = 0.0;
                                                J ( i + m, 6 ) = 1.0;
						J(i + m, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);
                                        }

                                        //Oblique aspect: latp, lonp, lat0
                                        else
                                        {
                                                //Upper part of the matrix: R, latp, lonp, lat0, lon0=0, dx, dy, c
						J(i, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
						J(i, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J(i, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
                                                J ( i, 5 ) = 1;
                                                J ( i, 6 ) = 0;
						J(i, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getXEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);

                                                //Lower part of the matrix: R, latp, lonp, lat0, lon0=0, dx, dy
						J(i + m, 0) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 1) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 2) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, print_exceptions);
						J(i + m, 3) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, print_exceptions);
                                                J ( i + m, 5 ) = 0;
                                                J ( i + m, 6 ) = 1;
						J(i + m, 7) = NumDerivative::getDerivative(FProjEquationDerivative6Var <T>(proj->getYEquatPostfix(), proj->getFThetaEquatPostfix(), proj->getTheta0EquatPostfix(),  pl_reference[i]->getLat(), pl_reference[i]->getLon(), proj->getA(), proj->getB(), proj->getLat1(), proj->getLat2(), trans_lon_dir), XT, FirstDerivative, VariableX8, NUM_DERIV_STEP, print_exceptions);
                                        }
                                }

                                //Decrease amount of corrected points
                                catch ( Exception & error )
                                {
                                        correct_derivatives --;
                                }
                        }

                        //Not enough points
                        if ( correct_derivatives < 3 )
                        {
                                throw BadDataException ( "BadDataException: not enough correct partial derivatives, maybe error in equation. ", "Can not compute Jacobi matrix." );
                        }

			iter++;
                }

};

#endif
