// Description: Functor, create Jacobian matrix J
// Method  M7 (6 determined parameters)
// Jacobian matrix elements computed numerically using the Stirling method

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


#ifndef FJM7_H
#define FJM7_H

#include "libalgo/source/types/TVector.h"

#include "libalgo/source/algorithms/numderivative/NumDerivative.h"
#include "libalgo/source/algorithms/numderivative/FDiffM7.h"
#include "libalgo/source/algorithms/matrixoperations2/MatrixOperations.h"


template <typename T>
class FJM7
{
	private:
		const TVector <Point3DCartesian <T> > & test_points;				//List of test points
		const TVector <Point3DGeographic <T> > & reference_points;			//List of reference points
		const p_coord_function <T> getX, getY;						//Pointer to the coordinate functions
		const TTransformedLongitudeDirection trans_lon_dir;				//Transformed longitude direction


	public:
		FJM7( const TVector <Point3DCartesian <T> > & test_points_, const TVector <Point3DGeographic <T> > & reference_points_, const p_coord_function <T> &pX, const p_coord_function <T> &pY, const TTransformedLongitudeDirection trans_lon_dir_) :
			test_points(test_points_), reference_points(reference_points_), getX(pX), getY(pY), trans_lon_dir (trans_lon_dir_) {}
			
			
		void operator () (const Matrix <T> &X, Matrix <T> &J)
		{
			//Evaluate members of the Jacobian matrix, method M7
			Matrix <T> J_T = J;		//Temporary Jacobian matrix
			
			//Create matrix XT (1, 7) from X (7, 1)
			Matrix <T> XT = MatrixOperations::trans(X);

			//Process all points: compute Jacobian matrix of partial derivatives
			unsigned int i = 0, m = test_points.size();
			for (const auto p : reference_points)
			{
				//Get coordinates of the point
				const T lat = p.getLat();
				const T lon = p.getLon();
				
				//Upper part of the Jacobian matrix: R, latp, lonp, lat0, lon0=0
				J_T(i, 0) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, false);
				J_T(i, 1) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, false)/*  * 180 / M_PI */;
				J_T(i, 2) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i, 3) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i, 4) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i, 5) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i, 6) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getX, trans_lon_dir), XT, FirstDerivative, VariableX7, NUM_DERIV_STEP, false);

				//Lower part of the Jacobian matrix: R, latp, lonp, lat0, lon0=0
				J_T(i + m, 0) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX1, NUM_DERIV_STEP, false);
				J_T(i + m, 1) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX2, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i + m, 2) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX3, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i + m, 3) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX4, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i + m, 4) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX5, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i + m, 5) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX6, NUM_DERIV_STEP, false)/* * 180 / M_PI */;
				J_T(i + m, 6) = NumDerivative::getDerivative(FDiffM7 <T> (lat, lon, getY, trans_lon_dir), XT, FirstDerivative, VariableX7, NUM_DERIV_STEP, false);
				
				//Increment index
				i++;
			}

			//Compute column sums of the Jacobian matrix
			T s0x = 0.0, s1x = 0.0, s2x = 0.0, s3x = 0.0, s4x = 0.0, s5x = 0.0, s6x = 0.0, s0y = 0.0, s1y = 0.0, s2y = 0.0, s3y = 0.0, s4y = 0.0, s5y = 0.0, s6y = 0.0;

			for (i = 0; i < m; i++)
			{
				//Sums of X derivatives
				s0x += J_T(i, 0); 
				s1x += J_T(i, 1); 
				s2x += J_T(i, 2); 
				s3x += J_T(i, 3);
				s4x += J_T(i, 4);
				s5x += J_T(i, 5); 
				s6x += J_T(i, 6);

				//Sums of Y derivatives
				s0y += J_T(i + m, 0); 
				s1y += J_T(i + m, 1); 
				s2y += J_T(i + m, 2); 
				s3y += J_T(i + m, 3);
				s4y += J_T(i + m, 4); 
				s5y += J_T(i + m, 5);
				s6y += J_T(i + m, 6);
			}

			//Compute Jacobian matrix
			for (i = 0; i < m; i++)
			{
				//X derivatives
				J(i, 0) = J_T(i, 0) - s0x / m; 
				J(i, 1) = J_T(i, 1) - s1x / m; 
				J(i, 2) = J_T(i, 2) - s2x / m;
				J(i, 3) = J_T(i, 3) - s3x / m;
				J(i, 4) = J_T(i, 4) - s4x / m; 
				J(i, 5) = J_T(i, 5) - s5x / m;
				J(i, 6) = J_T(i, 6) - s6x / m;

				//Y derivatives
				J(i + m, 0) = J_T(i + m, 0) - s0y / m; 
				J(i + m, 1) = J_T(i + m, 1) - s1y / m; 
				J(i + m, 2) = J_T(i + m, 2) - s2y / m;
				J(i + m, 3) = J_T(i + m, 3) - s3y / m; 
				J(i + m, 4) = J_T(i + m, 4) - s4y / m; 
				J(i + m, 5) = J_T(i + m, 5) - s5y / m;
				J(i + m, 6) = J_T(i + m, 6) - s6y / m;
			}
		}

		
};

#endif