// Description: Downhill simplex optimalization method

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


#ifndef SimplexMethod_HPP
#define SimplexMethod_HPP

#include <time.h>
 
//Include  C++98/C++11 version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/algorithms/nonlinearleastsquares/NonLinearLeastSquares.h"
#else
	#include "libalgo/source/algorithms/nonlinearleastsquares2/NonLinearLeastSquares.h"
#endif

        
//Set namespace
using namespace MatrixOperations;


template <typename T, typename Function>
T SimplexMethod::NelderMead ( Function function, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, const Matrix <T> &XMIN, const Matrix <T> &XMAX, unsigned int &iterations, const T max_error, const unsigned int max_iterations, const bool add_x0, std::ostream * output )
{
        //Compute Nelder-Mead method for a function
        //Algorithm based on Lagarias, Reeds, Wright, 1998, SIAM
	const unsigned int n = XMIN.cols(), m = n + 1, m1 = W.rows();
	const T Rho = 1.0, Xi = 2.0, Gama = 0.5, Sigma = 0.5;

	Matrix <T> XX(m, n);

	//Add the initial solution to the simplex
	if (add_x0)
		XX(X, 0, 0);

	//Create random simplex
        createRandSimplex ( XMIN, XMAX, XX, add_x0 );

        //Create matrices for the simplex operations
        Matrix <unsigned int> IX ( m, 1 );
        Matrix <T> VV2 ( m, 1 ), VR ( m1, 1 ), VS ( m1, 1 ), VE ( m1, 1 ), VCO ( m1, 1 ), VCI ( m1, 1 ), YR ( m1, 1 ), YS ( m1, 1 ),
		YE(m1, 1), YCO(m1, 1), YCI(m1, 1), WR(m1, m1, 0, 1), WS(m1, m1, 0, 1), WE(m1, m1, 0, 1), WCO(m1, m1, 0, 1), WCI(m1, m1, 0, 1);
	Matrix <T> YBEST(m1, 1), VBEST(1, 1), WBEST(m1, m1, 0, 1);

        //Set iterations to 0
        iterations = 0;

        //Compute initial function values
	//Matrix <T> VVV = VV2;
	for (unsigned int i = 0; i < m; i++)
	{
		Matrix <T> VVI(m1, 1), YI(m1, 1), WI(m1, m1, 0, 1);
		Matrix <T> XXI = XX.row(i);

		function(XXI, YI, VVI, WI);
		const T res2 = norm(trans(VVI) * W * VVI);
		VV2(i, 0) = res2;
	}

        //Sort residuals in ascending order
        sort ( VV2, IX );

        //Change order of XX rows: row i to row IX (i)
        Matrix <T> XS ( m, n );

        for ( unsigned int i = 0; i < m; i++ )
                XS.row ( XX.row ( IX ( i, 0 ) ), i );

        //Assign sorted matrix to XX
        XX = XS;

        //Perform Nelder-Mead algorithm
        do
        {
                //Compute centroid
                Matrix <T> XC = sumCols ( XX ( 0, n - 1, 0, n - 1 ) ) * 1.0 / n;

                //Compute reflection point and residuals
                Matrix <T> XR = XC * ( 1.0 + Rho ) - XX ( n, n, 0, n - 1 ) * Rho;

                //Reflection into the search space
                reflection ( XR, XMIN, XMAX);

                //Compute residuals
                function ( XR, YR, VR, WR );
                const T  fr = norm(trans(VR) * W * VR);

                //A reflection point acceptable
                if ( ( VV2 ( 0, 0 ) <= fr ) && ( fr < VV2 ( n - 1, 0 ) ) )
                {
                        XX ( XR, n, 0 );
                        VV2 ( n, 0 ) = fr;
                }

                //Expansion of the simplex
                else if ( fr < VV2 ( 0, 0 ) )
                {
                        //Compute expanded point
                        Matrix <T> XE = XC * ( 1.0 + Rho * Xi ) - XX ( n, n, 0, n - 1 ) * Rho * Xi;

                        //Reflection into the search space
                        reflection (XE, XMIN, XMAX );

                        //Compute residuals
                        function ( XE, YE, VE, WE );
                        const T  fe = norm(trans(VE) * W * VE);

                        //An expanded point is acceptable
                        if ( fe < fr )
                        {
                                XX ( XE, n, 0 );
                                VV2 ( n, 0 ) = fe;
                        }

                        //An expanded point is not acceptable (use a reflected one)
                        else
                        {
                                XX ( XR, n, 0 );
                                VV2 ( n, 0 ) = fr;
                        }
                }

                //Outside contraction of the simplex
                else if ( ( VV2 ( n - 1, 0 ) <= fr ) && ( fr < VV2 ( n , 0 ) ) )
                {
                        //Compute outside contracted point
                        Matrix <T> XCO = XC * ( 1.0 + Rho * Gama ) - XX ( n, n, 0, n - 1 ) * Rho * Gama;

                        //Reflection into the search space
                        reflection ( XCO, XMIN, XMAX);

                        //Compute residuals
                        function ( XCO, YCO, VCO, WCO );
                        const T  fco = norm(trans(VCO) * W * VCO);

                        //An outside contracted point is acceptable
                        if ( fco < VV2 ( n , 0 ) )
                        {
                                XX ( XCO, n, 0 );
                                VV2 ( n, 0 ) = fco;
                        }

                        //An outside contracted point is not acceptable: shrink a simplex
                        else
                        {
                                shrink ( function, W, XX, XMIN, XMAX, Y, VV2, Sigma );

                                //Increment iterations
                                iterations ++;
                        }
                }

                //Inside contraction of the simplex
                else
                {
                        //Compute outside contracted point
                        Matrix <T> XCI = XC * ( 1.0 - Gama ) + XX ( n, n, 0, n - 1 ) * Gama;

                        //Reflection into the search space
                        reflection (XCI, XMIN, XMAX);

                        //Compute residuals
                        function ( XCI, YCI, VCI, WCI );
                        const T  fci = norm(trans(VCI) * W * VCI);

                        //An inside contracted point is acceptable
                        if ( fci < VV2 ( n , 0 ) )
                        {
                                XX ( XCI, n, 0 );
                                VV2 ( n, 0 ) = fci;
                        }

                        //An inside contracted point is not acceptable: shrink a simplex
                        else
                        {
                                shrink ( function, W, XX, XMIN, XMAX, Y, VV2,  Sigma );

                                //Increment iterations
                                iterations ++;
                        }
                }

                //Sort residuals in ascending order
                sort ( VV2, IX );

                //Change order of XX rows: row i to row IX (i)
                for ( unsigned int i = 0; i < m; i++ )
                        XS.row ( XX.row ( IX ( i, 0 ) ), i );

                //Assign sorted matrix to XX
                XX = XS;

                //Increment iterations
                iterations ++;

                if ( iterations % 100 == 0 )
                {
                        std::cout.flush();
                        std::cout << ".";
                }
        }
        while ( ( iterations < max_iterations ) && ( fabs ( VV2 ( 0, 0 ) - VV2 ( n, 0 ) ) > max_error ) );

        //Get minimum
        X = XX ( 0, 0, 0, n - 1 );

        //Compute residuals for the found solution
        function ( X, Y, V, W );
	const T fx_min = norm(trans(V) * W * V);

	std::cout << " [" << iterations << " it., fmin = " << fx_min << "]" << '\n';
        //Print residuals
        //*output << " X:"; X.print ( output );
        //*output << "Residuals: "; V.print ( output );
        //*output << "Iterations: " << iterations << '\n';

        //Return squares of residuals
	return fx_min;
}


template <typename T>
void SimplexMethod::createRandSimplex ( const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &XX, const bool add_x0)
{
        //Create random simplex
        const unsigned int dim = XMIN.cols();

        //Initialize random number generator
        srand ( ( unsigned ) time ( 0 ) );

	//Compute difference
	const Matrix <T> DX = XMAX - XMIN;

        //Create random simplex
        for ( unsigned int i = 0; i < dim + 1; i++ )
        {
		if ((i != 0) && (add_x0) || (!add_x0))
		{
			for (unsigned int j = 0; j < dim; j++)
			{
				const T rand_num = XMIN(0, j) + DX(0, j) * rand() / (RAND_MAX + 1.0);
				XX(i, j) = rand_num;
			}
		}
        }
}



template <typename T, typename Function>
void SimplexMethod::shrink ( Function function, Matrix <T> &W, Matrix <T> &XX, const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &Y,  Matrix <T> &VV2, const T Sigma )
{
        //Shrink a simplex
        const unsigned int m = XX.rows(), n = XX.cols(), m1 = W.rows();

        //Create matrices
        Matrix <T> V1 ( m1, 1 ), VSH ( m1, 1 ), Y1 ( m1, 1 ), YSH ( m1, 1 ), W1 ( m1, m1, 0, 1 ), WSH ( m1, m1, 0, 1 );

        //Get first point of the simplex (best)
        Matrix <T> X1 = XX ( 0, 0, 0, n - 1 );

        //Compute residuals
        function ( X1, Y1, V1, W1 );
	const T  fv1 = norm(trans(V1) * W1 * V1);

        //Actualize VV matrix
        VV2 ( 0, 0 ) = fv1;

        //Shrink remaining simplex points
        for ( unsigned int i = 1; i < n + 1 ; i++ )
        {
                //Compute shrink point
                Matrix <T> XSH = XX ( 0, 0, 0, n - 1 ) + ( XX ( i, i, 0, n - 1 ) - XX ( 0, 0, 0, n - 1 ) ) * Sigma;

		//Reflection into the search space
                reflection (XSH, XMIN, XMAX);

                //Compute residuals
                function ( XSH, YSH, VSH, WSH );
		const T  fvs= norm(trans(VSH) * WSH * VSH);

                //Set submatrix
                XX ( XSH, i, 0 );
                VV2 ( i, 0 ) = fvs;
        }
}


template <typename T>
void SimplexMethod::reflection ( Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX)
{
	//Set elements of the simplex into the search space represented by the n-dimensional cuboid

	//Transpose matrix
	Matrix <T> XT = trans(X);
	
	//Call reflection
	NonLinearLeastSquares::reflection(XT, trans(XMIN), trans(XMAX));

	//Assign transposed matrix
	X = trans(XT);
}

#endif
