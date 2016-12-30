// Description: Find global minimum using the differential evolution algorithm

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


#ifndef DifferentialEvolution_H
#define DifferentialEvolution_H

#include <ostream>
#include <iostream>

#include "libalgo/source/types/TMutationStrategy.h"
#include "libalgo/source/types/TAdaptiveControl.h"

//Include  C++98/C++11 version of the library
#if CPP11_SUPPORT == 0 
	#include "libalgo/source/structures/matrix/Matrix.h"
#else
	#include "libalgo/source/structures/matrix2/Matrix.h"
#endif

//Example: Schwefel function
template <typename T>
class FSchwefel
{
        public:
                T operator () ( const Matrix <T> &arg )
                {
                        return -arg ( 0, 0 ) * sin ( sqrt ( fabs ( arg ( 0, 0 ) ) ) ) - arg ( 0, 1 ) * sin ( sqrt ( fabs ( arg ( 0, 1 ) ) ) );
                }
};


//Find global minimum using the differential evolution algorithm.
//The minimized function is defined by the functor, it is placed in operator () (const Matrix <T> &arg)
class DifferentialEvolution
{
        public:

		template <typename T, typename Function>
		static T diffEvolution(Function function, const unsigned int population_size, const T epsilon, const unsigned int max_iterations, Matrix <T> F, T CR, const TMutationStrategy & mutation_strategy, const TAdaptiveControl &adaptive_control, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &RES, Matrix <T> &XMIN, Matrix <T> &XMAX, Matrix <T> &XAVER, T &aver_res, T &fx_max, unsigned int &iterations, const bool add_x0 = false, std::ostream * output = &std::cout);


        private:

		template <typename T, typename Function>
		static void createInitialPopulation(Function function, const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &W, Matrix <T> &Y, Matrix <T> &RES, const unsigned int population_size, const unsigned int dim, Matrix <T> &P_A, Matrix <T> &P_V, const bool add_x0 = false);

                template <typename T>
		static void mutationStrategyDERand1(const Matrix <T> &P_A, const unsigned int & i, const unsigned int & population_size, const Matrix <T> & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDERand2(const Matrix <T> &P_A, const unsigned int & i, const unsigned int & population_size, const Matrix <T> & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDERandDir1(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T> & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDERandDir2(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T>  & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDERandBest1(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T>  & F, Matrix <T> &U);

                template <typename T>
		static void mutationStrategyDERandBest2(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T>  & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDERandBestDir1(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T>  & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategyDETargetToBest1(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const Matrix <T>  & F, Matrix <T> &U);

		template <typename T>
		static void mutationStrategySACP(const Matrix <T> &P_A, const Matrix <T> &P_V, const unsigned int & i, const unsigned int & population_size, const T &aver_res, Matrix <T> & F, T &CR, Matrix <T> &U);

                template <typename T>
                static void crossover ( const Matrix <T> &U,  const T & CR, const unsigned int dim, Matrix <T> &V );

		template <typename T>
		static void crossover(const Matrix <T> &U, const unsigned int dim, const Matrix <T> &AS, Matrix <T> &V);

                template <typename T>
                static void reflection (Matrix <T> &X, const Matrix <T> &XMIN, const Matrix <T> &XMAX);

};

#include "DifferentialEvolution.hpp"

#endif
