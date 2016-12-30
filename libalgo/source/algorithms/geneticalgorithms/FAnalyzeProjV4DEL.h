// Description: Functor, compute matrix V of squares of residuals for cartometric analysis, determine lon0 (1D optimization)
// Method: Differential evolution, hybrid method, rotation determioned from similarity transformation

// Copyright (c) 2010 - 2015
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


#ifndef FAnalyzeProjV4DEL_H
#define FAnalyzeProjV4DEL_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"


//Forward declaration
template <typename T>
class Projection;


//Functor, compute matrix V of squares of residuals for cartometric analysis, determine lon0 (1D optimization)
template <typename T>
class FAnalyzeProjV4DEL
{
	private:

		Container <Node3DCartesian <T> *> &nl_test;
		Container <Point3DGeographic <T> *> &pl_reference;
		typename TMeridiansList <T> ::Type &meridians;
		typename TParallelsList <T> ::Type &parallels;
		const Container <Face <T> *> &faces_test;
		Projection <T> *proj;
		T &R;
		T &q1;
		T &q2;
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

		FAnalyzeProjV4DEL(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
			const Container <Face <T> *> &faces_test_, Projection <T> *proj_, T &R_est_, T &q1_, T &q2_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &res_evaluations_, const TMEstimatorsWeightFunction &me_function_, const T k_, Matrix <unsigned int> &I_, std::ostream * output_)
			: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R(R_est_), q1(q1_), q2(q2_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
			created_samples(created_samples_), res_evaluations (res_evaluations_), me_function(me_function_), k(k_), I(I_), output(output_) {}


		void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
		{

			//Compute parameters of the V matrix: residuals
			evaluateResidualsL(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res, created_samples, res_evaluations, me_function, k, I, output);
		}
};



template <typename T>
void evaluateResidualsL(const Matrix <T> &XL, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, T &R, T &q1, T &q2, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples,
	unsigned int &res_evaluation, const TMEstimatorsWeightFunction &me_function, const T k, Matrix <unsigned int> &I, std::ostream * output)
{
	//Simple wrapper calling method from FAnalyzeProjV4
	//const unsigned int m = nl_test.size(), n = X.rows();

	//Create matrix of deter
	Matrix <T> X(1, 5);
	Point3DGeographic <T> cart_pole;
	cart_pole = proj->getCartPole();

	X(0, 0) = cart_pole.getLat();
	X(0, 1) = cart_pole.getLon();
	X(0, 2) = proj->getLat0();
	X(0, 3) = XL(0, 0);
	X(0, 4) = proj->getC();

	//Call function from FAnalyzeProjV4
	evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);
}

#endif