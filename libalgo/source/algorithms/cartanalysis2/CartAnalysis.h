// Description: Performs cartometric analysis (i.e. estimation of the cartographic projection)

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

#ifndef CartAnalysis_H
#define CartAnalysis_H

#include <list>
#include <vector>
#include <map>
#include <memory>

#include "libalgo/source/types/TVector.h"
#include "libalgo/source/types/TListS.h"
#include "libalgo/source/types/TResult.h"
#include "libalgo/source/types/TResults.h"
#include "libalgo/source/types/TAnalysisMethod.h"

#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/projection2/Projection.h"
#include "libalgo/source/structures/matrix2/Matrix.h"


//Perform cartometric analysis
class CartAnalysis
{
	public:
		template <typename T>
		static void analyzeProjection(const TVector<Point3DCartesian<T> > &test_points, const TVector <Point3DGeographic<T> > &reference_points,
			TListS <Projection<T>> &projections, TResults <T> &results, const TAnalysisMethod &method = NLSM7, std::ostream &output = std::cout);

		template <typename T>
		static void printResults(TResults <T> &results, std::ostream &output = std::cout);
	
	private:
		
		template <typename T>
		static Matrix <T> X0M7(const unsigned int m, const std::shared_ptr<Projection <T> > proj, const TAnalysisMethod &method);

		template <typename T>
		static Matrix <T> X0M8(const unsigned int m, const std::shared_ptr<Projection <T> > proj, const TAnalysisMethod &method);

		template <typename T>
		static Matrix <T> AM7(const unsigned int m, const std::shared_ptr <Projection <T> > proj, const T &scale);

		template <typename T>
		static Matrix <T> BM7(const unsigned int m, const std::shared_ptr <Projection <T> > proj, const T &scale);

		template <typename T>
		static Matrix <T> AM8(const unsigned int m, const std::shared_ptr <Projection <T> > proj);

		template <typename T>
		static Matrix <T> BM8(const unsigned int m, const std::shared_ptr <Projection <T> > proj);

};


#include "CartAnalysis.hpp"

#endif