// Description: Sort all created cartographic samples by all ratios using geometric mean (not so sensitive as mean to outlying items)
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

#ifndef sortSamplesByAllRatios_H
#define sortSamplesByAllRatios_H

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"

//Forward declaration
template <class T>
class Sample;

//Sorter of the cartographic samples by all ratios using the geometric mean (not so sensitive as mean to outlying items)
template <typename T>
class sortSamplesByAllRatios
{
        private:
                typename TAnalysisParameters <T>::TAnalysisType analysis_type;

        public:
                sortSamplesByAllRatios ( const typename TAnalysisParameters <T>::TAnalysisType & analysis_type_ ) : analysis_type ( analysis_type_ ) {}

                bool operator() ( const Sample <T> &s1, const Sample <T> &s2 ) const
                {
                        //Compare average of positions relative to all criteria using geometric mean
                        return s1.getSampleCost ( analysis_type ) < s2.getSampleCost ( analysis_type );
                }
};

#endif
