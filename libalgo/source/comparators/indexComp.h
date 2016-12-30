// Description: Comparator, returns indices of sorted elements
// Used to store a permutation of indices of sorted elements

// Copyright (c) 2010 - 2016
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


#ifndef indexComp_H
#define indexComp_H

#include <algorithm>


//Comparator, returns indices of sorted elements
template <class random_iterator>
class indexComp
{
        private:
                random_iterator const i_begin, i_end;

        public:
                indexComp ( random_iterator i_begin_, random_iterator i_end_ )
                        : i_begin ( i_begin_ ), i_end ( i_end_ ) {}

                bool operator () ( const unsigned int i1, const unsigned int i2 ) const
                {
                        return * ( i_begin + i1 ) < * ( i_begin + i2 );
                }

};

#endif
