// Description: Sort meridians by longitude

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

#ifndef sortMeridiansByLon_H
#define sortMeridiansByLon_H

//Sort meridians by longitude (Generic sorter)
template <typename Meridian>
class sortMeridiansByLon
{

        public:

                bool operator() ( const Meridian & m1, const Meridian & m2 ) const
                {
                        return m1.getLon() < m2.getLon();
                }
};


//Sort meridians by longitude (Partial specialization for Meridian *)
template <typename Meridian>
class sortMeridiansByLon <Meridian *>
{
        public:

                bool operator() ( const Meridian * m1, const Meridian * m2 ) const
                {
                        return m1->getLon() < m2->getLon();
                }
};


#endif
