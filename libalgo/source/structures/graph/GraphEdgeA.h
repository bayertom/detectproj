// Description: Class representing graph edge e = (end, weight) used with adjacent map

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


#ifndef GraphEdgeA_H
#define GraphEdgeA_H

//Graph edge used with adjacent map
template <typename T>
class GraphEdgeA
{
        protected:

                unsigned int end;	//End point of the edge
                T weight;		//Edge weight

        public:

                GraphEdgeA() : end ( 0 ), weight ( 0 ) {}
                GraphEdgeA ( const unsigned int end_, T weight_ ) : end ( end_ ), weight ( weight_ ) {}

                bool operator < ( const GraphEdgeA <T> &g1 ) const
                {
                        return weight < g1.weight;
                }

        public:
                unsigned int getEndPoint() const {return end;}
                T getWeight () const {return weight;}
                void setEndPoint ( const unsigned int end_ ) {end = end_;}
                void setWeight ( const unsigned int weight_ ) {weight = weight_;}
};

#endif
