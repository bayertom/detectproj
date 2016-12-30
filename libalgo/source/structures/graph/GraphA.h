// Description: List representation of graph (stores the adjacency in std::map)

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


#ifndef GraphA_H
#define GraphA_H

#include <list>
#include <map>

#include "GraphEdgeA.h"
#include "Graph.h"


//Adjacency described using adjacency map and represented by std::list
template <typename T>
struct AdjacencyMap
{
        typedef std::map <int, std::list < GraphEdgeA <T> > > adjacency_map_t;
};


//List representation of graph (does not store the adjacency)
template <typename T>
class GraphA : public Graph
{
        protected:
                typename AdjacencyMap <T>::Type a_map;		//Adjacency map

        public:
                GraphA() : Graph (), a_map ( 0 ) {}

                GraphA ( TIndexList & vertices, typename AdjacencyMap <T>::Type & a_map_, const bool graph_type = 0 ) :
                        Graph ( vertices, graph_type ), a_map ( a_map_ ) {}
                virtual ~GraphA() {}

        public:
                typename AdjacencyMap <T>::Type const & getAMap () const {return a_map;}
                typename AdjacencyMap <T>::Type & getAMap () {return a_map;}
                virtual void print() const {};
                virtual void clear() {}
};

#endif
