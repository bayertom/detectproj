// Description: //List representation of graph (does not store the adjacency)

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


#ifndef GraphL_H
#define GraphL_H

#include <vector>

#include "Graph.h"
#include "GraphEdge.h"


//List of graph edges
template <typename T>
struct GraphEdges
{
        typedef std::vector < GraphEdge <T> > Type;
};


//List representation of graph (does not store the adjacency)
template <typename T>
class GraphL : public Graph
{
        protected:
                typename GraphEdges <T>::Type edges;		//List of graph edges

        public:
                GraphL() : Graph () {}

                GraphL ( TIndexList & vertices, typename GraphEdges <T>::Type & edges_, const bool graph_type = 0 ) :
                        Graph ( vertices, graph_type ), edges ( edges_ ) {}
                virtual ~GraphL() {}

        public:
                typename GraphEdges <T>::Type const & getEdges () const {return edges;}
                typename GraphEdges <T>::Type & getEdges () {return edges;}
                virtual void print() const {};
                virtual void clear() {}
};

#endif
