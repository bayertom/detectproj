// Description: Abstract class with graph definition

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


#ifndef Graph_H
#define Graph_H

#include <iostream>
#include <ostream>

#include "libalgo/source/structures/list/IndexLists.h"


//Abstract class representing graph
class Graph
{
        protected:
                TIndexList vertices;
                bool graph_type;

        public:
                Graph ( const bool graph_type_ = 0 ) : graph_type ( graph_type_ ) {}
                Graph ( TIndexList & vertices_, bool graph_type_ ) : vertices ( vertices_ ), graph_type ( graph_type_ ) {}
                virtual ~Graph() = 0;
                TIndexList::iterator begin () { return vertices.begin(); }
                TIndexList::iterator end() { return vertices.end(); }

        public:
                TIndexList const & getVertices() const {return vertices;}
                TIndexList & getVertices() {return  vertices;}
                bool getGraphType() const {return graph_type;}
                void setGraphType ( const bool & graph_type_ ) {graph_type = graph_type_;}

        public:
                //virtual void print ( std::ostream * output = & std::cout ) const = 0;
                //friend void operator << ( std::ostream & output, const Graph & g ) {g.print ( &output );}
                //virtual void clear() = 0;
};

inline Graph::~Graph() {}

#endif
