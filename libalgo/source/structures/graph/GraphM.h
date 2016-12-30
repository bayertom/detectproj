// Description: Matrix representation of graph

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


#ifndef GraphM_H
#define GraphM_H

#include "Graph.h"
#include "libalgo/source/structures/matrix/Matrix.h"


//Matrix representation of graph
template <typename T>
class GraphM : public Graph
{
        protected:
                Matrix <unsigned short> V;	//Adjacency matrix
                Matrix <T> W;			//Weight matrix

        public:
                GraphM ( const unsigned int n_vertices, const bool graph_type = 0 ) :
                        Graph ( graph_type ), V ( n_vertices, n_vertices ), W ( n_vertices, n_vertices ) {}

                GraphM ( TIndexList & vertices, const Matrix <unsigned short> &V_, const Matrix <T> W_,
                         const bool graph_type = 0 ) : Graph ( vertices, graph_type ), V ( V_ ), W ( W_ ) {}

                virtual ~GraphM() {}

        public:
                Matrix <unsigned short> const & getV () const {return V;}
                Matrix <unsigned short> & getV () {return V;}
                Matrix <T> const & getW () const {return W;}
                Matrix <T> & getW () {return W;}

                void setV ( const Matrix <unsigned short> &  V_ ) {V = V_;}
                void setW ( const Matrix <T> &  W_ ) {W = W_;}

        public:
                virtual void print() const {};
                virtual void clear() {}
};

#endif
