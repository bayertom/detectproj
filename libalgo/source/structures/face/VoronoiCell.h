// Description: Structure representing Voronoi cell

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

#ifndef VoronoiCell_H
#define VoronoiCell_H

#include <iostream>
#include <ostream>

//Forward declarations
template <typename T>
class Node3DCartesian;

template <typename T>
class HalfEdge;

template <typename T>
class Face;

//Structure representing Voronoi cell using circular list
template <typename T>
class VoronoiCell: public Face <T>
{
        private:
                Node3DCartesian <T>  *generator;
                bool bounded;

        public:

                VoronoiCell ( Node3DCartesian <T> *generator_, HalfEdge <T> *edge, bool bounded_ ) : Face <T> ( edge ), generator ( generator_ ), bounded ( bounded_ ) {}
                virtual ~VoronoiCell();

        public:
                Node3DCartesian <T> *getGenerator() const {return generator;}
                bool getBounded() const {return bounded;}

                void setGenerator ( Node3DCartesian <T> *generator_ ) {generator = generator_;}
                void setBounded ( const bool bounded_ ) {bounded = bounded_;}

        public:
                virtual void removeAdjacency();
                virtual void print ( std::ostream * output = &std::cout ) const;
};

#include "VoronoiCell.hpp"

#endif

