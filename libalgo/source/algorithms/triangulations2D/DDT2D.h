// Description: 2D data depending triangulations using LOP and Simulated annealing

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


#ifndef DDT2D_H
#define DDT2D_H

#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class HalfEdge;

template <typename T>
class Point3DCartesian;


//Data depending triangulation
class DDT2D
{
        private:
                static double ( * pcriterion ) ( const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> * );  //Pointer to local swap criteria

        public:
                template <typename T>
                static void DDTLOP ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges, unsigned short swap_criterion_selected, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static void DDTSimulatedAnnealing ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges, unsigned short swap_criterion_selected, const bool print_message = false, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static float getInitialTemperature ( Container <HalfEdge <T> *> &half_edges );


        private:
                static void setSwapCriterion ( unsigned short swap_criterion_selected );

                template <typename T>
                static void setNewState ( float tk, unsigned int & good_swap, unsigned int n, T & global_cost, Container <HalfEdge <T> *> &half_edges );

                template <typename T>
                static T globalCostFunction ( Container <HalfEdge <T> *> &half_edges );

                template <typename T>
                static void getCostDifference ( HalfEdge <T> *e, T & cost_difference, const bool simplex, const bool convex );

                static bool endIteration ( unsigned int & index, unsigned int & good_swaps, unsigned int glimit_multiplier, unsigned int nlimit_multiplier );
};

#include "DDT2D.hpp"

#endif
