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


#ifndef DDT2D_HPP
#define DDT2D_HPP

#include <vector>
#include <list>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/line/HalfEdge.h"

#include "libalgo/source/algorithms/convexquadrilateral/ConvexQuadrilateral.h"
#include "libalgo/source/algorithms/swappingcriteria/SwappingCriteria.h"
#include "DT2D.h"

inline void DDT2D:: setSwapCriterion ( unsigned short swap_criterion_selected )
{
        //Set swap criterion
        switch ( swap_criterion_selected )
        {
                case 0:
                        //Angle between normals (ABN)
                        pcriterion = &SwappingCriteria::getAbn <double>;
                        break;

                case 1:
                        //Smooth of contours (SCO)
                        pcriterion = &SwappingCriteria::getSco <double>;
                        break;

                default:
                        //Angle between normals (ABN)
                        pcriterion = &SwappingCriteria::getAbn <double>;

        }
}


template <typename T>
void DDT2D::DDTLOP ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges, unsigned short swap_criterion_selected, const bool print_message, const bool print_exception, std::ostream * output )
{
        // Data Depending triangulation using selected local criterion
        bool swap_exist = true;

        //Set iterations
        unsigned iterations = 0;

        //Create Delaunay triangulation
        DT2D::DT ( nl, half_edges, print_message );

        //Get number of HalfEdges
        const unsigned int n = half_edges->size();

        //Set local swap criterion
        setSwapCriterion ( swap_criterion_selected );

        //Global cost before swapping
        const T global_cost_old = globalCostFunction ( half_edges );
        T global_cost = global_cost_old;

        //Print info
        if ( print_message )
        {
                *output << "> Starting DDTLOP... " ;
        }

        //Initialize counters
        unsigned int counter = MAX_INT, counter_old = MAX_INT;

        try
        {
                //Run until swap exists or decrease number of swaps between two loops
                do
                {
                        //We suppose ordered set of triangles, no swap will be required
                        swap_exist = false;

                        //Remember old counter
                        counter_old = counter;

                        //Assign new counter value
                        counter = 0;

                        //Loop all edges
                        for ( unsigned int i = 0; i < n; i++ )
                        {

                                //Take half edge
                                HalfEdge <T> *e = ( *half_edges ) [i];

                                //Use simplex flag to eliminate T processing of the edge
                                if ( !e->isSimplexEdge() )
                                {
                                        //Does twin edge exist?
                                        if ( e->getTwinEdge() )
                                        {
                                                // Test of convexity for quadrilateral
                                                HalfEdge <T> *e12 = e->getNextEdge();
                                                HalfEdge <T> *e13 = e12->getNextEdge();
                                                HalfEdge <T> *e21 = e->getTwinEdge();
                                                HalfEdge <T> *e22 = e21->getNextEdge();
                                                HalfEdge <T> *e23 = e22->getNextEdge();

                                                // Get nodes, counterclockwise set of nodes
                                                const Node3DCartesian <T> *p1 = e->getPoint();
                                                const Node3DCartesian <T> *p2 = e23->getPoint();
                                                const Node3DCartesian <T> *p3 = e12->getPoint();
                                                const Node3DCartesian <T> *p4 = e13->getPoint();

                                                //Is convex (non convex can not be swapped)
                                                if ( ConvexQuadrilateral::isStrictlyConvex ( p1, p2, p3, p4 ) == 1 )
                                                {

                                                        //Set twin edge to be processed
                                                        e->getTwinEdge()->setEdgeAsSimplex ( true );

                                                        //Get first triangle
                                                        e12 = e->getNextEdge();
                                                        e13 = e12->getNextEdge();

                                                        //Get coordinates (cast to parent using static_cast)
                                                        const Point3DCartesian <T> *p11 = e->getPoint();
                                                        const Point3DCartesian <T> *p12 = e12->getPoint();
                                                        const Point3DCartesian <T> *p13 = e13->getPoint();

                                                        //Get second triangle
                                                        e21 = e->getTwinEdge();
                                                        e22 = e21->getNextEdge();
                                                        e23 = e22->getNextEdge();

                                                        //Get coordinates
                                                        const Point3DCartesian <T> *p21 = p12;
                                                        const Point3DCartesian <T> *p22 = p11;
                                                        const Point3DCartesian <T> *p23 = e23->getPoint();

                                                        //Compute local criterion
                                                        const T c1 = ( *pcriterion ) ( p11, p12, p13, p21, p22, p23 );

                                                        //Calculation local criterion from swapped diagonal
                                                        const T c2 = ( *pcriterion ) ( p13, p11, p23, p12, p13, p23 );

                                                        //Swap diagonal
                                                        if ( c2 < c1 )
                                                        {
                                                                //Swap exists
                                                                swap_exist = true;

                                                                //Swap diagonal
                                                                DT2D::swapDiagonal ( e, e12, e13, e21, e22, e23 );

                                                                counter++;
                                                        }

                                                        //Set both edges to be processed
                                                        e->setEdgeAsSimplex ( true );
                                                        e21->setEdgeAsSimplex ( true );
                                                }
                                        }
                                }
                        }

                        //Set all edges as unused
                        for ( unsigned int i = 0; i < n; i++ )
                        {


                                //Take half edge
                                HalfEdge <T> *e = ( *half_edges ) [i];

                                //Set edge not to be a simplex
                                e->setEdgeAsSimplex ( false );
                        }

                        //Number of iterations
                        iterations ++;
                }
                while ( swap_exist );

                //Global cost after swapping
                global_cost = globalCostFunction ( half_edges );

                //Compute change ratio
                T ddt_ratio = 0;

                if ( global_cost  != global_cost_old )
                {
                        ddt_ratio = 100 * ( global_cost - global_cost_old ) / global_cost_old;
                }

                //Print info
                if ( print_message )
                {
                        *output << "Completed..." << std::endl;
                        *output << "> DT global cost old: " << global_cost_old << "  DT global cost new: " << global_cost << ", ratio: " << ddt_ratio << "%" << std::endl;
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                if ( print_exception )
                {
                        error.printException ( output );
                }

                //Delete lists
                half_edges->clear();

                *output << "DDT construction cancelled..." << std::endl;

                throw;
        }
}


template <typename T>
void DDT2D::DDTSimulatedAnnealing ( Container <Node3DCartesian <T> *> &nl, Container <HalfEdge <T> *> &half_edges, unsigned short swap_criterion_selected, const bool print_message, const bool print_exception, std::ostream * output )
{
        //DelaunayTriangulation2D using simulated annealing algorithm
        try
        {
                //Initial and minimal temperatures
                float t0 = 200;     		   					//Initial temperature
                const float tmin = 5; 							//Minimum temperature
                const float r = 0.95f;							//Annealing ratio r = (0.9, 0.95)

                //Do LOP as the first iteration
                DDTLOP ( nl, half_edges, swap_criterion_selected, print_message );

                const unsigned int glimit_multiplier = 5;  				//Good swaps glimit ratio, multiplier = 5 * n													//Number of good swaps glimit= 5 * n
                const unsigned int nlimit_multiplier = 10;				//Multiplier of iterations for each temperature  = 10 * n
                float k = 0;								//Power of the Boltzmann equation

                //Iterations and global cost
                unsigned int iterations = 0;						//Counter of iterations
                T global_cost_old = globalCostFunction ( half_edges );
                T global_cost = global_cost_old;

                //Automatic t0 temperature set
                bool automatic_t0_set = true;						//Automatic setting of the temperature

                //Total half edges
                const unsigned int n = half_edges.size();

                //Set new glimit and nlimit
                const unsigned int nlimit = nlimit_multiplier * n;
                const unsigned int glimit = glimit_multiplier * n;

                //Set initial temperature
                if ( automatic_t0_set )
                {

                        //Get initial temperature from analysis of cost differences
                        t0 = getInitialTemperature ( half_edges );
                }

                //Assign tk to t0
                float tk = t0;

                //Compute number of temperature states
                unsigned int total_tk_states = ( unsigned int ) ( 1 + ( log ( tmin ) - log ( t0 ) ) / ( log ( 0.95 ) ) );

                //Initialize random number generator
                srand ( ( unsigned ) time ( 0 ) );

                //Print info
                if ( print_message )
                {
                        *output << "> Starting DDT, Simulated Annealing, please wait..." << std::endl;
                }

                //Do until the temperature is frozen (tk > MIN_TEMPERATURE)
                do
                {
                        //Good swap indicator
                        unsigned int good_swap = 0;

                        //Total iterations for temperature tk
                        unsigned int iterations_tk = 0;

                        // Set new stat for temperature tk (MAX_ITERATIONS or GODD SWAP found)
                        do
                        {
                                //Set state
                                setNewState ( tk, good_swap, n, global_cost, half_edges );

                                //Increase iteration for the temperature tk (compare with n limit for each tk)
                                iterations_tk ++;

                                //Global iterations count (only for information)
                                iterations ++;
                        }

                        // Repeat if good_swaps < glimit
                        while ( good_swap < glimit );

                        //Set good swap to 0
                        good_swap = 0;

                        //Set new decreased temperature
                        tk = pow ( r, ( int ) k ) * t0;

                        //Increase power of the Boltzmann equation
                        k++;

                        // Update counter
                        if ( ( ( unsigned short ) ( 100.0 * ( k / total_tk_states ) ) ) % 1 == 0 )
                        {
                                *output << ".";
                        }
                }
                while ( tk > tmin );

                //Global cost after swapping
                const T global_cost_new = globalCostFunction ( half_edges );

                //Compute change ratio
                T ddt_ratio = 0;

                if ( global_cost  != global_cost_old )
                {
                        ddt_ratio = 100 * ( global_cost - global_cost_old ) / global_cost_old;
                }

                //Print info
                if ( print_message )
                {
                        *output << std::endl << "Finished..." << "Total iterations: " << iterations << std::endl;
                        *output << "> DT global cost old: " << global_cost_old << "  DT global cost new: " << global_cost << ", ratio: " << ddt_ratio << "%" << std::endl;
                }
        }

        //Throw exception
        catch ( Exception & error )
        {
                if ( print_exception )
                {
                        error.printException ( output );
                }

                //Delete lists
                half_edges->clear();

                *output << "DDT canceled..." << std::endl;

                throw;
        }
}


template <typename T>
void DDT2D::setNewState ( float tk, unsigned int & good_swap, unsigned int n, T & global_cost, Container <HalfEdge <T> *> &half_edges )
{
        //Set new state or recover old state
        unsigned int i = int ( n * rand() / ( RAND_MAX + 1.0 ) );

        //i can not be grater than n
        if ( i >= n )
        {
                i = n - 1;
        }

        //Get random edge
        HalfEdge <T> *e = ( *half_edges ) [i];

        //Cost difference
        T cost_difference = 0.0001 * MAX_FLOAT;

        //Compute cost difference (convex quadrilateral), else set difference to MAX_FLOAT
        getCostDifference ( e, cost_difference, 0, 1 );

        //Set new state, perform swap of the edge
        if ( cost_difference < 0 )
        {

                //First triangle
                HalfEdge <T> *e12 = e->getNextEdge();
                HalfEdge <T> *e13 = e12->getNextEdge();

                //Get the second triangle T2
                HalfEdge <T> *e21 = e->getTwinEdge();
                HalfEdge <T> *e22 = e21->getNextEdge();
                HalfEdge <T> *e23 = e22->getNextEdge();

                //Swap diagonal
                DT2D::swapDiagonal ( e, e12, e13, e21, e22, e23 );

                //Increment good swap
                good_swap ++;

                //Assign old global cost
                global_cost += cost_difference;
        }

        //Decide if set a new state or old state
        else
        {

                //Generate random number (0, 1)
                const float theta = ( T ) rand() / ( T ) RAND_MAX;

                //Compare with Boltzmann function and decide about new state
                if ( theta < exp ( - cost_difference / tk ) )
                {

                        //Set new state (acceptable increasing of the cost), perform the swap

                        //First triangle
                        HalfEdge <T> *e12 = e->getNextEdge();
                        HalfEdge <T> *e13 = e12->getNextEdge();

                        //Get the second triangle T2
                        HalfEdge <T> *e21 = e->getTwinEdge();
                        HalfEdge <T> *e22 = e21->getNextEdge();
                        HalfEdge <T> *e23 = e22->getNextEdge();

                        //Swap diagonal
                        DT2D::swapDiagonal ( e, e12, e13, e21, e22, e23 );

                        //Assign old global cost
                        global_cost += cost_difference;
                }
        }
}


inline bool DDT2D::endIteration ( unsigned int & index, unsigned int & good_swaps, unsigned int glimit, unsigned int nlimit )
{
        //End/continue iteration on the selected temperature tk
        if ( index < nlimit )
        {

                //Continue in calculation
                return true;
        }

        //Some good swaps were found?
        if ( good_swaps > glimit )
        {

                //Run again on the same temperature
                good_swaps = 0;

                //Update index
                index = 0;

                //Run another nlimit iterations
                return true;
        }

        //No good swap, change the temperature
        return false;
}


template <typename T>
float DDT2D::getInitialTemperature ( Container <HalfEdge <T> *> &half_edges )
{
        // Compute cost of the triangulation using selected criterion
        const unsigned int n = half_edges.size();
        T max_difference = 0; ;

        //Loop all edges
        for ( unsigned int i = 0; i < n; i++ )
        {

                //Take half edge
                HalfEdge <T> *e = half_edges [i];

                if ( !e->isSimplexEdge() )
                {
                        //Get next edges
                        const HalfEdge <T> *e12 = e->getNextEdge();
                        const HalfEdge <T> *e13 = e12->getNextEdge();

                        //Does a twin edge exist?
                        if ( e->getTwinEdge() )
                        {
                                //Get coordinates of the first triangle (cast to parent using static_cast)
                                const Point3DCartesian <T> *p11 = e->getPoint();
                                const Point3DCartesian <T> *p12 = e12->getPoint();
                                const Point3DCartesian <T> *p13 = e13->getPoint();

                                //Get second triangle
                                HalfEdge <T> *e21 = e->getTwinEdge();
                                HalfEdge <T> *e22 = e21->getNextEdge();
                                HalfEdge <T> *e23 = e22->getNextEdge();

                                //Get coordinates
                                const Point3DCartesian <T> *p21 = p12;
                                const Point3DCartesian <T> *p22 = p11;
                                const Point3DCartesian <T> *p23 = e23->getPoint();

                                //Only strictly convex quadrilaterals
                                if ( ConvexQuadrilateral::isStrictlyConvex ( p11,  p23, p12,  p13 ) == 1 )
                                {
                                        //Compute local criterion  for adjacent triangles
                                        const T fi1 = ( *pcriterion ) ( p11, p12, p13, p21, p22, p23 );

                                        //Compute local criterion  for swapped triangles
                                        const T fi2 = ( *pcriterion ) ( p13, p23, p12, p23, p13, p11 );

                                        //Find max difference
                                        if ( fabs ( fi2 - fi1 ) > max_difference )
                                        {
                                                max_difference = fabs ( fi2 - fi1 );
                                        }
                                }

                                //Set both edges to be processed
                                e->setEdgeAsSimplex ( true );
                                e21->setEdgeAsSimplex ( true );
                        }
                }
        }

        //Reset attribute
        for ( unsigned int i = 0; i < n; i++ )
        {
                ( *half_edges ) [i]->setEdgeAsSimplex ( false );
        }

        //Return result
        return ( float ) ( 2 * max_difference );
}


template <typename T>
T DDT2D::globalCostFunction ( Container <HalfEdge <T> *> &half_edges )
{
        // Compute cost of the triangulation using selected criterion
        const unsigned int n = half_edges.size();
        T cost = 0;

        //Loop all edges
        for ( unsigned int i = 0; i < n; i++ )
        {

                //Take half edge
                HalfEdge <T> *e = half_edges [i];

                //Use simplex indentificator to eliminate T processing of the edge
                if ( !e->isSimplexEdge() )
                {
                        //Get the first triangle
                        const HalfEdge <T> *e12 = e->getNextEdge();
                        const HalfEdge <T> *e13 = e12->getNextEdge();

                        //Does a twin edge exist?
                        if ( e->getTwinEdge() )
                        {
                                //Get coordinates of the first triangle
                                const Point3DCartesian <T> *p11 = e->getPoint();
                                const Point3DCartesian <T> *p12 = e12->getPoint() ;
                                const Point3DCartesian <T> *p13 =  e13->getPoint() ;

                                //Get second triangle
                                HalfEdge <T> *e21 = e->getTwinEdge();
                                const HalfEdge <T> *e22 = e21->getNextEdge();
                                const HalfEdge <T> *e23 = e22->getNextEdge();

                                //Get coordinates
                                const Point3DCartesian <T> *p21 = p12;
                                const Point3DCartesian <T> *p22 = p11;
                                const Point3DCartesian <T> *p23 = e23->getPoint();

                                //Compute local criterion
                                T fi = ( *pcriterion ) ( p11, p12, p13, p21, p22, p23 );

                                //Global cost
                                cost += fi;

                                //Set both edges to be processed
                                e->setEdgeAsSimplex ( true );
                                e21->setEdgeAsSimplex ( true );
                        }
                }
        }

        //Reset attribute
        for ( unsigned int i = 0; i < n; i++ )
        {
                ( *half_edges ) [i]->setEdgeAsSimplex ( false );
        }

        //Return global cost of the triangulation
        return cost;
}


template <typename T>
void DDT2D::getCostDifference ( HalfEdge <T> *e, T & cost_difference,  const bool simplex, const bool convex )
{
        /*Compute cost of the triangles T1, T2, T3, T4, T5, T6 adjacent to half edge and swapped variant
        It takes into account, if we want to process all edges or all triangles in DT
        It takes into account if we want to process only convex quadrilaterals */

        T c12 = 0, c12s = 0, c13 = 0, c13s = 0, c14 = 0, c14s = 0, c25 = 0, c25s = 0, c26 = 0, c26s = 0;
        T local_cost_swap = 0, local_cost = 0;

        //Use simplex indentificator to eliminate T processing of the edge
        if ( ( !simplex ) || !e->isSimplexEdge() )
        {
                //Does twin edge of this edge exist?
                if ( e->getTwinEdge() )
                {
                        const HalfEdge <T> *e12 = e->getNextEdge();
                        const HalfEdge <T> *e13 = e12->getNextEdge();

                        //Vertices of the first triangle
                        const Point3DCartesian <T> *p11 =  e->getPoint();
                        const Point3DCartesian <T> *p12 =  e12->getPoint();
                        const Point3DCartesian <T> *p13 =  e13->getPoint();

                        //Get the second triangle T2
                        HalfEdge <T> *e21 = e->getTwinEdge();
                        const HalfEdge <T> *e22 = e21->getNextEdge();
                        const HalfEdge <T> *e23 = e22->getNextEdge();

                        //Last vertex of the second triangle
                        const Point3DCartesian <T> *p23 =  e23->getPoint();

                        //Is a quadrilateral strictly convex or non-convex (both can not be swapped)
                        if ( ( !convex ) || ( ConvexQuadrilateral::isStrictlyConvex ( p11, p23, p12, p13 ) == 1 ) )
                        {
                                //Compute fi_12 (T1 and T2)
                                c12 = ( *pcriterion ) ( p11, p12, p13, p12, p11, p23 );
                                c12s = ( *pcriterion ) ( p12, p13, p23, p11, p23, p13 );

                                //Get the third triangle T3, adjacent to the first triangle T1
                                if ( e12->getTwinEdge() )
                                {

                                        //Triangle T3 exists
                                        const HalfEdge <T> *e31 = e12->getTwinEdge();
                                        const HalfEdge <T> *e32 = e31->getNextEdge();
                                        const HalfEdge <T> *e33 = e32->getNextEdge();

                                        //Last vertex of the third triangle
                                        const Point3DCartesian <T> *p33 = e33->getPoint();

                                        //Compute fi_13 (T1 and T3)
                                        c13 = ( *pcriterion ) ( p11, p12, p13, p13, p12, p33 );
                                        c13s = ( *pcriterion ) ( p13, p23, p12, p13, p12, p33 );
                                }

                                //Get the fourth triangle T4, adjacent to the first triangle T1
                                if ( e13->getTwinEdge() )
                                {
                                        //Triangle T4 exists
                                        const HalfEdge <T> *e41 = e13->getTwinEdge();
                                        const HalfEdge <T> *e42 = e41->getNextEdge();
                                        const HalfEdge <T> *e43 = e42->getNextEdge();

                                        //Last vertex of the fourth triangle
                                        const Point3DCartesian <T> *p43 = e43->getPoint();

                                        //Compute fi_14 (T1 and T4)
                                        c14 = ( *pcriterion ) ( p11, p12, p13, p11, p13, p43 );
                                        c14s = ( *pcriterion ) ( p13, p11, p23, p11, p13, p43 );
                                }

                                //Get the fifth triangle T5, adjacent to the second triangle T2
                                if ( e22->getTwinEdge() )
                                {

                                        //Triangle T5 exists
                                        const HalfEdge <T> *e51 = e22->getTwinEdge();
                                        const HalfEdge <T> *e52 = e51->getNextEdge();
                                        const HalfEdge <T> *e53 = e52->getNextEdge();

                                        //Last vertex of the fifth triangle
                                        const Point3DCartesian <T> *p53 = e53->getPoint();

                                        //Compute fi_25 (T2 and T5)
                                        c25 = ( *pcriterion ) ( p11, p23, p12, p11, p53, p23 );
                                        c25s = ( *pcriterion ) ( p13, p11, p23, p11, p53, p23 );
                                }

                                //Get the sixth triangle T6, adjacent to the second triangle T2
                                if ( e23->getTwinEdge() )
                                {

                                        //Triangle T6 exists
                                        const HalfEdge <T> *e61 = e23->getTwinEdge();
                                        const HalfEdge <T> *e62 = e61->getNextEdge();
                                        const HalfEdge <T> *e63 = e62->getNextEdge();

                                        //Last vertex of the sixth triangle
                                        const Point3DCartesian <T> *p63 =  e63->getPoint();

                                        //Compute fi_26 (T2 and T6)
                                        c26 = ( *pcriterion ) ( p11, p23, p12, p12, p23, p63 );
                                        c26s = ( *pcriterion ) ( p13, p23, p12, p12, p23, p63 );
                                }

                                //Compute cost criterion before swap
                                local_cost = c12 + c13  + c14  + c25  + c26 ;

                                //Compute cost criterion after swap
                                local_cost_swap = c12s + c13s + c14s + c25s + c26s;

                                //Set both edges to be processed
                                if ( simplex )
                                {
                                        e->setEdgeAsSimplex ( true );
                                        e21->setEdgeAsSimplex ( true );
                                }

                                //Cost difference
                                cost_difference = local_cost_swap - local_cost;
                        }
                }
        }
}

#endif
