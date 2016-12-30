// Description: Face overlay operations: union, difference, intersection, xor
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


#ifndef FaceOverlay_HPP
#define FaceOverlay_HPP

#include <iterator>
#include <algorithm>

#include "libalgo/source/structures/point/Node3DCartesian.h"

#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/algorithms/pointfaceposition/PointFacePosition.h"
#include "libalgo/source/algorithms/simpleface/SimpleFace.h"

#include "libalgo/source/const/Const.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename T>
void FaceOverlay::createOverlay ( const Face <T> *f1, const Face <T> *f2, const TOverlayType overlay_type, Container <Face <T> *> &result, Container <Node3DCartesian <T> * > &intersections, Container <HalfEdge <T> * > &hl, bool & f1_simple, bool & f2_simple )
{
        //Compute overlay of 2 faces using boolean operations
        Container <Node3DCartesian <T> *, NonDestructable > f1_nodes;
        Container <Node3DCartesian <T> *, NonDestructable > f2_nodes;

        //Check if both input faces are simple: At least one face is simple
        if ( ! ( f1_simple = SimpleFace::isSimpleFace ( f1 ) ) || ! ( f2_simple = SimpleFace::isSimpleFace ( f2 ) ) )
        {
                throw BadDataException ( "BadDataException: can not perform faces overlay, ", " input faces are not simple. " );
        }

        //Convert face to undestructable list of nodes
        f1->toNodesList ( f1_nodes );
        f2->toNodesList ( f2_nodes );

        //Add first node as the last node of the list (close the face)
        f1_nodes.push_back ( f1_nodes[0] );
        f2_nodes.push_back ( f2_nodes[0] );

        //Compute intersections
        typename TIntersectionsList <T>::Type il1 ( f1_nodes.size() ), il2 ( f2_nodes.size() );
        computeLineSegmentsIntersections ( f1_nodes, f2_nodes, il1, il2, intersections );

        //Compute weights of edges (inner or outer)
        TIndexList weights1 ( f1_nodes.size() ), weights2 ( f2_nodes.size() );
        setWeights ( f1_nodes, f2_nodes, weights1 );
        setWeights ( f2_nodes, f1_nodes, weights2 );

        //Create inner and outer segments for both polygons
        typename TSegmentsType <T>::Type f1_inner, f1_outer, f2_inner, f2_outer;
        splitSegmentsByWeights ( f1_nodes, weights1, f1_inner, f1_outer );
        splitSegmentsByWeights ( f2_nodes, weights2, f2_inner, f2_outer );

        //Intersection
        if ( overlay_type == Intersection )
        {
                createFaceOverlay ( f1_inner, f2_inner, Intersection, result, hl );
        }

        //Union
        else if ( overlay_type == Union )
        {
                createFaceOverlay ( f1_outer, f2_outer, Union, result, hl );
        }

        //Difference
        else
        {
                createFaceOverlay ( f1_outer, f2_inner, Difference,  result, hl );
        }
}


template <typename T>
void FaceOverlay::computeLineSegmentsIntersections ( Container <Node3DCartesian <T> *, NonDestructable >  &f1_nodes, Container <Node3DCartesian <T> *, NonDestructable > &f2_nodes, typename TIntersectionsList <T>::Type & il1, typename TIntersectionsList <T>::Type & il2, Container <Node3DCartesian <T> * > &intersections )
{
        //Compute intersection of each edge od the first face with second face
        for ( unsigned int i = 0; i < f1_nodes.size() - 1; i++ )
        {
                //Create temporary containers
                for ( unsigned int j = 0; j < f2_nodes.size() - 1; j++ )
                {
                        //Test intersections if d(p1,p2) && d(p3, p4) > POSITION_ROUND_ERROR
                        if ( EuclDistance::getEuclDistance2D ( f1_nodes [i]->getX(), f1_nodes [i]->getY(), f1_nodes [i + 1]->getX(), f1_nodes [i + 1]->getY() ) > POSITION_ROUND_ERROR &&
                                        EuclDistance::getEuclDistance2D ( f2_nodes [j]->getX(), f2_nodes [j]->getY(), f2_nodes [j + 1]->getX(), f2_nodes [j + 1]->getY() ) > POSITION_ROUND_ERROR )
                        {
                                //Find intersections
                                T x_int, y_int;

                                unsigned short t = LineLinePosition::get2LineSegmentsPosition ( f1_nodes [i], f1_nodes [i + 1], f2_nodes [j], f2_nodes [j + 1], x_int, y_int );

                                //Segments intersect in one endpoint or in a common point
                                if ( ( t == 1 ) || ( t == 2 ) )
                                {
                                        //Compute distances of the intersection to the start point
                                        const T d1_sq = ( ( f1_nodes [i]->getX() - x_int ) * ( f1_nodes [i]->getX() - x_int ) +
                                                          ( f1_nodes [i]->getY() - y_int ) * ( f1_nodes [i]->getY() - y_int ) );
                                        const T d2_sq = ( ( f2_nodes [j]->getX() - x_int ) * ( f2_nodes [j]->getX() - x_int ) +
                                                          ( f2_nodes [j]->getY() - y_int ) * ( f2_nodes [j]->getY() - y_int ) );

                                        //Compute distances of the intersection to the end point
                                        const T d3_sq = ( ( f1_nodes [i + 1]->getX() - x_int ) * ( f1_nodes [i + 1]->getX() - x_int ) +
                                                          ( f1_nodes [i + 1]->getY() - y_int ) * ( f1_nodes [i + 1]->getY() - y_int ) );
                                        const T d4_sq = ( ( f2_nodes [j + 1]->getX() - x_int ) * ( f2_nodes [j + 1]->getX() - x_int ) +
                                                          ( f2_nodes [j + 1]->getY() - y_int ) * ( f2_nodes [j + 1]->getY() - y_int ) );

                                        //Compute distances of the segment
                                        const T d5_sq = ( ( f1_nodes [i]->getX() - f1_nodes [i + 1]->getX() ) * ( f1_nodes [i]->getX() - f1_nodes [i + 1]->getX() ) +
                                                          ( f1_nodes [i]->getY() - f1_nodes [i + 1]->getY() ) * ( f1_nodes [i]->getY() - f1_nodes [i + 1]->getY() ) );
                                        const T d6_sq = ( ( f2_nodes [j]->getX() - f2_nodes [j + 1]->getX() ) * ( f2_nodes [j]->getX() - f2_nodes [j + 1]->getX() ) +
                                                          ( f2_nodes [j]->getY() - f2_nodes [j + 1]->getY() ) * ( f2_nodes [j]->getY() - f2_nodes [j + 1]->getY() ) );

                                        //Create new point only if intersection is different from end points of segments and lying between end points (additional test)
                                        const bool b1 = ( d1_sq > MIN_POSITION_DIFF * MIN_POSITION_DIFF ) && ( d1_sq + MIN_POSITION_DIFF * MIN_POSITION_DIFF < d5_sq );
                                        const bool b2 = ( d2_sq > MIN_POSITION_DIFF * MIN_POSITION_DIFF ) && ( d2_sq + MIN_POSITION_DIFF * MIN_POSITION_DIFF < d6_sq );
                                        const bool b3 = ( d3_sq > MIN_POSITION_DIFF * MIN_POSITION_DIFF ) && ( d3_sq + MIN_POSITION_DIFF * MIN_POSITION_DIFF < d5_sq );
                                        const bool b4 = ( d4_sq > MIN_POSITION_DIFF * MIN_POSITION_DIFF ) && ( d4_sq + MIN_POSITION_DIFF * MIN_POSITION_DIFF < d6_sq );

                                        if ( ( b1 && b3 ) || ( b2 && b4 ) )
                                        {
                                                //Create new node
                                                Node3DCartesian <T> *inters = new Node3DCartesian <T> ( x_int, y_int );

                                                //Add intersection to the first list
                                                if ( b1 && b3 )
                                                {
                                                        il1 [i].insert ( std::pair <T, Node3DCartesian <T> *> ( d1_sq, inters ) );
                                                }

                                                //Add intersection to the second list
                                                if ( b2 && b4 )
                                                {
                                                        il2 [j].insert ( std::pair <T, Node3DCartesian <T> *> ( d2_sq, inters ) );
                                                }

                                                //Add intersection point to the list
                                                intersections.push_back ( inters ) ;
                                        }
                                }
                        }
                }
        }

        //Add intersections of the first face with the second one to the list of vertices
        typename TItemsList <Node3DCartesian <T> *>::Type ::iterator i_nodes1 = f1_nodes.begin();

        for ( unsigned int i = 0; i_nodes1 != f1_nodes.end() && i < il1.size(); i_nodes1++, i++ )
        {
                //Process all intersections
                typename std::map <T, Node3DCartesian <T> *> ::iterator i_map = il1 [i].begin();

                for ( ; i_map != il1 [i].end(); i_map ++ )
                {
                        //Insert intersection between old vertices
                        i_nodes1 = f1_nodes.insert ( ++i_nodes1, i_map->second );
                }
        }

        //Add intersections of the second face with the first one to the list of vertices
        typename TItemsList <Node3DCartesian <T> *>::Type ::iterator i_nodes2 = f2_nodes.begin();

        for ( unsigned int i = 0; i_nodes2 != f2_nodes.end() && i < il2.size(); ++i_nodes2, i++ )
        {
                //Process all intersections
                typename std::map <T, Node3DCartesian <T> *> ::iterator i_map = il2 [i].begin();

                for ( ; i_map != il2 [i].end(); i_map ++ )
                {
                        //Insert intersection between old vertices
                        i_nodes2 = f2_nodes.insert ( ++i_nodes2, ( i_map->second ) );
                }
        }
}


template <typename T>
void FaceOverlay::setWeights ( const Container <Node3DCartesian <T> *, NonDestructable > &f1_nodes, Container <Node3DCartesian <T> *, NonDestructable > &f2_nodes, TIndexList & weights )
{
        //Set weights of edges
        const unsigned int n = f1_nodes.size() ;

        //Process all edges
        for ( unsigned int i = 1; i < n; i++ )
        {

                //Compute mid-point
                Node3DCartesian <T>  mid ( ( f1_nodes [i]->getX() + f1_nodes [i - 1]->getX() ) / 2.0,
                                           ( f1_nodes [i]->getY() + f1_nodes [i - 1]->getY() ) / 2.0 );

                //Additional test of the mid-point
                unsigned short t_add = PointFacePosition::getPointFacePosition ( &mid, f2_nodes );

                //Outer edge: mid point outside polygon
                if ( t_add == 1 )
                {
                        weights[i - 1] = 1;
                }

                //Otherwise inner edge: mod point inside polygon, on ege or identical with polygon vertex
                else
                {
                        weights[i - 1] = 0;
                }
        }

        //Add weight of the last point
        weights[n - 1] = weights[0];
}


template <typename T>
void FaceOverlay::splitSegmentsByWeights ( const Container <Node3DCartesian <T> *, NonDestructable > &f_nodes, const TIndexList &weights, typename TSegmentsType <T>::Type & inner,  typename TSegmentsType <T>::Type & outer )
{
        //Split segments according to weight to inner and outer
        int k1 = -1, k2 = -1;

        //First edge is inner
        if ( weights[0] == 0 )
        {
                inner.push_back ( std::vector <Node3DCartesian <T> *> () );
                inner[++k1].push_back ( f_nodes [0] );
        }

        //First edge is outer
        else
        {
                outer.push_back ( std::vector <Node3DCartesian <T> *> () );
                outer[++k2].push_back ( f_nodes [0] );
        }

        //Process other edges
        for ( unsigned int i = 1; i < weights.size(); i++ )
        {
                //Inner segments
                if ( weights[i] == 0 )
                {
                        //Edge has different weight as previous edge
                        if ( weights[i] != weights[i - 1] )
                        {
                                //Add to the actual list: end point of the last edge
                                outer[k2].push_back ( f_nodes [i] );

                                //Increment k1
                                k1++;

                                //Create new list
                                inner.push_back ( std::vector <Node3DCartesian <T> *> () );
                        }

                        //Add to the new list: start point of the first edge
                        inner[k1].push_back ( f_nodes [i] );
                }

                //Outer segments
                else
                {
                        //Edge has different weight as previous edge
                        if ( weights[i] != weights[i - 1] )
                        {
                                //Add to the actual list: end point of the last edge
                                inner[k1].push_back ( f_nodes [i] );

                                //Increment k2
                                k2++;

                                //Create new list
                                outer.push_back ( std::vector <Node3DCartesian <T> *> () );
                        }

                        //Add to the new list: start point of the first edge
                        outer[k2].push_back ( f_nodes [i] );
                }
        }

        //Test last point of the last inner segment with first point of the first inner segment: inner segment was splitted
        if ( k1 > 0 ) //At least one inner segment found
        {
                unsigned int inner_last_size = inner[k1].size();

                if ( * ( inner[k1][inner_last_size - 1 ] ) == * ( inner[0][0] ) )
                {
                        inner[0].insert ( inner[0].begin(), inner[k1].begin(), inner[k1].end() - 1 );
                        inner.pop_back();
                }
        }

        //Test last point of the last outer segment with first point of the first outer segment: outer segment was splitted
        if ( k2 > 0 )  //At least one outer segment found
        {
                unsigned int outer_last_size = outer[k2].size();

                if ( * ( outer[k2][outer_last_size - 1] ) == * ( outer[0][0] ) )
                {
                        outer[0].insert ( outer[0].begin(), outer[k2].begin(), outer[k2].end() - 1 );
                        outer.pop_back();
                }
        }
}


template <typename T>
void FaceOverlay::createFaceOverlay ( typename TSegmentsType <T>::Type & f1,  typename TSegmentsType <T>::Type & f2, const TOverlayType & overlay_type, Container <Face <T> *> &result, Container <HalfEdge <T> * > &hl )
{
        //Create overlay depending on weights of points
        unsigned int i = 0, j = 0, empty_loops = 0;
        bool first = true, second = false, t1 = false, t2 = false;
        T min_dist = MIN_POSITION_DIFF;
        Container <Node3DCartesian <T> *, NonDestructable > face_3;

        //First list is empty
        if ( f1.size() == 0 ) { return; }

        //Remember first point of the face and start point of the next segment
        Node3DCartesian <T> * face_start_point = f1[i].front();
        Node3DCartesian <T> * next_edge_start_point = f1[i].front();

        //Merge lists, continue until both list are empty
        while ( ( !f1.empty() || !f2.empty() ) )
        {
                //Perform tests of the possible start point to both list (must be non empty)
                t1 = ( !f1.empty() ) && ( ( f1[i].front()->getX() - next_edge_start_point->getX() ) * ( f1[i].front()->getX() - next_edge_start_point->getX() ) +
                                          ( f1[i].front()->getY() - next_edge_start_point->getY() ) * ( f1[i].front()->getY() - next_edge_start_point->getY() ) < min_dist * min_dist ) ;
                overlay_type != Difference ? t2 = ( !f2.empty() ) && ( ( f2[j].front()->getX() - next_edge_start_point->getX() ) * ( f2[j].front()->getX() - next_edge_start_point->getX() ) +
                                                  ( f2[j].front()->getY() - next_edge_start_point->getY() ) * ( f2[j].front()->getY() - next_edge_start_point->getY() ) < min_dist * min_dist ) :
                                                  t2 = ( !f2.empty() ) && ( ( f2[j].back()->getX() - next_edge_start_point->getX() ) * ( f2[j].back()->getX() - next_edge_start_point->getX() ) +
                                                                  ( f2[j].back()->getY() - next_edge_start_point->getY() ) * ( f2[j].back()->getY() - next_edge_start_point->getY() ) < min_dist * min_dist );

                //First point of the first segment is a start point
                if ( t1 )
                {
                        //Merge both lists
                        face_3.insert ( face_3.end(), f1[i].begin(), f1[i].end() - 1 );

                        //Remember start point of the next edge
                        next_edge_start_point = f1[i].back();

                        //Erase sequence
                        f1.erase ( f1.begin() + i );

                        //Set indicators
                        first = true; second = false;

                        //Correct index i
                        if ( i == f1.size() ) i = 0;

                        //Reset empty loops and set min_dist to min value
                        empty_loops = 0;
                        min_dist = MIN_POSITION_DIFF;
                }

                //First point of the second segment is a start point
                else if ( t2 )
                {
                        //Merge both list
                        overlay_type != Difference ? face_3.insert ( face_3.end(), f2[j].begin(), f2[j].end() - 1 ) :
                        face_3.insert ( face_3.end(), f2[j].rbegin(), f2[j].rend() - 1 );

                        //Remember start point of the next edge
                        overlay_type != Difference ? next_edge_start_point = f2[j].back() : next_edge_start_point = f2[j].front();

                        //Erase sequence
                        f2.erase ( f2.begin() + j );

                        //Set indicators
                        first = false; second = true;

                        //Correct index j
                        if ( j == f2.size() ) j = 0;

                        //Reset empty loops and set min_dist to min value
                        empty_loops = 0;
                        min_dist = MIN_POSITION_DIFF;
                }

                //We close the face: we found start point or both lists are empty (to avoid round errors, after adding last segment min_dist can not be increased )
                if ( ( next_edge_start_point->getX() - face_start_point ->getX() ) * ( next_edge_start_point->getX() - face_start_point ->getX() ) +
                                ( next_edge_start_point->getY() - face_start_point ->getY() ) * ( next_edge_start_point->getY() - face_start_point ->getY() ) < min_dist * min_dist ||
                                ( f1.empty() && f2.empty() ) )
                {
                        //If area of the polygon is not zero or there is enough points
                        if ( ( face_3.size() > 2 ) && ( FaceArea::getFaceArea ( &face_3 ) > AREA_ROUND_ERROR ) )
                        {
                                //Create new face
                                Face <T> *face_new = new Face <T> ( face_3, hl );

                                //Add face to the list
                                result.push_back ( face_new );
                        }

                        //Clear list of vertices of the resulted polygon
                        face_3.clear();

                        //First list is not empty, continue
                        if ( !f1.empty() )
                        {
                                //Remeber start point and edge
                                face_start_point = f1[i].front();

                                //Remember start point of the next edge
                                next_edge_start_point = f1[i] .front();
                        }

                        //Second list is not empty (but first is empty), continue
                        else if ( !f2.empty() )
                        {
                                //Remember start point and edge
                                face_start_point = f2[j].front();

                                //Remember start point of the next edge
                                next_edge_start_point = f2[j].front();
                        }

                        //Jump to next iteration
                        continue;
                }

                //First point the first segment is not a start point
                if ( ( first ) && ( !t2 ) )
                {
                        if ( overlay_type != Difference )
                        {
                                j == f2.size() - 1  ?  j = 0 : j++;
                        }

                        else
                        {
                                j == 0  ?  j = f2.size() - 1 : j--;
                        }
                }

                //First point of the second segmnent is not a start point
                else if ( ( second ) && ( !t1 ) )
                {
                        i == f1.size() - 1 ? i = 0 : i++;
                }

                //Increment min_dist: increase treshold to find a suitable segment
                if ( ++empty_loops > std::max ( f1.size(), f2.size() ) ) min_dist *= 2;
        }
}


#endif

