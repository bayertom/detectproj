// Description: KD-tree implementation (2D)

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


#ifndef KDTree2D_HPP
#define KDTree2D_HPP

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "libalgo/source/const/Const.h"

#include "KDNode.h"
#include "libalgo/source/structures/point/Point3DCartesian.h"
#include "libalgo/source/structures/list/IndexLists.h"

#include "libalgo/source/comparators/sortPointsByX.h"
#include "libalgo/source/comparators/sortPointsByY.h"
#include "libalgo/source/comparators/sortPointsByID.h"

#include "libalgo/source/exceptions/BadOutputException.h"


template <typename Point>
void KDTree2D <Point> ::createKDTree2D ( Container <Point *> &pl, const bool print_exception )
{
        //Create KD tree
        const unsigned int n = pl.size();

        //Create indices lists for sorted sets
        TIndexList kd_index_list_x, kd_index_list_y;

        //Sort points by x
        std::sort ( pl.begin(), pl.end(), sortPointsByX() );
        pl.toIndexList ( kd_index_list_x );

        //Sort points by y
        std::sort ( pl.begin(), pl.end(), sortPointsByY () );
        pl.toIndexList ( kd_index_list_y );

        //Sort points by point id
        std::sort ( pl.begin(), pl.end(), sortPointsByID <Point>() );

        //Build KD Tree
        root = buildKDTree2D ( pl, kd_index_list_x, kd_index_list_y, n, 1 );
}


template <typename Point>
KDNode <Point> * KDTree2D <Point>::buildKDTree2D ( Container <Point *> &pl, TIndexList & kd_index_list_x, TIndexList & kd_index_list_y, const unsigned int n, const unsigned int depth )
{
        //Build 2D KD Tree using recursive approach

        //No leaf will be created
        if ( n == 0 )
        {
                return NULL;
        }

        //Only one point: create leaf of KD Tree
        else if ( n == 1 )
        {
                //Increment nodes count
                nodes_count++;

                //Create one leaft: even depth of recursion
                if ( depth % 2 == 0 )
                {
                        return new KDNode <Point> ( pl [ kd_index_list_x [0]] , depth );
                }

                //Create one leaft: odd depth of recursion
                else
                {
                        return new KDNode <Point> ( pl [ kd_index_list_y [0]] , depth );
                }
        }

        //At least 2 points: create one leaf, split tree into left and right subtree
        else
        {
                //New KD node
                KDNode <Point> *node_median = NULL;

                //Increment nodes count
                nodes_count++;

                //Get median index
                const unsigned int median_index = n / 2;

                //Create new KD Lists
                TIndexList kd_index_list1_x, kd_index_list2_x, kd_index_list1_y, kd_index_list2_y;

                //Initialize median point
                Point * p_median = NULL;

                //The depth is even, process by x coordinate
                unsigned int median_index_point = 0;

                if ( depth % 2 == 0 )
                {
                        //Get median point index
                        median_index_point = kd_index_list_x [median_index];

                        //Get median point
                        p_median = pl [median_index_point];
                }

                //The depth is odd, process by y coordinates
                else
                {
                        //Get median point index
                        median_index_point = kd_index_list_y [median_index];

                        //Get median point
                        p_median = pl [median_index_point] ;
                }

                //Create median node
                node_median = new KDNode <Point> ( p_median, depth );

                //Set split value for the node
                const typename Point::Type split_value = p_median->getCoordinate ( depth % 2 );
                node_median->setSplitValue ( split_value );

                //Split list (n itemd) into two sublists (n1, n2 items)
                unsigned int n1 = 0, n2 = 0;
                splitLists ( median_index_point, split_value, pl, kd_index_list_x, kd_index_list_y, kd_index_list1_x, kd_index_list1_y, kd_index_list2_x, kd_index_list2_y, n, n1, n2, depth );

                //Build left subtree
                node_median->setLeft ( buildKDTree2D ( pl,  kd_index_list1_x, kd_index_list1_y, n1, depth + 1 ) );

                //Build right subtree
                node_median->setRight ( buildKDTree2D ( pl, kd_index_list2_x, kd_index_list2_y, n2, depth + 1 ) );

                //Return new node
                return node_median;
        }
}


template <typename Point>
void KDTree2D <Point>::splitLists ( const unsigned int median_index_point, const typename Point::Type median_split_value, Container <Point*>  &pl,
                                    const TIndexList & kd_index_list_x, const TIndexList & kd_index_list_y, TIndexList & kd_index_list1_x, TIndexList & kd_index_list1_y,
                                    TIndexList & kd_index_list2_x, TIndexList & kd_index_list2_y, const unsigned int n, unsigned int & n1, unsigned int & n2, const unsigned int depth )
{
        //Split KD lists:
        //		depth = even: split according to x coordinate
        //		depth = odd: split according to y coordinate
        //

        //Split list according to median coordinate
        for ( unsigned int i = 0; i < n; i++ )
        {
                //***** Process first list *****

                //Get actual point from the list sorted by x
                const unsigned int point_index_x = kd_index_list_x [i];
                const Point * p_x = pl [point_index_x];

                //Do not process median point again
                if ( point_index_x != median_index_point )
                {
                        //Add point to the first list: x <= median.x (for old depth of recursion)
                        if ( p_x->getCoordinate ( depth % 2 ) <= median_split_value )
                        {
                                kd_index_list1_x.push_back ( point_index_x );
                                n1++;
                        }

                        //Add point to the second list: x < median.x (for old depth of recursion)
                        else if ( p_x->getCoordinate ( depth % 2 ) > median_split_value )
                        {
                                kd_index_list2_x.push_back ( point_index_x );
                                n2++;
                        }
                }

                // ***** Process second list *****

                //Get actual point from the list sorted by y
                const unsigned int point_index_y = kd_index_list_y [i];
                const Point * p_y = pl [point_index_y];

                //Do not process median point again
                if ( point_index_y != median_index_point )
                {
                        //Add point to the first list: y < median.y (for old depth of recursion)
                        if ( p_y->getCoordinate ( depth % 2 ) <= median_split_value )
                        {
                                kd_index_list1_y.push_back ( point_index_y );
                        }

                        //Add point to the second list: y < median.y (for old depth of recursion)
                        else if ( p_y->getCoordinate ( depth % 2 ) > median_split_value )
                        {
                                kd_index_list2_y.push_back ( point_index_y );
                        }
                }
        }

}


template <typename Point>
KDNode <Point> *KDTree2D <Point>:: findParentKDNode ( const Point * point ) const
{
        //Find parent node in KD Tree
        KDNode <Point> *parent = NULL;
        KDNode <Point> *found = findNode ( point, root, &parent, 1 );

        //Point is not present in KD tree
        if ( found == NULL )
        {
                return parent;
        }

        //Point is present in KD Tree
        else
        {
                return found;
        }
}


template <typename Point>
template <typename Point2>
KDNode <Point> *KDTree2D <Point>:: findKDNode ( const Point2 * point, KDNode <Point> **parent ) const
{
        //Find node identical with tested point and its predecessor
        *parent = NULL;

        //No KD-tree constructed
        if ( root == NULL )
        {
                throw BadOutputException ( "BadOutputException: can not find KD-Node, ", "no KD-Tree has been constructed." );
        }

        //Start finding
        return findNode ( point, root, parent, 1 );
}


template <typename Point>
template <typename Point2>
KDNode <Point> *KDTree2D <Point>::findNode ( const Point2 * point, KDNode <Point> *node, KDNode <Point> **parent, const unsigned int depth ) const
{
        //Find node identical with tested point and its predecessor: recursive procedure
        //Point is not present in KD Tree
        if ( node == NULL )
        {
                return NULL;
        }

        //Tested point has same coordinates as found node
        if ( ( point->getX() - node->getData()->getX() ) * ( point->getX() - node->getData()->getX() ) +
                        ( point->getY() - node->getData()->getY() ) * ( point->getY() - node->getData()->getY() ) < MIN_POSITION_DIFF * MIN_POSITION_DIFF )
        {
                return node;
        }

        //Remember parent node
        *parent = node;

        //Search left subtree
        if ( point->getCoordinate ( depth % 2 ) <= node->getData()->getCoordinate ( depth % 2 ) )
        {
                return  findNode ( point, node->getLeft(), parent, depth + 1 );
        }

        //Search right subtree
        else
        {
                return findNode ( point, node->getRight(), parent, depth + 1 );
        }
}


template <typename Point>
template <typename Point2>
Point * KDTree2D <Point>::findNN ( const Point2 * point ) const
{
        //Find nearest Neighbour node and return containing point
        typename Point::Type min_dist_q_node = MAX_FLOAT;
        KDNode <Point> *nearest_neighbour = NULL;
        findNNNode ( point, root, &nearest_neighbour, 1, min_dist_q_node );
        return nearest_neighbour->getData();
}


template <typename Point>
template <typename Point2>
void KDTree2D <Point> ::findNNNode ( const Point2 * point, KDNode <Point> *node, KDNode <Point> **nearest_neighbour, const unsigned int depth, typename Point::Type & min_dist_q_node ) const
{
        //Find node identical with tested point and its predecessor: recursive procedure
        //Point is not present in KD Tree
        if ( node == NULL )
        {
                return;
        }

        //Search left subtree
        if ( point->getCoordinate ( depth % 2 ) <= node->getData()->getCoordinate ( depth % 2 ) )
        {
                findNNNode ( point, node->getLeft(), nearest_neighbour, depth + 1, min_dist_q_node );
        }

        //Search right subtree
        else
        {
                findNNNode ( point, node->getRight(), nearest_neighbour, depth + 1, min_dist_q_node );
        }

        //Compute test distance (actual node and point) and remember actual nearest neighbour and nearest distance
        typename Point::Type dist_q_node = ( node->getData()->getX() - point->getX() ) * ( node->getData()->getX() - point->getX() ) +
                                           ( node->getData()->getY() - point->getY() ) * ( node->getData()->getY() - point->getY() );

        if ( dist_q_node < min_dist_q_node )
        {
                min_dist_q_node = dist_q_node;
                *nearest_neighbour = node;
        }

        //Compute distance from point to intersection of the hypersphere and hyperrectangle
        typename Point::Type dist_q_node_straight = ( point->getCoordinate ( node->getDepth() % 2 ) - node->getData()->getCoordinate ( node->getDepth() % 2 ) ) *
                        ( point->getCoordinate ( node->getDepth() % 2 ) - node->getData()->getCoordinate ( node->getDepth() % 2 ) ) ;

        //There is an intersection of the hypershpere and hyperrectangle
        if ( dist_q_node_straight < min_dist_q_node )
        {
                //Point is in the left halfplane (subtree)
                if ( point->getCoordinate ( node->getDepth() % 2 ) <= node->getData()->getCoordinate ( node->getDepth() % 2 ) )
                {
                        //Continue with right subtree
                        findNNNode ( point, node->getRight(), nearest_neighbour, depth + 1, min_dist_q_node );
                }

                //Point is in the right halfplane (subtree)
                else
                {
                        //Continue with left subtree
                        findNNNode ( point, node->getLeft(), nearest_neighbour, depth + 1, min_dist_q_node );
                }
        }
}


template <typename Point>
template <typename Point2>
void KDTree2D <Point>::findAllKNN ( const Point2 * point, Container <Point *, NonDestructable > &knn, unsigned int k ) const
{
        //Find all k-nearest neighbours to the specified point
        typename TNNeighboursList <Point>::Type nn;

        //No KD-tree constructed
        if ( root == NULL )
        {
                throw BadOutputException ( "BadOutputException: can not find all k-nearest neighbour nodes, ", "no KD-Tree has been constructed." );
        }

        //Set new k
        k = std::min ( k, nodes_count - 1 );

        //Start from root
        findAllKNNNodes ( point, nn, k, root, 1 );

        //Copy nearest neighbours to the list
        typename TNNeighboursList <Point>::Type ::reverse_iterator i_knn = nn.rbegin();

        for ( unsigned int i = 0; i < k; ++i_knn, ++i )
        {
                knn.push_back ( ( ( i_knn )->node->getData() ) );
        }
}


template <typename Point>
template <typename Point2>
void KDTree2D <Point>::findAllKNNNodes ( const Point2 * point, typename TNNeighboursList <Point>::Type & knn, unsigned int k, KDNode <Point> *node, const unsigned int depth ) const
{
        //Find all k nearest neighbours using recursive search and store in Priority Queue (PQ)

        //Point is not present in KD Tree
        if ( node == NULL )
        {
                return;
        }

        //Search left subtree
        if ( point->getCoordinate ( depth % 2 ) <= node->getData()->getCoordinate ( depth % 2 ) )
        {
                findAllKNNNodes ( point, knn, k, node->getLeft(), depth + 1 );
        }

        //Search right subtree
        else
        {
                findAllKNNNodes ( point, knn, k, node->getRight(), depth + 1 );
        }

        //Compute test distance
        typename Point::Type dist_q_node = ( node->getData()->getX() - point->getX() ) * ( node->getData()->getX() - point->getX() ) +
                                           ( node->getData()->getY() - point->getY() ) * ( node->getData()->getY() - point->getY() );

        //PQ does not have any free position: remove item on the top or do not do anything?
        if ( knn.size() == k )
        {
                //We add to PQ only a node closer to q than node on the top
                if ( dist_q_node < knn.begin()->priority )
                {
                        //Remove element from the top: removed node is further to q than added
                        knn.erase ( knn.begin() );

                        //Add node closer to Q to PQ
                        knn.insert ( TKDNodePriority <Point> ( node,  dist_q_node ) );
                }
        }

        //There are free positions in PQ
        else
        {
                //Add node to PQ
                knn.insert ( TKDNodePriority <Point> ( node,  dist_q_node ) );
        }

        //Test intersection of the hypersphere and hyperrectangle
        typename Point::Type dist_q_node_straight = ( point->getCoordinate ( node->getDepth() % 2 ) - node->getData()->getCoordinate ( node->getDepth() % 2 ) ) *
                        ( point->getCoordinate ( node->getDepth() % 2 ) - node->getData()->getCoordinate ( node->getDepth() % 2 ) ) ;

        //We stil have not found enough neighbours or there is an itersection of the largest hypersphere (top of the PQ) and hyperrectangle
        typename Point::Type top_priority =  knn.begin()->priority;

        if ( ( knn.size() < k ) || ( dist_q_node_straight ) < ( top_priority ) )
        {
                //Point is in the left halfplane (subtree)
                if ( point->getCoordinate ( node->getDepth() % 2 ) <= node->getData()->getCoordinate ( node->getDepth() % 2 ) )
                {
                        //Continue with right subtree
                        findAllKNNNodes ( point, knn, k, node->getRight(), depth + 1 );
                }

                //Point is in the right halfplane (subtree)
                else
                {
                        //Continue with left subtree
                        findAllKNNNodes ( point, knn, k, node->getLeft(), depth + 1 );
                }
        }
}


template <typename Point>
void KDTree2D <Point> ::clearKDTree2D()
{
        //Clear KD tree
        clearAllNodes ( root );
}


template <typename Point>
void KDTree2D <Point> ::clearAllNodes ( const KDNode <Point> *node )
{
        //Clear all nodes in KD Tree
        if ( node != NULL )
        {
                //Clear left sub tree
                clearAllNodes ( node->getLeft() );

                //Clear right sub tree
                clearAllNodes ( node->getRight() );

                //Call destructor for KD node
                delete node;
                node = NULL;
        }
}


template <typename Point>
void KDTree2D <Point>:: printKDTree2D ( const TPrintKDTreeMethod & method, std::ostream * output ) const
{
        // Print KD tree:
        //	0: pre-order
        //	1: in-order
        //	2: post-order
        //	oher: post-order
        //

        //No KD-tree constructed
        if ( root == NULL )
        {
                throw BadOutputException ( "BadOutputException: can not print KD-tree, ", "no KD-Tree has been constructed." );
        }

        //Set method
        switch ( method )
        {
                        //Print tree pre-order
                case PreOrder:
                        {
                                *output << "Print KD-tree (pre-order)" << std::endl;
                                preOrder ( root, output );
                        }
                        break;

                        //Print tree in-order
                case InOrder:
                        {
                                *output << "Print KD-tree (in-order)" << std::endl;
                                inOrder ( root, output );
                        }
                        break;

                        //Print tree post-order
                default:
                        {
                                *output << "Print KD-tree (post-order)" << std::endl;
                                postOrder ( root, output );
                        }
        }
}


template <typename Point>
void KDTree2D <Point> :: preOrder ( const KDNode <Point> *node, std::ostream * output ) const
{
        //Print KD tree pre-order
        if ( node == NULL )
        {
                return;
        }

        //Recursively process subtrees
        else
        {
                *output << node->getData();
                preOrder ( node->getLeft() );
                preOrder ( node->getRight() );
        }
}


template <typename Point>
void KDTree2D <Point>:: inOrder ( const KDNode <Point> *node, std::ostream * output ) const
{
        //Print KD tree in-order
        if ( node == NULL )
        {
                return;
        }

        //Recursively process subtrees
        else
        {
                inOrder ( node->getLeft() );
                *output << node->getData();
                inOrder ( node->getRight() );
        }
}


template <typename Point>
void KDTree2D <Point> :: postOrder ( const KDNode <Point> *node, std::ostream * output ) const
{
        //Print KD tree post-order
        if ( node == NULL )
        {
                return;
        }

        //Recursively process subtrees
        else
        {
                *output << node->getData();
                postOrder ( node->getLeft() );
                postOrder ( node->getRight() );
        }
}

#endif
