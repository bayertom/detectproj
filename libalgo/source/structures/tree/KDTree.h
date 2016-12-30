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


#ifndef KDTree2D_H
#define KDTree2D_H

#include <set>
#include <ostream>
#include <iostream>

#include "libalgo/source/structures/list/Container.h"

//Forward declarations
template <typename T>
class KDNode;

class TIndexList;

//Structure for the BFS/DFS processing of the KD tree, storing KD node and its priority
template <typename Point>
struct TKDNodePriority
{
        KDNode <Point> *node;
        typename Point::Type priority;

        TKDNodePriority() : node ( NULL ), priority ( 0 ) {}
        TKDNodePriority ( KDNode <Point> *node_, typename Point::Type priority_ ) : node ( node_ ), priority ( priority_ ) {}

        //Descending sort
        bool operator < ( const TKDNodePriority <Point> &n1 ) const
        {
                return priority > n1.priority;
        }
};


//New user type, multiset of sorted nearest points ( some points could have same distances from q )
template <typename Point>
struct TNNeighboursList
{
        typedef std::multiset < TKDNodePriority <Point> > Type;
};


//Print KD-tree: preorder, inorder, postorder
typedef enum
{
        PreOrder = 0,
        InOrder,
        PostOrder,
} TPrintKDTreeMethod;


//Class storing KD tree 2D: Partial specialization for Point*
template <typename Point>
class KDTree2D
{
        private:
                KDNode <Point> *root;				//Root of the KD tree
                unsigned int nodes_count;       		//Nodes count of the KD tree

        public:
                KDTree2D () : root ( NULL ), nodes_count ( 0 ) {}
                ~KDTree2D() { clearKDTree2D(); }

        public:
                void createKDTree2D ( Container <Point *> &pl, const bool print_exception = true );
                KDNode <Point> *findParentKDNode ( const Point * point ) const;
                void clearKDTree2D();
                void printKDTree2D ( const TPrintKDTreeMethod & method = PostOrder, std::ostream * output = &std::cout ) const;
                unsigned int getNodesCount() const { return nodes_count;}

        public:
                template <typename Point2>
                KDNode <Point> *findKDNode ( const Point2 * point, KDNode <Point> **nearest_neighbour ) const;

                template <typename Point2>
                Point * findNN ( const Point2 * point ) const;

                template <typename Point2>
                void findAllKNN ( const Point2 * point, Container <Point *, NonDestructable> &knn, unsigned int k ) const;


        private:

                KDNode <Point> *buildKDTree2D ( Container <Point *> &pl, TIndexList & kd_index_list_x, TIndexList & kd_index_list_y, const unsigned int n, const unsigned int depth );

                void splitLists ( const unsigned int median_index_point, const typename Point::Type median_split_value, Container <Point *> &pl,
                                  const TIndexList & kd_index_list_x, const TIndexList & kd_index_list_y, TIndexList & kd_index_list1_x, TIndexList & kd_index_list1_y,
                                  TIndexList & kd_index_list2_x, TIndexList & kd_index_list2_y, const unsigned int n, unsigned int & n1, unsigned int & n2, const unsigned int depth );

                void clearAllNodes ( const KDNode <Point> *node );

                void inOrder ( const KDNode <Point> *node, std::ostream * output ) const;
                void preOrder ( const KDNode <Point> *node, std::ostream * output ) const;
                void postOrder ( const KDNode <Point> *node, std::ostream * output ) const;

        private:

                template <typename Point2>
                KDNode <Point> *findNode ( const Point2 * point, KDNode <Point> *node, KDNode <Point> **parent, const unsigned int depth ) const;

                template <typename Point2>
                void findNNNode ( const Point2 * point, KDNode <Point> *node, KDNode <Point> **parent, const unsigned int depth, typename Point::Type & min_dist_q_node ) const; // nove

                template <typename Point2>
                void findAllKNNNodes ( const Point2 * point, typename TNNeighboursList <Point>::Type & knn, unsigned int k, KDNode <Point> *node, const unsigned int depth ) const;
};

#include "KDTree.hpp"

#endif
