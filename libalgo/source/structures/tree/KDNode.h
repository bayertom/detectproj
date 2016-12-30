// Description: Class storing KD node (2D)

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


#ifndef KDNode_H
#define KDNode_H

#include <stdio.h>

//Class storing KDNode (2D)
template <typename Point>
class KDNode
{
        private:
                Point * data;				//Stored any type of point
                KDNode <Point> *left;			//Pointer to the left node
                KDNode <Point> *right;			//Pointer to the right node
                unsigned int depth;                     //Depth of the recursion
                typename Point::Type split_value;       //Splitting value

        public:
                KDNode() : data ( NULL ), left ( NULL ), right ( NULL ), depth ( 0 ), split_value ( 0 ) {}
                KDNode ( Point * data_ ) : data ( data_ ), left ( NULL ), right ( NULL ), depth ( 0 ), split_value ( 0 )  {}
                KDNode ( Point * data_, const unsigned int depth_ ) : data ( data_ ), left ( NULL ), right ( NULL ), depth ( depth_ ), split_value ( 0 )  {}
                KDNode ( Point * data_, KDNode <Point> *left_, KDNode <Point> *right_ ) : data ( data_ ), left ( left_ ), right ( right_ ), depth ( 0 ), split_value ( 0 )  {}
                KDNode ( Point * data_, KDNode <Point> *left_, KDNode <Point> *right_, const unsigned int depth_, const typename Point::Type split_value_ )
                        : data ( data_ ), left ( left_ ), right ( right_ ), depth ( depth_ ), split_value ( split_value_ ) {}
                ~KDNode() {data = NULL; left = NULL; right = NULL;}

        public:
                Point * getData() const {return data;}
                KDNode <Point> *getLeft() const {return left;}
                KDNode <Point> *getRight() const {return right;}
                unsigned int getDepth() const {return depth;}
                typename Point::Type getSplitValue() const {return split_value;}

        public:
                void setLeft ( KDNode <Point> *left_ ) {left = left_;}
                void setRight ( KDNode <Point> *right_ ) {right = right_;}
                void setDepth ( const unsigned int depth_ ) {depth = depth_;}
                void setSplitValue ( const typename Point::Type split_value_ ) {split_value = split_value_;}
};

#endif

