// Description: Various NN-distances

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


#ifndef NNDistance_HPP
#define NNDistance_HPP

#include <algorithm>

#include "libalgo/source/structures/tree/KDTree.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"

#include "libalgo/source/exceptions/BadDataException.h"


template <typename Point1, typename Point2>
typename Point1::Type NNDistance::getCrossNearestNeighbourDistance ( const Container <Point1 *> &list1, const Container <Point2 *> &list2 )
{
        //Compute cross nearest neighbour distance for two lists of points
        typename Point1::Type dist = 0;
        const unsigned int n1 = list1.size();
        const unsigned int n2 = list2.size();

        //Throw exception
        if ( n1 != n2 || n1 == 0 || n2 == 0 )
        {
                throw BadDataException ( "BadDataException: can not compute cross nearest diastance,", " both lists have different size." );
        }

        //Create KD trees for neighbours searching
        KDTree2D <Point1> kd_tree1;
        KDTree2D <Point2> kd_tree2;
        kd_tree1.createKDTree2D ( const_cast < Container <Point1 *> &> ( list1 ) );
        kd_tree2.createKDTree2D ( const_cast < Container <Point2 *> &> ( list2 ) );

        //Process all points
        for ( unsigned int i = 0; i < list1.size(); i++ )
        {
                Point1 * point1 = list1 [i];
                Point2 * point2 = list2 [i];

		/*
				//Find nearest points
				typename Point1::Type nearest_dist12 = MAX_FLOAT, nearest_dist21 = MAX_FLOAT;
				Point1 * nearest12 = list2[0], *nearest21 = list1[0];

				for (unsigned int j = 0; j < list1.size(); j++)
				{
					const typename Point1::Type dist12 = EuclDistance::getEuclDistance2D(list1[i]->getX(), list1[i]->getY(), list2[j]->getX(), list2[j]->getY());
					const typename Point1::Type dist21 = EuclDistance::getEuclDistance2D(list1[j]->getX(), list1[j]->getY(), list2[i]->getX(), list2[i]->getY());

					//Find shorter distance d12
					if ( dist12 < nearest_dist12)
					{
						nearest_dist12 = dist12;
						nearest12 = list2[j];
					}

					//Find shorter distance d21
					if (dist21 < nearest_dist21)
					{
						nearest_dist21 = dist21;
						nearest21 = list1[j];
					}
				}
		*/

                //Find nearest point to point from first dataset in KD tree of the second dataset
                Point1 * nearest21 = kd_tree1.findNN ( point2 );

                //Find nearest point to point from second dataset in KD tree of the first dataset
                Point2 * nearest12 = kd_tree2.findNN ( point1 );

                //Compute distances
                dist += EuclDistance::getEuclDistance2D ( nearest21->getX(), nearest21->getY(), list2 [i]->getX(), list2 [i]->getY() ) +
                        EuclDistance::getEuclDistance2D ( nearest12->getX(), nearest12->getY(), list1 [i]->getX(), list1 [i]->getY() );			
        }

        return dist / ( 2 * n1 );
}



template <typename Point1, typename Point2>
typename Point1::Type NNDistance::getAverageNearestNeighbourDistance ( const KDTree2D <Point1> &tree, const Point2 * point,  unsigned int k )
{
        //Compute cross nearest neighbour distance for two lists
        typename Point1::Type average_nn_dist = 0;

        //Throw exceptions: no KD-tree
        if ( tree == NULL )
        {
                throw BadDataException ( "BadDataException: can not compute average nearest neighbour distance,", " KD-tree is empty." );
        }

        //All k-nearest neighbours container, non-destructable
        Container <Point1 *, NonDestructable > knn;

        //Find all k-narest neighbours
        tree.findAllKNN ( point, &knn, k );

        //Compute average distance
        for ( unsigned int i = 0; i < k; i++ )
        {
                average_nn_dist += EuclDistance::getEuclDistance2D ( point, knn[i] ) ;
        }

        return average_nn_dist / k;
}


template <typename Point1, typename Point2>
typename Point1::Type NNDistance::compare2DatasetsUsingAverageNearestNeighbourDistance ( const Container <Point1 *> &list1, const Container <Point2 *> &list2, unsigned int k )
{
        //Compare 2 datasets using ANND ratio
        KDTree2D <Point1> kd_tree1;
        KDTree2D <Point2> kd_tree2;

        //Get total items of the dataset
        const unsigned int n1 = list1.size();
        const unsigned int n2 = list2.size();

        //Throw exception
        if ( n1 != n2 )
        {
                throw BadDataException ( "BadDataException: can not compare 2 datasets using ANND ratio, ", "datasets have different sizes." );
        }

        //Correct k: if there are less points than k
        k = std::min ( k, n1 - 2 );

        //Create KD-trees
        kd_tree1.createKDTree2D ( const_cast < Container <Point1 *> *> ( list1 ) );
        kd_tree2.createKDTree2D ( const_cast < Container <Point2 *> *> ( list2 ) );

        //Compute ANND ratio for both datasets
        typename Point1::Type annd_ratio = 0;

        for ( unsigned int i = 0; i < n1; i++ )
        {
                //Compute annd for each data set
                typename Point1::Type annd1 = getAverageNearestNeighbourDistance ( kd_tree1, ( *list1 ) [i], k );
                typename Point2::Type annd2 = getAverageNearestNeighbourDistance ( kd_tree2, ( *list2 ) [i], k );

                //Compute ANND ratio
                annd_ratio += fabs ( annd1 - annd2 );
        }

        return annd_ratio / n1;
}

#endif
