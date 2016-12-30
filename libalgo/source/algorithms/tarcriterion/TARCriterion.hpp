// Description: Compute TAR criterion

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


#ifndef TARCriterion_HPP
#define TARCriterion_HPP

#include <algorithm>
#include <vector>
#include <math.h>


#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"
#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/trianglearea/TriangleArea.h"
#include "libalgo/source/algorithms/faceperimeter/FacePerimeter.h"
#include "libalgo/source/algorithms/graphalgorithms/GraphAlgorithms.h"


template <typename T>
T TARCriterion::compare2FacesUsingTARCriterion ( const Face <T> *p1, const Face <T> *p2 )
{
        //Compare two Faces using normalized TAR criterion: cost of the bipartite matching of the TAR'S distances
        T matching_cost = 0;
        T sc1 = 0, sc2 = 0;

        //Get HalfEdge
        HalfEdge <T> *e1 = p1->getHalfEdge();
        HalfEdge <T> *e2 = p2->getHalfEdge();

        //There is a valid edge
        if ( ( e1 != NULL ) && ( e2 != NULL ) )
        {
                //Get vertices count
                unsigned int n1 = p1->getVerticesCount(), n2 = p2->getVerticesCount();

                //Compute max vertices step
                const unsigned int max_vertices_step = std::min ( n1 + 1 , n2 + 1 ) / 3;

                //Matrix of TAR's for the first polygon
                Matrix <T> TAR1 ( n1, max_vertices_step ) ;

                //Create matrix of tar differences
                Matrix <T> D ( n1, n2 );

                //Assign first edge of the first face
                HalfEdge <T> *e_start = e1;

                //Proces all edges of the first face: compute tar criterion for different level given by vertices_step value
                unsigned int i = 0, n1_unique = 0;

                do
                {
                        //Jump duplicate point
                        if ( EuclDistance::getEuclDistance2D ( e1->getPoint(), e1->getPreviousEdge()->getPoint() ) > MIN_POSITION_DIFF )
                        {
                                //Compute TAR1 for different levels
                                T min1 = 1, max1 = 0, tar1 = 0;

                                for ( unsigned int vertices_step = 1; vertices_step <= max_vertices_step; vertices_step ++ )
                                {
                                        TAR1 ( i, vertices_step - 1 ) = ( tar1 = 1000 * getTARCriterion ( e1, vertices_step ) );

                                        if ( tar1 < min1 ) min1 = tar1;

                                        if ( tar1 > max1 ) max1 = tar1;
                                }

                                //Compute sc1
                                sc1 += fabs ( max1 - min1 );

                                n1_unique++;
                        }

                        //Increment edge
                        e1 = e1->getNextEdge();

                        //Increment index
                        i++;

                }
                while ( e1 != e_start );

                //Assign first edge of the second face
                e_start = e2;

                //Proces all edges of the first face: compute tar criterion for different level given by vertices_step value
                unsigned int j = 0, n2_unique = 0;

                do
                {
                        //Jump duplicate point
                        if ( EuclDistance::getEuclDistance2D ( e2->getPoint(), e2->getPreviousEdge()->getPoint() ) > MIN_POSITION_DIFF )
                        {
                                //Compute tar differences
                                T min2 = 1, max2 = 0, tar2 = 0;
                                Matrix <T> TAR2 ( 1, max_vertices_step );

                                for ( unsigned int vertices_step = 1; vertices_step <= max_vertices_step; vertices_step ++ )
                                {
                                        TAR2 ( 0, vertices_step - 1 ) = ( tar2 = 1000 * getTARCriterion ( e2, vertices_step ) );

                                        if ( tar2 < min2 ) min2 = tar2;

                                        if ( tar2 > max2 ) max2 = tar2;
                                }

                                //Compute sc2
                                sc2 += fabs ( max2 - min2 );

                                //Compute TAR difference matrix D
                                for ( unsigned int k = 0; k < n1; k++ )
                                {
                                        for ( unsigned int vertices_step = 1; vertices_step <= max_vertices_step; vertices_step ++ )
                                        {
                                                D ( k, j ) += fabs ( TAR1 ( k, vertices_step - 1 )  - TAR2 ( 0, vertices_step - 1 ) );
                                        }
                                }

                                n2_unique++;
                        }

                        //Increment edge
                        e2 = e2->getNextEdge();

                        //Increment index
                        j++;

                }
                while ( e2 != e_start );

                //Compute DT using min of its 4 neighbours: applicable only to raster data
                Matrix <T> DT ( D );
                /*
                for ( unsigned int i = 1; i < n1; i++ )
                {
                        for ( unsigned int j = 1; j < n2; j ++ )
                        {
                                DT ( i, j ) = D ( i, j ) +  min ( min ( D ( i - 1 , j ), D ( i - 1, j - 1 ) ), D ( i, j - 1 ) );
                        }
                }
                */
                //Create bipartite matching matrix DTMIN
                Matrix <unsigned short> DTMIN ( n1, n2 );

                //Best bipartite matching using Hungarian algorithm
                GraphAlgorithms::bestBipartiteMatching ( DT, DTMIN, matching_cost );

                //Compute sc1, sc2
                sc1 /= n1_unique; sc2 /= n2_unique;
        }

        //Get matching cost of both Faces
        return matching_cost / ( 1 + sc1 + sc2 );
}


template <typename T>
T TARCriterion::getTARCriterion ( HalfEdge <T> *e, unsigned int vertices_step )
{
        //Get TAR criterion of the point with the specified vertex step
        T tar_normalized = 0;

        if ( e != NULL )
        {
                //Initialize previous and next edge
                HalfEdge <T> *e_prev = e, *e_next = e ;

                //Get Face
                const Face <T> *f = const_cast < Face <T> * > ( e->getFace() );

                if ( f != NULL )
                {
                        //Get i-th previous edge of  the face in one for cycle
                        unsigned int i = 0;

                        for ( ; i < vertices_step; )
                        {
                                e_prev = e_prev->getPreviousEdge();

                                if ( EuclDistance::getEuclDistance2D ( e_prev->getNextEdge()->getPoint(), e_prev->getPoint() ) > MIN_POSITION_DIFF ) i++;
                        }

                        //Get i-th previous edge of  the face in one for cycle
                        unsigned int j = 0;

                        for ( ; j < vertices_step; )
                        {
                                e_next = e_next->getNextEdge();

                                if ( EuclDistance::getEuclDistance2D ( e_next->getPreviousEdge()->getPoint(), e_next->getPoint() ) > MIN_POSITION_DIFF ) j++;
                        }

                        //Get perimeter of the cace
                        const T face_perimeter  = FacePerimeter::getFacePerimeter ( f );

                        //Compute normalized TAR
                        tar_normalized =   TriangleArea::getTriangleAreaSigned ( e_prev->getPoint(), e->getPoint(), e_next->getPoint() )
                                           / ( face_perimeter * face_perimeter ) ;
                }
        }

        //Return normalized TAR
        return  tar_normalized;
}

#endif
