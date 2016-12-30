// Description: Compute inner distance of the face

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


#ifndef InnerDistance_HPP
#define InnerDistance_HPP

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/face/Face.h"
#include "libalgo/source/structures/graph/GraphM.h"

#include "libalgo/source/algorithms/linelineangle/LineLineAngle.h"
#include "libalgo/source/algorithms/bisector/Bisector.h"
#include "libalgo/source/algorithms/shapecontext/ShapeContext.h"
#include "libalgo/source/algorithms/graphalgorithms/GraphAlgorithms.h"


template <typename T>
T InnerDistance::compare2FacesUsingInnerDistances ( const Face <T> *f1, const Face <T> *f2 )
{
        //Compare two Faces using NN distances
        const unsigned int m1 = f1->getVerticesCount();
        const unsigned int m2 = f2->getVerticesCount();

        //Create matrix represenatation of graphs
        GraphM <T> g1 ( m1, false );
        GraphM <T> g2 ( m2, false );

        //Create empty matrices of  w-distances
        Matrix <T> D1 ( m1, m1, MAX_FLOAT, 0 );
        Matrix <T> D2 ( m2, m2, MAX_FLOAT, 0 );

        //Compute inner distances: convert face to w-distance matrix, normalize edges
        GraphAlgorithms::FaceToInnerDistanceGraph ( f1, g1, true );
        GraphAlgorithms::FaceToInnerDistanceGraph ( f2, g2, true );

        //Matrix of predecessors
        Matrix <unsigned int> P1 ( m1, m1 );
        Matrix <unsigned int> P2 ( m2, m2 );

        /*
        g1.getW().print();
        g2.getW().print();
        */

        //Process matrices using Floyd-Warshall algorithm
        GraphAlgorithms::floydWarshall ( g1, D1, P1 );
        GraphAlgorithms::floydWarshall ( g2, D2, P2 );

        //P1.print();
        //P2.print();

        //Create empty matrices for theta angles
        Matrix <unsigned short> T1 ( m1, m1 );
        Matrix <unsigned short> T2 ( m2, m2 );

        //Compute theta angle for each start segment of the inner distance
        computeTheta ( f1, D1, P1, T1 );
        computeTheta ( f2, D2, P2, T2 );
        /*
        D1.print();
        D2.print();

        T1.print();
        T2.print();
        */

        //List of shape context matrices
        typedef std::vector <Matrix <unsigned short> > TSCList;
        TSCList sc_list;

        //Compute shape context for each point ( given by index i ) of the first Face
        for ( unsigned int i = 0; i < m1; i++ )
        {
                //Create empty shape context matrix
                Matrix <unsigned short> SC1 ( 5, 12 );

                //Compute shape context of the first face
                ShapeContext::computeShapeContext ( D1, T1, SC1, i );

                //Add shape context to the list
                sc_list.push_back ( SC1 );
        }

        //Create shape context difference matrix
        Matrix <T> SCD ( m1, m2 );

        //Compute shape context for each point ( given by index j) of the second Face
        for ( unsigned int j = 0; j < m2; j++ )
        {
                //Create empty shape context matrix
                Matrix <unsigned short> SC2 ( 5, 12 );

                //Compue shape context of the second face
                ShapeContext::computeShapeContext ( D2, T2, SC2, j );

				//SC2.print();

                //Compare only non-zero shape contexts
                if ( MatrixOperations::sum ( SC2 ) > 0 )
                {
                        Matrix <int> DSC ( 5, 12 );

                        //Compute shape contexts difference using X quadrat criterion
                        for ( unsigned int i = 0; i < m1; i++ )
                        {
                                //Compare only non-zero shape contexts
                                if ( MatrixOperations::sum ( sc_list[i] ) > 0 )
                                {
                                        SCD ( i, j ) = ShapeContext::compare2ShapeContexts <T> ( sc_list[i], SC2 );
                                }
                        }
                }
        }

        //Create bipartite matching matrix
        Matrix <unsigned short> MSC ( m1, m2 );

        //Best bipartite matching using Hungarian algorithm
        T matching_cost = 0;
	GraphAlgorithms::bestBipartiteMatching(SCD, MSC, matching_cost);

        //Return cost of the matching
        /*
        if (matching_cost > 0)
        {
        	Container <Face <double> *, NonDestructable> faces1, faces2;
                faces1.push_back(const_cast <Face <double> *> (f1));
        	Container <Point3DCartesian <double > > l1, l2;
        	f1->toPointsList(l1);
        	f2->toPointsList(l2);
        	l1.print();
        	l2.print();
        	D1.print();
        	D2.print();
        	T1.print();
        	T2.print();
        	P1.print();
        	P2.print();
        	//Matrix <int> DTT ( m1, m2 );
        	//DTT = T1 - T2;
        	//DTT.print();
        	//g1.getW().print()
        	SCD.print();
        	faces2.push_back(const_cast <Face <double> *> (f2));
                //DXFExport::exportFacesToDXF ( "D:\\Tomas\\Cpp\\DetectProj\\DetectProj\\out\\f1.dxf", &faces1 );
        	//DXFExport::exportFacesToDXF ( "D:\\Tomas\\Cpp\\DetectProj\\DetectProj\\out\\f2.dxf", &faces2 );

        }
        */

        return matching_cost;
}



template <typename T>
void InnerDistance::computeTheta ( const Face <T> *f, const Matrix <T> &D, const Matrix <int> &P, Matrix <unsigned short> &TH )
{
        //Compute angle of the first edge (pj, pl) of the inner distance to tangent of the Face on point pj
        const unsigned int m = P.cols();

        //Convert actual Face to list of points
        Container <Point3DCartesian <T> > pl;
        f->toPointsList ( pl );

        //Get triplet of points
        Point3DCartesian <T> *pii = &pl[0];

        //Find point pi that is not closer to pii than MIN_POSITION_DIFF: jump closer vertices
        Point3DCartesian <T> *pi = &pl[m - 1];

        for ( unsigned int i = m - 1; EuclDistance::getEuclDistance2D ( pi, pii ) < MIN_POSITION_DIFF; pi = &pl[--i] ) {}

        //Compute theta for each non-duplicate vertex of the face
        for ( unsigned int i = 0; ; )
        {
                //Find point piii not closer to pii than MIN_POSITION_DIFF: jump closer vertices
                Point3DCartesian <T> *piii = NULL;
                unsigned int i_new = ( i == m - 1 ? 0 : i + 1 );

                for ( piii = &pl[i_new]; EuclDistance::getEuclDistance2D ( pii, piii ) < MIN_POSITION_DIFF; ( i_new == m - 1  ? i_new = 0 : i_new++ ), piii = &pl[i_new] ) {}

                //Compute normalized bisector
                T n_bis_x = 0, n_bis_y = 0;
                Bisector::getNormalizedBisector2D ( pi, pii, piii, n_bis_x, n_bis_y );

                //Gompute tangent vector (normal to normal
                const T t_x = -1.0 * n_bis_y;
                const T t_y = n_bis_x;

                //Compute theta angle
                for ( unsigned int k = 0; k < m; k++ )
                {
                        //There is a valid predecessor and edge ( p[i], p[i-1] ) > MIN_POSITION_DIFF
                        if ( ( P ( i, k ) != 0 ) /*&& ( D ( k, ( k == 0 ? m - 1: k - 1 ) ) < 1 )*/ )
                        {
                                //Get first point of the shortest path between points p(i) and p(j) as the last point of the
                                //shortest path from between points p(j) and p(i) (i.e. predecessor of p[j]), item P[j][i] of the matrix
                                const unsigned int pred_index = P ( k, i ) - 1;

                                //First segment of the shortest path has a non-zero length
                                if ( D ( i, pred_index ) < 1 )
                                {
                                        //Set item to theta matrix
                                        Point3DCartesian <T> p_temp ( pii->getX() + t_x , pii->getY() + t_y);
                                        TH ( i, k ) = LineLineAngle::getLineLineAngle ( &p_temp, pii, pii, &pl[pred_index] ) + 0.5 ;
                                }
                        }
                }

                //Testm, if we processed all poits of the and face
                if ( i_new < i ) break;

                //Set new index i
                i = i_new;

                //Assign points
                pi = pii;
                pii = piii;
        }
}

#endif
