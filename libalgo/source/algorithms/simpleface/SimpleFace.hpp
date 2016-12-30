// Description: Check whether a face is simple or not

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


#ifndef SimpleFace_HPP
#define SimpleFace_HPP


#include "libalgo/source/structures/point/Node3DCartesian.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/eucldistance/EuclDistance.h"
#include "libalgo/source/algorithms/facearea/FaceArea.h"
#include "libalgo/source/algorithms/linelineposition/LineLinePosition.h"


template <typename T>
bool SimpleFace::isSimpleFace ( const Face <T> *face )
{
        //Test, if face is simple
        if ( face != NULL )
        {
                const HalfEdge <T> *e_start = face->getHalfEdge();

                //There is a valid edge
                if ( e_start != NULL )
                {
                        //Face is a triangle
                        if ( face->getVerticesCount() == 3 )
                        {
                                //Compute area of the triangle
                                if ( FaceArea::getFaceArea ( face ) > AREA_ROUND_ERROR )
                                {
                                        return true;
                                }

                                return false;
                        }

                        //Get start point of the first edge
                        const Node3DCartesian <T> *p1 = e_start->getPoint();

                        //Proces all edges of the Face
                        for ( HalfEdge <T> *e1 = const_cast <HalfEdge <T> *> ( e_start ); e1 != e_start ->getPreviousEdge()->getPreviousEdge() ; e1 = e1->getNextEdge () )
                        {
                                //Get end point of the first edge
                                const Node3DCartesian <T> *p2 = e1->getNextEdge()->getPoint();

                                //First segment does not have to be too short
                                if ( EuclDistance::getEuclDistance2D ( p1->getX(), p1->getY(), p2->getX(), p2->getY() ) > MIN_POSITION_DIFF )
                                {
                                        //Proces all edges of the Face
                                        for ( HalfEdge <T> *e2 = e1->getNextEdge()->getNextEdge(); ( e2 != e_start ) && ( e2 != e1->getPreviousEdge() ); e2 = e2->getNextEdge () )
                                        {
                                                //Get start point of the second edge
                                                const Node3DCartesian <T> *p3 = e2->getPoint();

                                                //Start point of the second segment is not to close to end point of the first segment
                                                if ( EuclDistance::getEuclDistance2D ( p2->getX(), p2->getY(), p3->getX(), p3->getY() ) > MIN_POSITION_DIFF )
                                                {
                                                        //Get end point of the second edge
                                                        const Node3DCartesian <T> *p4 = e2->getNextEdge()->getPoint();

                                                        //Second segment does not have to be too short, end point of the second segment is not to close to first point of the first segment
                                                        if ( EuclDistance::getEuclDistance2D ( p3->getX(), p3->getY(), p4->getX(), p4->getY() ) > MIN_POSITION_DIFF &&
                                                                        EuclDistance::getEuclDistance2D ( p4->getX(), p4->getY(), p1->getX(), p1->getY() ) > MIN_POSITION_DIFF )
                                                        {
                                                                //Test intersection: do not round results
                                                                T x_int, y_int;
                                                                unsigned int t = LineLinePosition::get2LineSegmentsPosition ( p1, p2, p3, p4, x_int, y_int );

                                                                if ( ( t > 0 ) && ( t < 4 ) || ( t == 5 ) )
                                                                {
                                                                        //A face is not simple
                                                                        /*
                                                                        std::cout << "\n Face not simple: \n";
                                                                        Container <Face <T> *, NonDestructable> faces;
                                                                        faces.push_back(const_cast <Face <T> *> (face));
                                                                        Container <Node3DCartesian <T> *, NonDestructable > f_nodes;
                                                                        face->toNodesList ( &f_nodes );
                                                                        f_nodes.printItems();
                                                                        DXFExport::exportFacesToDXF ( "D:\\Tomas\\Cpp\\GridDetection\\GridDetection\\out\\simple_face_error.dxf", &faces );

                                                                        std::cout<< "First segment: \n";
                                                                        p1->print();
                                                                        p2->print();

                                                                        std::cout << "Second segment: \n";
                                                                        p3->print();
                                                                        p4->print();
                                                                        */
                                                                        return false;
                                                                }
                                                        }
                                                }
                                        }
                                }

                                //Assign points
                                p1 = p2;
                        }
                }

                //Face is simple
                return true;
        }

        //Throw exception
        else
        {
                throw BadDataException ( "BadDataException: can not test, if a face is simple,", " a line segment = NULL." );
        }
}

#endif
