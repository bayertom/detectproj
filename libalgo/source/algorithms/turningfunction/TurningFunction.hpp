// Description: Compute turning function

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


#ifndef TurningFunction_HPP
#define TurningFunction_HPP


#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/angle3points/Angle3Points.h"
#include "libalgo/source/algorithms/bearing/Bearing.h"
#include "libalgo/source/algorithms/faceperimeter/FacePerimeter.h"
#include "libalgo/source/algorithms/polylinelength/PolyLineLength.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"

template <typename T>
T TurningFunction::compare2FacesUsingTurningFunction ( const Face <T> *f1, const Face <T> *f2, const TTurningFunctionRotationMethod & rotation_method,
                const TTurningFunctionScaleMethod & scale_method )
{
        //Compare 2 faces using difference of turning functions
        //Result does not dependend on the order of vertices:  perform cyclic map rotation
        T rotation1, rotation2;
        typename TTurningFunction <T>::Type tf1, tf2;

        computeTurningFunctionFace ( f1, tf1, rotation1, rotation_method, scale_method );
        computeTurningFunctionFace ( f2, tf2, rotation2, rotation_method, scale_method );

        //Find minimum weight difference of turning functions: unrotated first map and rotated second map
        T difference = MAX_FLOAT, min_difference = MAX_FLOAT;

        //Perform cyclic rotation of the second map
        for ( unsigned int i = 0; i < tf2.size() ; i++ ) // perform n-1 cyclic rotations
        {
                //Zero difference: stop computing
                if ( ( difference = compareTurningFunctions <T> ( tf1, tf2, rotation_method ) ) < ARGUMENT_ROUND_ERROR )
                        return difference;

                //Remember min value
                if ( difference < min_difference ) min_difference = difference;

                //Remember first value in the second map
                const T second_key = ( ++ tf2.begin() ) ->first;
                const T first_value = tf2.begin()->second;

                //Rotate all keys in the second map
                typename TTurningFunction <T>::Type ::iterator i_tf = tf2.begin();

                for ( typename TTurningFunction <T>::Type ::iterator i_tf_next = ++tf2.begin(); i_tf_next != tf2.end(); ++i_tf_next )
                {
                        //Get actual key
                        typename TTurningFunction <T>::Type ::key_type * actual_key = const_cast <typename TTurningFunction <T>::Type ::key_type *> ( &i_tf->first ) ;

                        //Change key and value: perform cycling rotation
                        *actual_key = i_tf_next->first - second_key; i_tf->second = i_tf_next->second;

                        //Subtract first turning angle from all turning angles of the rotated map
                        if ( rotation_method == RotationInvariant ) i_tf->second -= first_value;

                        //Assign old value
                        i_tf = i_tf_next;
                }

                //Correct element at the end of the map:
                // angle[n-1] = angle[0] or angle[n-1] = angle[0] + angle[n-2]
                i_tf->second = ( rotation_method == RotationDependent ? tf2.begin()->second : ( ++tf2.rbegin() )->second + tf2.begin()->second );
        }

        return min_difference;
}



template <typename Point>
typename Point::Type TurningFunction::compare2PolyLinesUsingTurningFunction ( const Container <Point > &pl1, const Container <Point> &pl2, const TTurningFunctionRotationMethod & rotation_method,
                const TTurningFunctionScaleMethod & scale_method )
{
        //Compare 2 polylines using difference of turning functions
        typename TTurningFunction <typename Point::Type>::Type tf1, tf2;

        //Not enough points, throw exception
        if ( pl1.size() < 3 || pl2.size() < 3 )
        {
                throw BadDataException ( "BadDataException: not enough points. ", "Can not compare 2 polylines using Turning function." );
        }

        //Compute turning functions
        computeTurningFunctionPolyLine ( pl1 , tf1, rotation_method, scale_method );
        computeTurningFunctionPolyLine ( pl2 , tf2, rotation_method, scale_method );

        //Compute difference of turning functions
        return compareTurningFunctions <typename Point::Type> ( tf1, tf2, rotation_method );
}


template <typename T>
void TurningFunction::computeTurningFunctionFace ( const Face <T> *face, typename TTurningFunction <T>::Type & tf, T & rotation, const TTurningFunctionRotationMethod & rotation_method,
                const TTurningFunctionScaleMethod & scale_method )
{
        //Compute normalized turning function for the Face  given by HalfEdge
        T turning_angle_sum = 0, normalize_length_sum = 0;
        const HalfEdge <T> *e_start = face->getHalfEdge();
        HalfEdge <T> *e = const_cast <HalfEdge <T> *> ( e_start ) ;

        //Get triplet of points
        Node3DCartesian <T> *pi, *piii = NULL;
        Node3DCartesian <T> *pii = e_start->getPoint();

        //Jump short segments
        for ( pi = e->getPreviousEdge()->getPoint(); EuclDistance::getEuclDistance2D ( pi->getX(), pi->getY(), pii->getX(), pii->getY() ) < MIN_POSITION_DIFF ; e = e->getPreviousEdge(), pi = e->getPreviousEdge()->getPoint() );

        //Get third point for initial edge
        for ( piii = e->getNextEdge()->getPoint(); EuclDistance::getEuclDistance2D ( pii->getX(), pii->getY(), piii->getX(), piii->getY() ) < MIN_POSITION_DIFF ; e = e->getNextEdge(), piii = e->getNextEdge()->getPoint() );

        //Compute perimeter of the Face
        const T face_perimeter = FacePerimeter::getFacePerimeter ( e_start );

        //Throw exception
        if ( fabs ( face_perimeter ) < MIN_FLOAT )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: turning function, can not normalize function (1 / perimeter), ", "perimeter = 0.", face_perimeter );
        }

        //Compute initial angle: rotation invariant
        if ( rotation_method == RotationInvariant )
        {
                turning_angle_sum += 180 - Angle3Points::getAngle3Points ( pi, pii, piii );
        }

        //Compute initial angle: rotation dependent (do not compute the sum)
        else
        {
                //Create temporary point P (x_temp, y_temp) on x axis
                const T dx = std::max ( fabs ( piii->getX() - pii->getX() ), fabs ( piii->getY() - pii->getY() ) );
                Node3DCartesian <T> n_temp ( pii->getX() - 2 * dx , pii->getY() );

                turning_angle_sum += Bearing::getBearing ( pii, &n_temp );
        }

        //Add turning angle computed for normalized length of initial edge to map
        tf[normalize_length_sum] = turning_angle_sum;

        //Assign points from initial edge
        pi = pii;
        pii = piii;

        //Increment initial edge
        e = e_start->getNextEdge();

        //Process all edges of the Face (2, n-1)
        do
        {
                //Get third point
                piii = e->getNextEdge()->getPoint();

                //Segment is too short
                if ( EuclDistance::getEuclDistance2D ( pii->getX(), pii->getY(), piii->getX(), piii->getY() ) > MIN_POSITION_DIFF )
                {
                        //Compute angle: rotation invariant
                        turning_angle_sum += 180.0 - Angle3Points::getAngle3Points ( pi, pii, piii );

                        //Compute normalized length sum
                        normalize_length_sum += ( scale_method == ScaleInvariant ? EuclDistance::getEuclDistance2D ( pi->getX(), pi->getY(), pii->getX(), pii->getY() ) / face_perimeter :
                                                  EuclDistance::getEuclDistance2D ( pi->getX(), pi->getY(), pii->getX(), pii->getY() ) );

                        //Add turning angle computed for normalized length to map
                        tf[normalize_length_sum] = turning_angle_sum;

                        //Assign points
                        pi = pii;
                        pii = piii;
                }

                //Increment edge
                e = e->getNextEdge();

        }
        while ( e != e_start );

        //Add last value: we visit first point again
        const T last_key = ( scale_method == ScaleInvariant ? 1.0 : face_perimeter );
        tf[last_key] = ( rotation_method == RotationInvariant ? turning_angle_sum + tf.begin()->second : tf.begin()->second );
}



template <typename Point>
void TurningFunction::computeTurningFunctionPolyLine ( const Container <Point> &points, typename TTurningFunction <typename Point::Type>::Type & tf, const TTurningFunctionRotationMethod & rotation_method,
                const TTurningFunctionScaleMethod & scale_method )
{
        //Compute normalized turning function for the polyline
        const unsigned int n = points.size();
        typename Point::Type normalize_length_sum = 0, turning_angle_sum = 0;

        //Get polylline length
        const typename Point::Type polyline_length = PolyLineLength::getPolyLineLength ( points );

        //Compute initial angle: rotation invariant
        if ( rotation_method == RotationInvariant )
        {
                turning_angle_sum += 180 - Angle3Points::getAngle3Points ( & points [0], & points [1],  & points [2] );
        }

        //Compute initial angle: rotation dependent
        else
        {
                //Create temporary point on x axis
                const typename Point::Type dx = std:: max ( fabs ( points [2].getX() - points [1].getX() ), fabs ( points [2].getY() - points [1].getY() ) );
                Point p_temp ( points [1].getX() - 2 * dx, points [1].getY() );

                turning_angle_sum = 360 - Angle3Points::getAngle3Points ( & p_temp,  & points [1], & points [2] );
        }

        //Compute normalized length
        normalize_length_sum += EuclDistance::getEuclDistance2D ( points [1].getX(), points [1].getY(),
                                points [2].getX(), points [2].getY() ) / polyline_length;

        //Add turning angle computed for normalized length to map
        tf[normalize_length_sum] = turning_angle_sum;

        //Process all inernal vertices of the polyline
        for ( unsigned int i = 2; i < n - 1; i++ )
        {
                //Compute angle: rotation invariant
                if ( rotation_method == RotationInvariant )
                {
                        turning_angle_sum += 180 - Angle3Points::getAngle3Points ( & points [i - 1], & points [i],  & points [i + 1] );
                }

                else
                {
                        //Create temporary point on x axis
                        const typename Point::Type dx = std::max ( fabs ( points [i + 1].getX() - points [i].getX() ), fabs ( points [i + 1].getY() - points [i].getY() ) );
                        Point p_temp ( points [i].getX() - 2 * dx, points [i].getY() );

                        turning_angle_sum = 360 - Angle3Points::getAngle3Points ( &p_temp,  & points [i], & points [i + 1] );
                }

                //Compute normalized length
                normalize_length_sum += ( scale_method == ScaleInvariant ? EuclDistance::getEuclDistance2D ( points [i].getX(), points [i].getY(), points [i + 1].getX(), points [i + 1].getY() ) / polyline_length :
                                          EuclDistance::getEuclDistance2D ( points [i].getX(), points [i].getY(), points [i + 1].getX(), points [i + 1] .getY() ) );

                //Add turning angle computed for normalized length to map
                tf[normalize_length_sum] = turning_angle_sum;
        }

        //Add last item
        scale_method == ScaleInvariant ? tf[1.0] = turning_angle_sum : tf[polyline_length] = turning_angle_sum;
}



template <typename T>
T TurningFunction::compareTurningFunctions ( const typename TTurningFunction <T>::Type & tf1, const typename TTurningFunction <T>::Type & tf2, const TTurningFunctionRotationMethod & rotation_method )
{
        //Compare two normalized turning functions for faces f1, f2 given by turning functions and returns area difference ratio
        //Method analogous to the 2-way merge sort
        T key_old = 0, key = 0, area_diff = 0;
        const unsigned int n = tf1.size(), m = tf2.size();

        //Iterators
        typename TTurningFunction <T>::Type ::const_iterator i_tf1 = tf1.begin(), i_tf1_i1;
        typename TTurningFunction <T>::Type ::const_iterator i_tf2 = tf2.begin(), i_tf2_i1;

        //Initialize
        T val1 = i_tf1->second, val2 = i_tf2->second;
        i_tf1++; i_tf2++;

        //Merge both maps and compute differencies
        for ( unsigned int i = 2; i < m + n; i++ )
        {
                //Turning angles
                val1 = ( -- ( i_tf1_i1 = i_tf1 ) )->second;
                val2 = ( -- ( i_tf2_i1 = i_tf2 ) )->second;

                //Remember old key
                key_old = key;

                //Last element in the first map
                if ( i_tf1 == tf1.end() )
                {
                        //Get key
                        key = i_tf2->first;

                        //Increment iterator
                        i_tf2 ++;

                        //Jump to the next iteration
                        continue;
                }

                //Last element in the second map
                if ( i_tf2 == tf2.end() )
                {
                        //Get key
                        key = i_tf1->first;

                        //Increment iterator
                        i_tf1 ++;

                        //Jump to the next iteration
                        continue;
                }

                //Map 1 value < Map 2 value
                if ( i_tf1->first < i_tf2 -> first )
                {
                        //Get key
                        key = i_tf1->first;

                        //Increment iterator
                        i_tf1 ++;
                }

                //Map 2 value < Map 1 value
                else
                {
                        //Get key
                        key = i_tf2->first;

                        //Increment iterator
                        i_tf2 ++;
                }

                T val_diff = std::min ( std::min ( fabs ( val1 - val2 ), fabs ( val1 - 360 - val2 ) ), fabs ( val1 - val2 + 360 ) );

                //Compute difference
                T area = fabs ( ( key - key_old ) * val_diff );
                area_diff += area;
        }

        return area_diff;
}



template <typename T>
void TurningFunction::print ( std::ostream * output, const typename TTurningFunction <T>::Type & tf )
{
        typename TTurningFunction <T>::Type::const_iterator i_tf;

        for ( i_tf = tf.begin(); i_tf != tf.end(); ++i_tf )
        {
                *output << "key: " << i_tf->first << "\t" << "val: " << i_tf->second << std::endl;
        }

        *output << std::endl;
}


#endif
