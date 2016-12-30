// Description: Compute tangent function

// Copyright (c) 2010 - 2011
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


#ifndef TangentFunction_HPP
#define TangentFunction_HPP


#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/line/HalfEdge.h"
#include "libalgo/source/structures/face/Face.h"

#include "libalgo/source/algorithms/angle3points/Angle3Points.h"
#include "libalgo/source/algorithms/faceperimeter/FacePerimeter.h"
#include "libalgo/source/algorithms/polylinelength/PolyLineLength.h"

#include "libalgo/source/exceptions/MathZeroDevisionException.h"


template <typename T>
T TangentFunction::compare2FacesUsingTangentFunction ( const Face <T> *f1, const Face <T> *f2, const TTangentFunctionRotationMethod & rotation_method,
	const TTangentFunctionScaleMethod & scale_method )
{
        //Compare 2 faces using difference of tangent functions
        //Result does not dependend on the order of vertices:  perform cyclic map rotation
        T rotation1, rotation2;
        typename TTangentFunction <T>::Type tf1, tf2;

        computeTangentFunctionFace ( f1, tf1, rotation1, rotation_method, scale_method );
        computeTangentFunctionFace ( f2, tf2, rotation2, rotation_method, scale_method );

        //Find minimum weight difference of tangent functions: unrotated map and cyclic rotated map
        //Compute initial difference of 2 tangent functions
        T difference = compareTangentFunctions <T> ( tf1, tf2, rotation_method );

        //Initial difference is zero, do not perform a cyclic rotation
        T min_difference = difference;
        if ( min_difference < ARGUMENT_ROUND_ERROR ) return min_difference;

        //Otherwise perform cyclic rotation of the first map according to the second map
        for ( unsigned int i = 0; i < tf1.size() ; i++ ) // perform n-1 cyclic rotations
        {
                //Remember first value in the second map
                T first_value = tf1.begin()->second;

                //Rotate all keys in the second map
                rotateKeys <T> ( tf1 );

                //Subtract first turning angle from all turning angles of the rotated map
                if ( rotation_method == RotationInvariant )
                {
                        for ( typename TTangentFunction <T>::Type ::iterator i_tf1 = tf1.begin();  i_tf1 != tf1.end(); ++i_tf1 )
                        {
                                i_tf1->second -= first_value;
                        }

                        //Compute last turning angle of the rotated map: angle[n-1] = angle[0] + angle[n-2]
                        tf1.rbegin()->second = ( ++tf1.rbegin() ) ->second + tf1.begin()->second;
                }

                //Copy first turning angle of the rotated map to the last one: angle[n-1] = angle[0]
                else
                {
                        tf1.rbegin()->second = tf1.begin()->second;
                }

                //Zero difference: stop computing
                if ( ( difference = compareTangentFunctions <T> ( tf1, tf2, rotation_method ) ) < ARGUMENT_ROUND_ERROR ) return difference;

                //Remember min value
                if ( difference < min_difference ) min_difference = difference;
        }

        //For rotation_method = RotationDependent do not perform cyclic rotation of the second map , return result
        if ( rotation_method == RotationDependent ) return min_difference;

        //Otherwise perform cyclic rotation of the second map according to the first
        for ( unsigned int i = 0; i < tf2.size() ; i++ ) // perform n-1 cyclic rotations
        {
                //Remember first value in the second map
                T first_value = tf2.begin()->second;

                //Rotate all keys in the second map
                rotateKeys <T> ( tf2 );

                //Subtract first turning angle from all turning angles of the rotated map
                for ( typename TTangentFunction <T>::Type ::iterator i_tf2 = tf2.begin();  i_tf2 != tf2.end(); ++i_tf2 )
                {
                        i_tf2->second -= first_value;
                }

                //Compute last turning angle of the rotated map: angle[n-1] = angle[0] + angle[n-2]
                tf2.rbegin()->second = ( ++tf2.rbegin() ) ->second + tf2.begin()->second;

                //Zero difference: stop computing
                if ( ( difference = compareTangentFunctions <T> ( tf1, tf2, rotation_method ) ) < ARGUMENT_ROUND_ERROR ) return difference;

                //Remember min value
                if ( difference < min_difference ) min_difference = difference;
        }

        return min_difference;
}


template <typename Point>
typename Point::Type TangentFunction::compare2PolyLinesUsingTangentFunction ( const Container <Point > *pl1, const Container <Point> *pl2, const TTangentFunctionRotationMethod & rotation_method,
	const TTangentFunctionScaleMethod & scale_method )
{
        //Compare 2 polylines using difference of tangent functions
        typename TTangentFunction <typename Point::Type>::Type tf1, tf2;

        //Compute tangent functions
        computeTangentFunctionPolyLine ( pl1 , tf1, rotation_method, scale_method );
        computeTangentFunctionPolyLine ( pl2 , tf2, rotation_method, scale_method );

        //Compute difference of tangent functions
        return compareTangentFunctions <typename Point::Type> ( tf1, tf2, rotation_method );
}


template <typename T>
void TangentFunction::computeTangentFunctionFace ( const Face <T> *face, typename TTangentFunction <T>::Type & tf, T & rotation, const TTangentFunctionRotationMethod & rotation_method,
	const TTangentFunctionScaleMethod & scale_method )
{
        //Compute normalized tangent function for the Face  given by HalfEdge
        T turning_angle_sum = 0, normalize_length_sum = 0;
        const HalfEdge <T> *e_start = face->getHalfEdge();
        HalfEdge <T> *e = const_cast <HalfEdge <T> *> ( e_start ) ;

        //Get triplet of points
        Node3DCartesian <T> *pi, *piii = NULL;
        Node3DCartesian <T> *pii = e_start->getPoint();

        //Jump short segments
        for ( pi = e->getPreviousEdge()->getPoint(); EuclDistance::getEuclDistance2D ( pi, pii ) < MIN_POSITION_DIFF ; e = e->getPreviousEdge(), pi = e->getPreviousEdge()->getPoint() ) {}

        //Get third point for initial edge
        for ( piii = e->getNextEdge()->getPoint(); EuclDistance::getEuclDistance2D ( pii, piii ) < MIN_POSITION_DIFF ; e = e->getNextEdge(), piii = e->getNextEdge()->getPoint() ) {}

        //Compute perimeter of the Face
        const T face_perimeter = FacePerimeter::getFacePerimeter ( e_start );

        //Throw exception
        if ( face_perimeter == 0 )
        {
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: tangent function, can not normalize function (1 / perimeter), ", "perimeter = 0." );
        }

        //Compute initial angle: rotation invariant
        if ( rotation_method == RotationInvariant )
        {
                turning_angle_sum += 180 - Angle3Points::getAngle3Points ( pi, pii, piii );
        }

        //Compute initial angle: rotation dependent (do not compute the sum)
        else
        {
                turning_angle_sum = 360 - Angle3Points::getAngle3Points ( &Node3DCartesian <T> ( pii->getX() - 1, pii->getY() ), pii, piii );
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
                if ( EuclDistance::getEuclDistance2D ( pii, piii ) > MIN_POSITION_DIFF )
                {
                        //Compute angle: rotation invariant
                        if ( rotation_method == RotationInvariant )
                        {
                                turning_angle_sum += 180 - Angle3Points::getAngle3Points ( pi, pii, piii );
                        }

                        //Compute initial angle: rotation dependent (do not compute the sum)
                        else
                        {
                                turning_angle_sum = 360 - Angle3Points::getAngle3Points ( &Node3DCartesian <T> ( pii->getX() - 1, pii->getY() ), pii, piii );
                        }

                        //Compute normalized length sum
                        normalize_length_sum += ( scale_method == ScaleDependent ? EuclDistance::getEuclDistance2D ( pi, pii ) / face_perimeter : EuclDistance::getEuclDistance2D ( pi, pii ) );

                        //Add turning angle computed for normalized length to map
                        tf[normalize_length_sum] = turning_angle_sum;

                        //Assign points
                        pi = pii;
                        pii = piii;
                }

                //Increment edge
                e = e->getNextEdge();

        } while ( e != e_start );

        //Add last value: we visit first point again
	const T last = ( scale_method == ScaleDependent ? 1.0 : face_perimeter );
        if ( rotation_method == RotationInvariant )
        {
		tf[ last] = turning_angle_sum + tf.begin()->second;
        }

        else
        {
                tf[last] = tf.begin()->second;
        }
}


template <typename Point>
void TangentFunction::computeTangentFunctionPolyLine ( const Container <Point> *points, typename TTangentFunction <typename Point::Type>::Type & tf, const TTangentFunctionRotationMethod & rotation_method,
	const TTangentFunctionScaleMethod & scale_method )
{
        //Compute normalized tangent function for the polyline
        const unsigned int n = points->size();
        typename Point::Type normalize_length_sum = 0, turning_angle_sum = 0;

        //Process all vertices
        const Point temp ( ( *points ) [0].getX() + 10, ( *points ) [0].getY() );

        //Get polyle=ine length
        const typename Point::Type polyline_length = PolyLineLength::getPolyLineLength ( points );

        //Compute initial angle: rotation invariant
        if ( rotation_method == RotationInvariant )
        {
                turning_angle_sum += 180 - Angle3Points::getAngle3Points ( & ( *points ) [0], & ( *points ) [1],  & ( *points ) [2] );
        }

        //Compute initial angle: rotation dependent
        else
        {
                turning_angle_sum = 360 - Angle3Points::getAngle3Points ( &Point ( ( *points ) [1].getX() - 1, ( *points ) [1].getY() ),  & ( *points ) [1], & ( *points ) [2] );
        }

        //Compute normalized length
        normalize_length_sum += EuclDistance::getEuclDistance2D ( & ( *points ) [1],  & ( *points ) [2] ) / polyline_length;

        //Add turning angle computed for normalized length to map
        tf[normalize_length_sum] = turning_angle_sum;

        //Process all inernal vertices of the polyline
        for ( unsigned int i = 2; i < n - 1; i++ )
        {
                //Compute angle: rotation invariant
                if ( rotation_method == RotationInvariant )
                {
                        turning_angle_sum += 180 - Angle3Points::getAngle3Points ( & ( *points ) [i-1], & ( *points ) [i],  & ( *points ) [i+1] );
                }

                else
                {
                        turning_angle_sum = 360 - Angle3Points::getAngle3Points ( &Point ( ( *points ) [i].getX() - 1, ( *points ) [i].getY() ),  & ( *points ) [i], & ( *points ) [i+1] );
                }

                //Compute normalized length
		normalize_length_sum += ( scale_method == ScaleDependent ? EuclDistance::getEuclDistance2D ( & ( *points ) [i],  & ( *points ) [i+1] ) / polyline_length : EuclDistance::getEuclDistance2D ( & ( *points ) [i],  & ( *points ) [i+1] ) );

                //Add turning angle computed for normalized length to map
                tf[normalize_length_sum] = turning_angle_sum;
        }

        //Add last item
	const typename Point::Type last = ( scale_method == ScaleDependent ? 1.0 : polyline_length );
        tf[last] = turning_angle_sum;
}


template <typename T>
T TangentFunction::compareTangentFunctions ( typename TTangentFunction <T>::Type & tf1, typename TTangentFunction <T>::Type & tf2, const TTangentFunctionRotationMethod & rotation_method )
{
        //Compare two normalized tangent functions for Faces p1 and p2 and returns difference ratio
        T key = 0, key_old = 0, value_1  = 0, value_2 = 0, area_diff = 0;
        const unsigned int n = tf1.size();
        const unsigned int m = tf2.size();

        //Iterators
        typename TTangentFunction <T>::Type ::iterator i_tf1 = tf1.begin();
        typename TTangentFunction <T>::Type ::iterator i_tf2 = tf2.begin();

        //Initialize area_diff
        area_diff = 0;

        //Merge both maps and compute differencies
        for ( unsigned int i = 0; i < m + n; i++ )
        {
                //Compute angles difference
                T val_diff =  fabs ( value_2 - value_1 );

                //Rotation dependent algorithm: find minimum value of both angles ( to avoid 359 deg  and 1deg difference, etc)
                if ( rotation_method == RotationDependent )
                {
                        val_diff = min ( min ( val_diff, fabs ( value_2 - 360 - value_1 ) ), fabs ( value_2 - value_1 + 360 ) );
                }

                //Last element in the first map
                if ( i_tf1 == tf1.end() )
                {
                        //Compute difference: using small heuristic
                        area_diff += fabs ( ( i_tf2->first - key_old ) * val_diff );

                        //Assign values
                        value_1 = i_tf2->second;
                        key_old = i_tf2->first;

                        //Increment iterator
                        i_tf2 ++;

                        //Jump to the next iteration
                        continue;
                }

                //Last element in the second map
                if ( i_tf2 == tf2.end() )
                {
                        //Compute difference
                        area_diff += fabs ( ( i_tf1->first - key_old ) * val_diff );

                        //Assign old values
                        value_2 = i_tf1->second;
                        key_old = i_tf1->first;

                        //Increment iterator
                        i_tf1 ++;

                        //Jump to the next iteration
                        continue;
                }

                //Map 1 value < Map 2 value
                if ( i_tf1->first < i_tf2 -> first )
                {
                        //Compute difference
                        if ( i_tf2 != tf2.begin() )
                        {
                                area_diff += fabs ( ( i_tf1->first - key_old ) * val_diff );
                        }

                        //Assign values
                        value_1 = i_tf1->second;
                        key_old = i_tf1->first;

                        //Increment iterator
                        i_tf1 ++;
                }

                //Map 2 value < Map 1 value
                else
                {
                        //Compute difference
                        if ( i_tf1 != tf1.begin() )
                        {
                                area_diff += fabs ( ( i_tf2->first - key_old ) * val_diff );
                        }

                        //Assign values
                        value_2 = i_tf2->second;
                        key_old = i_tf2->first;

                        //Increment iterator
                        i_tf2 ++;
                }
        }

        return area_diff;
}


template <typename T>
void TangentFunction::rotateKeys ( typename TTangentFunction <T>::Type & tf )
{
        //Cyclic rotation or the keys in map
        T first_value = tf.begin()->second; //Remember first value

        //Erase first element
        tf.erase ( tf.begin() );

        //Remeber second key (non zero key)
        T second_key = tf.begin()->first;

        //Change all keys of the actual map
        for ( typename TTangentFunction <T>::Type ::iterator i_tf = tf.begin(); i_tf != tf.end(); ++i_tf )
        {
                typename TTangentFunction <T>::Type ::key_type * actual_key = const_cast <typename TTangentFunction <T>::Type ::key_type *> ( &i_tf->first ) ;
                *actual_key = i_tf->first - second_key;
        }

        //Add new element at the end of the map: map is rotated
        tf[ tf.rbegin()->first + second_key ] = first_value;
}

#endif
