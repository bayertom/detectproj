// Description: Class representing a container of items

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

#ifndef Container_H
#define Container_H

#include "GenericContainer2.h"
#include "Container.h"

//Forward declarations
class TIndexList;

template < typename T>
class sortPointsByID;

template <typename T>
class Point3DGeographic;

template <typename T>
class Node3DCartesian;

template <typename T>
class Edge;

template <typename T>
class Projection;

template <typename T>
class Face;

template <typename T>
class VoronoiCell;

template < typename T, const TDestructable destructable = Destructable >
class Container;

//Dimension of loaded points: 2D or 3D
typedef enum
{
        Dim2D = 0,
        Dim3D,
} TDimension;


//Structure for Tag Dispatching: emulate Partial specialization
template <const TDimension dimension>
struct Dimension {};


//Container for points ands its partial specializations for other items
template <typename Point, const TDestructable destructable>
class Container : public GenericContainer2 <Point, destructable>
{
        public:
                Container() : GenericContainer2 <Point, destructable>() {}

                template <typename Point2>
                Container ( const Container <Point2> & source, const TIndexList & indices ) : GenericContainer2  <Point> ( source, indices ) {}

                Container ( const Container <Point> & source ) : GenericContainer2  <Point, destructable> ( source ) {}

                virtual ~Container() {}

        public:
                template <TDimension dim>
                void load ( const char * file, const bool print_exception = true, std::ostream * output = &std::cout );

                void loadFromVector ( const std::vector<Point>& data );

                void toIndexList ( TIndexList & il );

        private:
                template <TDimension dim>
                void loadPoints ( const char * file, Dimension <dim>, const bool print_exception, std::ostream * output );

                void loadPoints ( const char * file, Dimension <Dim2D>, const bool print_exception, std::ostream * output );

        protected:
                void updateIDOfItems();

        public:
                //Other functions
                template <typename CompSort, typename CompEqual>
                void removeDuplicateElements ( typename TItemsList <Point>::Type ::iterator it_begin, typename TItemsList <Point>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_equal );
};


//Partial specialization for pointers (Points *)
template <typename Point, const TDestructable destructable>
class Container <Point *, destructable> : public GenericContainer2 <Point *, destructable>
{
        public:
                Container() : GenericContainer2 <Point *, destructable> () {}
                Container ( const unsigned int n ) : GenericContainer2 <Point *, destructable> ( n ) {}

                template <const TDestructable destructable2>
                Container ( const Container <Point*, destructable2 > &source ) : GenericContainer2  <Point*, destructable> ( source ) {}

                template <typename Point2>
                Container ( const Container <Point2> & source, const TIndexList & indices ) : GenericContainer2  <Point*, destructable> ( source, indices ) {}

                Container <Point *, destructable > & operator = ( const Container <Point *, destructable > &source );

                template <const TDestructable destructable2>
                Container <Point *, destructable > & operator = ( const Container <Point *, destructable2 > &source );

                virtual ~Container() {}

        public:
                template <TDimension dim>
                void load ( const char * file, const bool print_exception = true, std::ostream * output = &std::cout );
                
                void loadFromVector ( const std::vector<Point*>& data );

                void toIndexList ( TIndexList & il );

        private:
                template <TDimension dim>
                void loadPoints ( const char * file, Dimension <dim>, const bool print_exception, std::ostream * output );

                void loadPoints ( const char * file, Dimension <Dim2D>, const bool print_exception, std::ostream * output );                

        public:
                //Other functions
                template <typename CompSort, typename CompEqual>
                void removeDuplicateElements ( typename TItemsList <Point *>::Type ::iterator it_begin, typename TItemsList <Point *>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_equal );
};


//Partial specialization for Edge: storing edges
template <typename Point, const TDestructable destructable>
class Container <Edge <Point> , destructable> : public GenericContainer2 <Edge <Point>, destructable>
{
        public:
                Container () : GenericContainer2 <Edge <Point>, destructable> () {}


        public:
                //Other functions
                void load ( const char * file, Container <Point, destructable> *pl, const bool print_exception = true );
};


//Partial specialization for Face: storing faces
template <typename T, const TDestructable destructable>
class Container <Face <T> *, destructable> : public GenericContainer2 <Face <T> *, destructable>
{
        public:
                Container () : GenericContainer2 <Face <T> *, destructable> () {}

                template <const TDestructable destructable2>
                Container ( const Container <Face <T>*, destructable2 > &source ) : GenericContainer2  <Face <T>*, destructable> ( source ) {}

                Container <Face <T> *, destructable > & operator = ( const Container <Face <T> *, destructable > &source );

                virtual ~Container() {}
};


//Partial specialization for Voronoi cell: storing Voronoi Cells
template <typename T, const TDestructable destructable>
class Container <VoronoiCell <T> *, destructable> : public GenericContainer2 <VoronoiCell <T> *, destructable>
{
        public:
                Container () : GenericContainer2 <VoronoiCell <T> *, destructable> () {}

                template <const TDestructable destructable2>
                Container ( const Container <VoronoiCell <T>*, destructable2 > &source ) : GenericContainer2  <VoronoiCell <T>*, destructable> ( source ) {}

                Container <VoronoiCell <T> *, destructable > & operator = ( const Container <VoronoiCell <T> *, destructable > &source );

                virtual ~Container() {}
};


//Container for cartographic projections: Partial specialization for pointers (Projection *): storing cartographic projection
template <typename T, const TDestructable destructable>
class Container <Projection <T> *, destructable> : public GenericContainer2 <Projection <T> *, destructable>
{
        public:
                Container() : GenericContainer2 <Projection<T> *, destructable> () {}

                virtual ~Container() {}

        public:
                virtual void load ( const char * current_dir, const bool print_exception = true );

};

#include "Container.hpp"

#endif
