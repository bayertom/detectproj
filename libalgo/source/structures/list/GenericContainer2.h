// Description: Abstract class representing generic container

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


#ifndef GenericContainer2_H
#define	GenericContainer2_H

#include <ostream>
#include <iostream>
#include <map>

#include "GenericContainer.h"

//Forward declaration
class TIndexList;

//Destructable or non-destructable container
typedef enum
{
        Destructable = 0,
        NonDestructable,
} TDestructable;


//Replacement of item numbers: ommit gaps in numbering to keep consistency of the container
typedef std::map < unsigned int, unsigned int> TReplaceItemsID;


//Forward declaration
template < typename T, const TDestructable destructable = Destructable >
class GenericContainer2;

//Abstract generic container for Items
template < typename Item, const TDestructable destructable >
class GenericContainer2 : public GenericContainer <Item>
{
        public:
                GenericContainer2 () : GenericContainer <Item> () {}
                GenericContainer2 ( const unsigned int n ) : GenericContainer <Item> () {for ( unsigned int i = 0; i < n; i++ ) this->items.push_back ( Item() );}

                template <typename Item2>
                GenericContainer2 ( const GenericContainer2 <Item2 > & source, const TIndexList & indices );

                virtual ~GenericContainer2() = 0;

        public:
                //Operators> safe version with .at()
                Item & operator [] ( const unsigned int index ) {return this->items.at ( index );}
                Item const & operator [] ( const unsigned int index ) const  {return this->items.at ( index );}

        public:
                //Member functions
                typename TItemsList <Item>::Type ::iterator insert ( typename TItemsList <Item>::Type ::iterator it_position, Item & p ) { return this->items.insert ( it_position, p );}
                typename TItemsList <Item>::Type ::iterator insert ( typename TItemsList <Item>::Type ::iterator it_position, const Item & p ) { return this->items.insert ( it_position, p );}

                template <typename InputIterator>
                void insert ( typename TItemsList <Item>::Type ::iterator it_position, InputIterator it_start, InputIterator it_end ) { this->items.insert ( it_position, it_start, it_end ); }

                Item & front() {return this->items.front();}
                Item const & front() const {return this->items.front();}
                Item & back() {return this->items.back();}
                Item const & back() const {return this->items.back();}
                void push_back ( Item & p ) {this->items.push_back ( p );}
                void push_back ( const Item & p ) {this->items.push_back ( p );}
                Item & at ( const unsigned int index ) {return this->items.at ( index );}
                Item const & at ( const unsigned int index ) const {return this->items.at ( index );}

                void clear();

        public:
                //Other functions
                virtual void print ( std::ostream * output = &std::cout ) const;
                void updateIDOfItems();
};


//Abstract generic container for Items*
template <typename Item, const TDestructable destructable>
class GenericContainer2 <Item *, destructable> : public GenericContainer <Item *>
{
        public:
                GenericContainer2 () : GenericContainer <Item *> () {}
                GenericContainer2 ( const unsigned int n ) : GenericContainer <Item *> ()
                {for ( unsigned int i = 0; i < n; i++ ) this->items.push_back ( new Item() );}

                GenericContainer2 ( const GenericContainer2 <Item *, destructable> & source );

                template <const TDestructable destructable2>
                GenericContainer2 ( const GenericContainer2 <Item *, destructable2> & source );

                template <typename Item2>
                GenericContainer2 ( const GenericContainer2 <Item2 *> & source, TIndexList & indices );

                GenericContainer2 <Item *, destructable > & operator = ( const GenericContainer2 <Item*, destructable > &source );

                template <const TDestructable destructable2>
                GenericContainer2 <Item *, destructable > & operator = ( const GenericContainer2 <Item*, destructable2 > &source );

                virtual ~GenericContainer2() = 0;

        public:
                //Operators: safe version with .at()
                Item * operator [] ( const unsigned int index ) {return this->items.at ( index );}
                Item * operator [] ( const unsigned int index ) const {return this->items.at ( index );}

        public:
                //Member functions
                typename TItemsList <Item*>::Type ::iterator insert ( typename TItemsList <Item*>::Type ::iterator it_position, Item * p ) { return this->items.insert ( it_position, p ); }
                typename TItemsList <Item*>::Type ::iterator insert ( typename TItemsList <Item*>::Type ::iterator it_position, const Item * p ) { return this->items.insert ( it_position, p ); }

                template <typename InputIterator>
                void insert ( typename TItemsList <Item*>::Type ::iterator it_position, InputIterator it_start, InputIterator it_end ) { this->items.insert ( it_position, it_start, it_end );}

                Item * front() {return this->items.front();}
                Item * front() const {return this->items.front();}
                Item * back() {return this->items.back();}
                Item * back() const {return this->items.back();}
                void push_back ( Item * p ) {this->items.push_back ( p );}
                void push_back ( const Item * p ) {this->items.push_back ( p );}
                Item * at ( const unsigned int index ) {return this->items.at ( index );}
                Item * at ( const unsigned int index ) const {return this->items.at ( index );}

                void clear();

        public:
                //Other functions
                virtual void print ( std::ostream * output = &std::cout ) const;

        protected:
                template <const TDestructable destructable2>
                void createCopy ( const GenericContainer2 <Item *, destructable2> &source );
                void updateIDOfItems();
};


template <typename T>
void operator << ( std::ostream & output, const GenericContainer2 <T> &gc )
{
        //Print container: overloaded operator, common for all derived class
        gc.print ( &output ) ;
}


#include "GenericContainer2.hpp"

#endif
