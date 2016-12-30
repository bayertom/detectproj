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


#ifndef GenericContainer_H
#define	GenericContainer_H

#include <vector>
#include <iterator>
#include <algorithm>

//New user defined type: list of dynamically allocated Items
template <typename Item>
struct TItemsList
{
        typedef std::vector <Item> Type;
};


//Abstract generic container for Items
template <typename Item>
class GenericContainer
{
        protected:
                typename TItemsList <Item>::Type items;

        public:
                //Store type of the object
                typedef Item Type;

        public:
                GenericContainer() : items ( 0 ) {}
                virtual ~GenericContainer() = 0;

        public:
                //Iterators
                typename TItemsList <Item>::Type ::iterator begin() { return this->items.begin(); }
                typename TItemsList <Item>::Type ::const_iterator begin() const { return this->items.begin(); }
                typename TItemsList <Item>::Type ::iterator end() { return this->items.end(); }
                typename TItemsList <Item>::Type ::const_iterator end() const { return this->items.end(); }
                std::back_insert_iterator <typename TItemsList <Item>::Type> back_inserter() {return std::back_inserter ( items ); }

        public:
                //Member functions
                void pop_back() {items.pop_back();}
                unsigned int size() const {return items.size();}
                typename TItemsList <Item>::Type ::iterator erase ( typename TItemsList <Item>::Type ::iterator i_position ) {return items.erase ( i_position );}
                bool empty () const {return items.empty();}
                void erase ( typename TItemsList <Item>::Type ::iterator i_start, typename TItemsList <Item>::Type ::iterator i_end ) {items.erase ( i_start, i_end );}
                void pop_front() {items.erase ( items.begin() );}

                template <typename InputIt, typename OutputIt>
                OutputIt copy ( InputIt it_first, InputIt it_last, OutputIt it_result )  { return std::copy ( it_first, it_last, it_result ); }
                //typename TItemsList <Item>::Type ::iterator copy (typename TItemsList <Item>::Type ::iterator it_first, typename TItemsList <Item>::Type ::iterator it_last, typename TItemsList <Item>::Type ::iterator it_result ) ;

                template <typename InputIt, typename OutputIt>
                OutputIt reverse_copy ( InputIt it_first, InputIt it_last, OutputIt it_result ) { return std::reverse_copy ( it_first, it_last, it_result ); }
                //typename TItemsList <Item>::Type ::iterator reverse_copy (typename TItemsList <Item>::Type ::iterator it_first, typename TItemsList <Item>::Type ::iterator it_last, typename TItemsList <Item>::Type ::iterator it_result ) ;

        public:
                template <typename CompSort>
                void sort ( typename TItemsList <Item>::Type ::iterator it_begin, typename TItemsList <Item>::Type ::iterator it_end, CompSort comp_sort );

        public:
                //Other functions
                template <typename CompSort, typename CompEqual>
                void removeDuplicateElements ( typename TItemsList <Item>::Type ::iterator it_begin, typename TItemsList <Item>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_eqaul );
};


#include "GenericContainer.hpp"

#endif
