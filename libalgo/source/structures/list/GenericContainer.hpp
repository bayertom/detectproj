// Description: Abstract class of the container

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


#ifndef GenericContainer_HPP
#define GenericContainer_HPP

#include <algorithm>

template <typename Item>
GenericContainer <Item>::~GenericContainer() {}


template <typename Item>
template <typename CompSort>
void GenericContainer <Item>::sort ( typename TItemsList <Item>::Type ::iterator it_begin, typename TItemsList <Item>::Type ::iterator it_end, CompSort comp_sort )
{
        //Sort items
        std::sort ( it_begin, it_end, comp_sort );
}


template <typename Item>
template <typename CompSort, typename CompEqual>
void GenericContainer <Item>::removeDuplicateElements ( typename TItemsList <Item>::Type ::iterator it_begin, typename TItemsList <Item>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_equal )
{
        //Remove duplicate items from the Generic Container
        std::sort ( it_begin, it_end, comp_sort );
        typename TItemsList <Item>::Type ::iterator i_new_end = std::unique ( it_begin, it_end, comp_equal );
        this->items.erase ( i_new_end, items.end() );
}


#endif
