// Description: Abstract class of the container with defined clear method

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


#ifndef GenericContainer2_HPP
#define	GenericContainer2_HPP

#include <algorithm>
#include <fstream>

#include "libalgo/source/structures/list/IndexLists.h"

template <typename Item, const TDestructable destructable>
GenericContainer2 <Item, destructable>::~GenericContainer2()
{
        //Destructor (Generic container)
        clear();
}


template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item, destructable>::clear()
{
        //Clear items (Generic container)
        this->items.clear();
}


template <typename Item, const TDestructable destructable>
template <typename Item2>
GenericContainer2 <Item, destructable> :: GenericContainer2 ( const GenericContainer2 <Item2> & source, const TIndexList & indices )
{
        //Create copy of selected items given by indices list from container
        for ( unsigned int i = 0; i < indices. size(); i++ )
        {
                //Get index
                const unsigned int index = indices [i];

                //Add point to the list
                this->items.push_back ( static_cast <Item &> ( * source [index] ) );
        }
}


template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item, destructable>::updateIDOfItems()
{
        //Renumber all items: ommit spaces in items numbers
        const unsigned int n = this->items.size();

        //Create container for replacement of item numbers
        TReplaceItemsID replace_items_id;

        if ( n > 0 )
        {
                //Create pairs: <item_id, index>
                for ( unsigned int i = 0; i < n; i++ ) replace_items_id [this->items[i].getPointID() ] = i;

                //Replace item numbers to ommit (possible) gaps
                for ( TReplaceItemsID::const_iterator i_replace_items_id = replace_items_id.begin();
                                i_replace_items_id != replace_items_id.end(); i_replace_items_id ++ )
                {
                        this->items[i_replace_items_id->second].updateID();
                }
        }
}


template <typename Item, const TDestructable destructable>
GenericContainer2 <Item *, destructable>::GenericContainer2 ( const GenericContainer2 <Item *, destructable> & source )
{
        //Copy constructor, deep / shallow copy (Explicit specialization for Item *): same default templates parameters
        if ( this != &source )
        {
                createCopy ( source );
        }
}


template <typename Item, const TDestructable destructable>
template <const TDestructable destructable2>
GenericContainer2 <Item *, destructable>::GenericContainer2 ( const GenericContainer2 <Item *, destructable2> & source )
{
        //Copy constructor, deep / shallow copy (Explicit specialization for Item *): different default templates parameters
        createCopy ( source );
}


template <typename Item, const TDestructable destructable>
GenericContainer2 <Item *, destructable > & GenericContainer2 <Item *, destructable>::operator = ( const GenericContainer2 <Item*, destructable > &source )
{
        //Operator = , deep / shallow copy (Explicit specialization for Item *) : same default templates parameters
        if ( this != &source )
        {
                //Create copy (deep or shallow)
                createCopy ( source );
        }

        return *this;
}


template <typename Item, const TDestructable destructable>
template <const TDestructable destructable2>
GenericContainer2 <Item *, destructable > & GenericContainer2 <Item *, destructable>::operator = ( const GenericContainer2 <Item*, destructable2 > &source )
{
        //Operator = , deep / shallow copy (Explicit specialization for Item *) : different default templates parameters
        createCopy ( source );

        return *this;
}


template <typename Item, const TDestructable destructable>
GenericContainer2 <Item *, destructable>::~GenericContainer2()
{
        //Destructor (Explicit specialization for Item *)
        clear();
}


template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item *, destructable>::updateIDOfItems()
{
        //Renumber all itemss: ommit gaps in items numbers (Partial specialization for Item*)
        const unsigned int n = this->items.size();

        //Create container for replacement of item numbers
        TReplaceItemsID replace_items_id;

        if ( n > 0 )
        {
                //Create pairs: <item_id, index>
                for ( unsigned int i = 0; i < n; i++ ) replace_items_id [this->items[i]->getPointID() ] = i;

                //Replace numebrs numbers with actual item_id to ommit (possible) gaps
                for ( TReplaceItemsID::const_iterator i_replace_items_id = replace_items_id.begin();
                                i_replace_items_id != replace_items_id.end(); i_replace_items_id ++ )
                {
                        this->items[i_replace_items_id->second] ->updateID();
                }
        }
}


template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item *, destructable>::clear()
{
        //Clear all items (Explicit specialization for Item *)
        const int n = this->items.size();

        if ( n > 0 )
        {
                //Call destructor for each item
                if ( destructable == Destructable )
                {
                        for ( typename TItemsList <Item*>::Type ::iterator i_items = this->items.begin(); i_items != this->items.end(); ++i_items )
                        {
                                //Delete each node
                                if ( * ( i_items ) != NULL )
                                {
                                        delete *i_items;
                                        ( *i_items ) = 0;
                                }
                        }
                }

                //Clear list
                this->items.clear();
        }
}



template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item, destructable>::print ( std::ostream * output ) const
{
        //Print all items (Generic container)
        *output << '\n';

        for ( typename TItemsList <Item>::Type ::const_iterator i_items = this->items.begin(); i_items != this->items.end(); ++i_items )
        {
                //*output << '\n' << "Item #" << ( ( i_items - this->items.begin() ) + 1 );
                ( *i_items ).print ( output );
        }
}


template <typename Item, const TDestructable destructable>
void GenericContainer2 <Item *, destructable>::print ( std::ostream * output ) const
{
        //Print all items (Explicit specialization for Item *)
        *output << '\n';

        for ( typename TItemsList <Item*>::Type ::const_iterator i_items = this->items.begin(); i_items != this->items.end(); ++i_items )
        {
                ( *i_items )->print ( output );
        }
}


template <typename Item, const TDestructable destructable>
template <const TDestructable destructable2>
void GenericContainer2 <Item *, destructable> :: createCopy ( const GenericContainer2 <Item *, destructable2> &source )
{
        //Create deep copy of the container (Explicit specialization for Item *)
        typename TItemsList <Item *>::Type ::const_iterator i_item_source;

        //Clear items
        clear();

        //Create new nodes as the copy of old nodes
        for ( i_item_source = source.begin(); i_item_source != source.end(); ++i_item_source )
        {
                //Create deep copy of the old item
                if ( destructable == Destructable )
                {
                        Item * temp = new Item ( *i_item_source );

                        //Add new node to the container
                        this->items.push_back ( temp );
                }

                //Create shallow copy of the old item
                else
                {
                        this->items.push_back ( *i_item_source );
                }
        }
}


template <typename Item, const TDestructable destructable>
template <typename Item2>
GenericContainer2 <Item *, destructable> :: GenericContainer2 ( const GenericContainer2 <Item2 *> & source, TIndexList & indices )
{
        //Create copy of selected items given by indices list from container (Explicit specialization for Item* )
        for ( unsigned int i = 0; i < indices. size(); i++ )
        {
                //Create index
                const unsigned int index = indices [i];

                //Create new Item
                if ( destructable )
                {
                        this->items.push_back ( new Item ( static_cast <Item2 *> ( source[index] ) ) );
                }

                //Copy pointers
                else
                {
                        this->items.push_back ( static_cast <Item2 *> ( source[index] ) );
                }
        }
}


#endif
