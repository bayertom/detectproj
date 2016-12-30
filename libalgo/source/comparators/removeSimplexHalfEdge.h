#ifndef removeSimplexHalfEdge_H
#define removeSimplexHalfEdge_H

#include "libalgo/source/structures/line/HalfEdge.h"

//Remove simplex half edges
class removeSimplexHalfEdge
{
        public:
                removeSimplexHalfEdge() {};

                template <typename T>
                bool operator() ( const HalfEdge <T> *h ) const
                {
                        //Set simplex half edge to be deleted
                        if ( h->isSimplexEdge() )
                        {
                                delete h;				//Call destructor
                                h = NULL;				//Initialize to null

                                return 1;				//Mark item as deleteable
                        }

                        else
                        {
                                return 0;				//Mark item as undeleteable
                        }
                }
};

#endif

