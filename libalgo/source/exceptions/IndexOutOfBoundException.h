// Description: Index ouf of bound, throw exception

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


#ifndef IndexOutOfBoundException_H
#define IndexOutOfBoundException_H

#include "Exception.h"


//Index out of bound error
class IndexOutOfBoundException : public Exception
{
        protected:

		std::string  index_text;

        public:

		IndexOutOfBoundException(const std::string & exception_text_, const std::string & index_text_) : Exception(exception_text_), index_text(index_text_) {}
		IndexOutOfBoundException(const IndexOutOfBoundException & e) : Exception(e.exception_text), index_text(e.index_text) {}
		virtual ~IndexOutOfBoundException() throw() {}

        public:


                virtual void printException ( std::ostream * output = &std::cout ) const
                {
			Exception::printException ( output );
                        *output << index_text << '\n';
                }

                virtual short getExceptionCode() const { return 8;}
};



#endif
