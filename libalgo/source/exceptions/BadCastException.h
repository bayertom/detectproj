// Description: Error bad casting class

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


#ifndef BadCastException_H
#define BadCastException_H


#include "Exception.h"

//Bad cast error
class BadCastException: public Exception
{
        protected:
                std::string type_text;

        public:
		BadCastException(const std::string  &exception_text_, const std::string  & type_text_) : Exception(exception_text_), type_text(type_text_) {}
		BadCastException(const BadCastException & e) : Exception(e.exception_text), type_text(e.type_text) {}
		virtual ~BadCastException() throw() {}
               

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
			Exception::printException ( output );
                        *output << type_text << '\n';
                }

                virtual short getExceptionCode() const { return 3;}
};

#endif
