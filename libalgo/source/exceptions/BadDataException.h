// Description: Bad input data, throw exception

// Copyright (c) 2010 - 2016
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


#ifndef BadDataException_H
#define BadDataException_H

#include "Exception.h"


//Bad data exception
class BadDataException : public Exception
{
        protected:
                std::string data_text;

        public:

		BadDataException(const std::string  & exception_text_, const std::string  & data_text_) : Exception(exception_text_), data_text(data_text_) {}
		BadDataException(const BadDataException & e) : Exception(e.exception_text), data_text(e.data_text) {}
		virtual ~BadDataException() throw() {}

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Exception::printException ( output );
                        *output << data_text << '\n';
                }

                virtual short getExceptionCode() const  { return 4;}
};


#endif
