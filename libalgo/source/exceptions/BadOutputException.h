// Description: Error bad Output class

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


#ifndef BadOutputException_H
#define BadOutputException_H


#include "Exception.h"


//Bad Output error
class BadOutputException : public Exception
{
        protected:
		std::string output_text;

        public:

		BadOutputException(const std::string & exception_text_, const std::string & output_text_) : Exception(exception_text_), output_text(output_text_) {}
		BadOutputException(const BadOutputException & e) : Exception(e.exception_text), output_text(e.output_text) {}
                virtual ~BadOutputException() throw(){}

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
			Exception::printException();
                        *output << output_text << '\n';
                }

                virtual short getExceptionCode() const  { return 5;}
};


#endif
