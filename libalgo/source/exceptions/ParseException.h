// Description: Parse equation error class

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


#ifndef ParseException_H
#define ParseException_H


#include "Exception.h"


//Parser error: infix to postfix notation
class ParseException : public Exception
{
        protected:
		const std::string function_text;

        public:

		ParseException(const std::string & exception_text_, const std::string & function_text_) : Exception(exception_text_), function_text(function_text_) {}
		ParseException(const ParseException & e) : Exception(e.exception_text), function_text(e.function_text) {}
		virtual ~ParseException() throw() {}

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
			Exception::printException ( output );
                        *output << function_text << '\n';
                }

                virtual short getExceptionCode() const { return 20;}
};

#endif
