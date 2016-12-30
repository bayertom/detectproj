// Description: General math error class (other math error classes are derived from this class), throw exception

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


#ifndef MathException_H
#define MathException_H

#include "Exception.h"


//Math error
template <typename T>
class MathException : public Exception
{
        protected:
		const std::string math_text;

        public:

		MathException(const std::string & exception_text_, const std::string & math_text_) : Exception(exception_text_), math_text(math_text_) {}
		MathException(const MathException & e) : Exception(e.exception_text), math_text(e.math_text) {}
 		virtual ~MathException() throw() {}

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Exception::printException ( output );
                        *output << math_text << '\n';
                }

                virtual T getArg() const = 0;

                virtual short getExceptionCode() const { return 9;}
};


#endif
