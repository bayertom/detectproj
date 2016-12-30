// Description: General exception class, other error classes derived from this class

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


#ifndef Exception_H
#define Exception_H


#include <exception>
#include <ostream>
#include <iostream>
#include <string>


//Exception, inherited from std::exception
class Exception : public std::exception
{
        protected:
                std::string exception_text;

        public:
                Exception() : exception_text ( "" ) {}
		Exception(const std::string & exception_text_) : exception_text(exception_text_) {}
		Exception(const Exception & e) : exception_text(e.exception_text) {}
                virtual ~Exception () throw() { }

        public:
                const std::string &getExceptionText() const {return exception_text;}

        public:
                virtual void printException ( std::ostream * output = &std::cout ) const { *output << exception_text << '\n'; }
                virtual short getExceptionCode() const { return 1;}
};


#endif

