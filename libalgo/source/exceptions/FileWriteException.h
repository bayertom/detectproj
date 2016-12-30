// Description: Exception file write class

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


#ifndef FileWriteException_H
#define FileWriteException_H


#include "Exception.h"

//File write error
class FileWriteException : public Exception
{
        protected:
		std::string  file_text;

        public:

		FileWriteException(const std::string & exception_text_, const std::string & file_text_) : Exception(exception_text_), file_text(file_text_) {}
		FileWriteException(const FileWriteException & e) : Exception(e.exception_text), file_text(e.file_text) {}       
		virtual ~FileWriteException() throw() {}

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Exception::printException ( output );
                        *output << file_text << '\n';
                }

                virtual short getExceptionCode() const { return 7;}
};

#endif
