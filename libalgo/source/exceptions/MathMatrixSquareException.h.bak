// Description: Square matrix error class

// Copyright (c) 2010 - 2012
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

#ifndef ErrorMathMatrixSquare_H
#define ErrorMathMatrixSquare_H


#include "ErrorMathMatrix.h"


//Matrix error: matrix not square
template <typename TMatrix>
class ErrorMathMatrixSquare : public ErrorMathMatrix <TMatrix>
{
        public:

                ErrorMathMatrixSquare ( const char * exception_text_, const char * function_text_, unsigned int rows_count_, unsigned int columns_count_ )
                        : ErrorMathMatrix <TMatrix> ( exception_text_, function_text_, rows_count_, columns_count_ ) {}

        public:
                virtual ~ErrorMathMatrixSquare() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        ErrorMath <TMatrix>::printException ( output );
                        *output << "Matrix B, rows count: " << this->rows_count << ", cols count: " << this->columns_count << '\n';
                }

                virtual short getErrorCode() const { return 15;}
};

#endif
