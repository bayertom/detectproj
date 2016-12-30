// Description: General error in matrix algebra, other matrix error classes derived from this class

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


#ifndef MathMatrixException_H
#define MathMatrixException_H

#include "MathException.h"


//Math matrix error
template <typename TMatrix>
class MathMatrixException : public MathException <TMatrix>
{
        protected:
                TMatrix M;

        public:

                MathMatrixException (const std::string & exception_text_, const std::string & function_text_, const TMatrix &M_ ) :
                        MathException <TMatrix> ( exception_text_, function_text_ ), M ( M_ ) {}

        public:
                virtual ~MathMatrixException() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        MathException <TMatrix>::printException ( output );
                        *output << "Matrix A, rows count: " << M.rows() << ", cols count: " << M.cols() << '\n';
                        M.print ( output );
                }

                virtual TMatrix getArg( ) const {return M;}
                virtual short getExceptionCode() const { return 12;}

};

#endif
