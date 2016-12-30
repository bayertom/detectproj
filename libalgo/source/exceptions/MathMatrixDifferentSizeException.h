// Description: Different size of matrices A, B in matrix algebra, throw exception

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


#ifndef MathMatrixDifferentSizeException_H
#define MathMatrixDifferentSizeException_H

#include "MathMatrixException.h"


//Matrix error: different size
template <typename TMatrix>
class MathMatrixDifferentSizeException: public MathMatrixException <TMatrix>
{
        protected:
                TMatrix M2;

        public:

                MathMatrixDifferentSizeException (const std::string & exception_text_, const std::string & function_text_, const TMatrix &M1_, const TMatrix &M2_ ) :
                        MathMatrixException <TMatrix> ( exception_text_, function_text_ , M1_ ), M2 ( M2_ ) {}

        public:
                virtual ~MathMatrixDifferentSizeException () throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        MathMatrixException <TMatrix>::printException ( output );
                        *output << "Matrix B, rows count: " << M2.rows() << ", cols count: " << M2.cols() << '\n';
                        M2.print ( output );
                }

                virtual TMatrix getArg() const {return M2;}
                virtual short getExceptionCode() const { return 13;}
};

#endif
