// Description: Matrix error: matrix not positive definite (throw exception)

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

#ifndef MathMatrixNotPositiveDefiniteException_H
#define MathMatrixNotPositiveDefiniteException_H

#include "MathMatrixException.h"


//Matrix error: matrix not positive definite
template <typename TMatrix>
class MathMatrixNotPositiveDefiniteException: public MathMatrixException <TMatrix>
{
        public:

                MathMatrixNotPositiveDefiniteException(const std::string & exception_text_, const std::string & function_text_, const TMatrix &M_ )
                        : MathMatrixException <TMatrix> ( exception_text_, function_text_, M_ ) {}

        public:
                virtual ~MathMatrixNotPositiveDefiniteException() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        MathMatrixException <TMatrix>::printException ( output );
                        *output << "Matrix A, rows count: " << this->M.rows() << ", cols count: " << this->M.cols() << '\n';
                }

                virtual TMatrix getArg() const { return MathMatrixException<TMatrix>::getArg(); }
                virtual short getExceptionCode() const { return 21;}
};

#endif
