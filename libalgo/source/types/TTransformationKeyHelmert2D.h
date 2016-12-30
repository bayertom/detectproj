// Description: Helmert 2D transformation key

// Copyright (c) 2015 - 2016
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
// along with this library. If not, see <http://www.gnu.org/licenses/>


#ifndef TTransformationKeyHelmert2D_H
#define TTransformationKeyHelmert2D_H


template <typename T>
struct TTransformationKeyHelmert2D
{
        T x_mass_local, y_mass_local;		//Centre of mass: local system
        T x_mass_global, y_mass_global;		//Centre of mass: global system
        T c1, c2;				//Both transformation coefficients
        T J;					//(X^T * X)^-1
        T k;					//Sum of weights

        TTransformationKeyHelmert2D () : x_mass_local ( 0.0 ), y_mass_local ( 0.0 ), x_mass_global ( 0.0 ), y_mass_global ( 0.0 ), c1 ( 1.0 ), c2 ( 1.0 ), J ( 0.0 ), k ( 0.0 ) {}
};


#endif