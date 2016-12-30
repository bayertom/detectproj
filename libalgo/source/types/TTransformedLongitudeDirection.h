// Description: Set a direction of the transformed longitude: positive value in east (normal) or west (reversed) direction
// from meridian passing through the cartographic pole

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


#ifndef TTransformedLongitudeDirection_H
#define TTransformedLongitudeDirection_H


//Set a direction of the transformed longitude: positive value in east (normal) or west (reversed) direction
//from meridian passing through the cartographic pole
typedef enum
{
	NoDirection = 0,	//Projection currently exists only in the normal aspect (ellipsoid -> sphere)
	NormalDirection,	//Positive direction lefts from the oriented line connecting the North Pole and transformed pole, measured at transformed pole
	ReversedDirection,	//Positive direction rights from the oriented line connecting the North Pole and transformed pole, measured at transformed pole
	NormalDirection2,	//Positive direction lefts from the  oriented line conncting the transformed pole and the North Pole, measured at transformed pole
	ReversedDirection2,	//Positive direction rights from the oriented line conncting the transformed pole and the North Pole, measured at transformed pole
	EquatorDirection,	//Positive values in the north direction from cartographic equator
} TTransformedLongitudeDirection;

#endif