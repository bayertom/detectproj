// Description: Mutation strategy used in the differential evolution

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


#ifndef TMutationStrategy_H
#define TMutationStrategy_H


typedef enum
{
        DERand1Strategy = 0,		//DE/Rand/1 mutation strategy
	DERand2Strategy,		//DE/Rand/2 mutation strategy
	DERandDir1Strategy,		//DE/Rand/Dir/1 mutation strategy
	DERandDir2Strategy,		//DE/Rand/Dir/2 mutation strategy
	DERandBest1Strategy,		//DE/Rand/Best/1 mutation strategy
        DERandBest2Strategy,		//DE/Rand/Best/2 mutation strategy
	DERandBestDir1Strategy,		//DE/Rand/Best/Dir/1 mutation strategy
	DETargetToBest1Strategy,	//DE/Target-to-best/Best/1 mutation strategy
	SACPStrategy,			//Self adapting control of parameters
} TMutationStrategy;

#endif