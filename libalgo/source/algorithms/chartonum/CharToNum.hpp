// Description: Convert char to number: better performance than atof()

// Copyright (c) 2010 - 2015
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


#ifndef CharToNum_HPP
#define CharToNum_HPP

#include <limits> 

template <typename T>
T CharToNum::atof2(const char *num)
{
	//Improvement of the atof() function
	//Original source from http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/
	// Tino Didriksen, 2014

	// Skip leading white space, if any.
	while (white_space(*num)) 
	{
		num += 1;
	}

	T r = 0.0;
	unsigned int c = 0; // counter to check how many numbers we got!

	// Get the sign!
	bool neg = false;
	if (*num == '-') 
	{
		neg = true;
		++num;
	}

	else if (*num == '+'){
		neg = false;
		++num;
	}

	// Get the digits before decimal point
	while (valid_digit(*num)) 
	{
		r = (r * 10.0) + (*num - '0');
		++num; 
		++c;
	}

	// Get the digits after decimal point
	if (*num == '.') {
		T f = 0.0;
		T scale = 1.0;
		++num;

		while (*num >= '0' && *num <= '9') 
		{
			f = (f*10.0) + (*num - '0');
			++num;
			scale *= 10.0;
			++c;
		}
		r += f / scale;
	}

	// FIRST CHECK:
	if (c == 0)
	{ 
		// We got no dezimal places! this cannot be any number!
		throw ("BadDataException: the strinh has no dezimal place.", "Conversion to number aborted.");

	} 

	// Get the digits after the "e"/"E" (exponenet)
	if (*num == 'e' || *num == 'E'){
		unsigned int e = 0;

		bool negE = false;
		++num;
		if (*num == '-') 
		{
			negE = true;
			++num;
		}

		else if (*num == '+')
		{
			negE = false;
			++num;
		}

		// Get exponent
		c = 0;
		while (valid_digit(*num)) 
		{
			e = (e * 10) + (*num - '0');
			++num; ++c;
		}

		if (!neg && e > std::numeric_limits<T>::max_exponent10)
		{
			e = std::numeric_limits<T>::max_exponent10;
		}

		else if (e < std::numeric_limits<T>::min_exponent10)
		{
			e = std::numeric_limits<T>::max_exponent10;
		}

		// SECOND CHECK:
		if (c == 0)
		{ 
			// We got no  exponent! this was not intended!!
			throw ("BadDataException: the string has no exponent.", "Conversion to number aborted.");
		} 

		T scaleE = 1.0;
		// Calculate scaling factor.

		while (e >= 50) 
		{ 
			scaleE *= 1E50; 
			e -= 50; 
		}

		while (e > 0) 
		{ 
			scaleE *= 10.0; 
			e -= 1; 
		}

		if (negE)
		{
			r /= scaleE;
		}

		else
		{
			r *= scaleE;
		}
	}

	// POST CHECK:
	// skip post whitespaces
	while (white_space(*num))
	{
		++num;
	}

	if (*num != '\0')
	{ 
		// If next character is not the terminating character
		throw ("BadDataException: the next character is not \n.", "Conversion to number aborted.");
	} 

	// Apply sign to number
	if (neg)
	{ 
		r = -r; 
	}

	return r;
}

#endif