// Description: Some formatting routines

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


#ifndef Format_HPP
#define Format_HPP


template <typename T>
void Format::modScientific ( const T number, char * number_text )
{
        //Modified scientific format (condensed): shorten mantisa to zero decimal places and remove non valid zeros from exponent
        //Example: 6378137.135 -> 6.38e+6
        sprintf ( number_text, "%.2e", number ) ;

        //Get exponent and jump the sign
        char * exponent = strchr ( number_text, 'e' ) + 2;

        //Trim exponent to one or two digits
        if ( ( exponent != NULL ) && ( exponent[0] == '0' ) )
        {
                exponent[0] = exponent[1];
                exponent[1] = exponent[2];
                exponent[2] = '\0';

                //Trim exponent to one digit
                if ( exponent[0] == '0' )
                {
                        exponent[0] = exponent[1];
                        exponent[1] = '\0';
                }
        }
}


template <typename T>
std::string Format::modScientific(const T number)
{
	//Modified scientific format (condensed): shorten mantisa to zero decimal places and remove non valid zeros from exponent
	//Example: 6378137.135 -> 6.38e+6
	char text[255];

	//Call method for the char
	modScientific(number, text);

	//Convert char to string
	return (std::string) text;
}

#endif
