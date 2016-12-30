// Description: Load text file, items will be denominated into cols

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
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef File_HPP
#define File_HPP

#include <cstring>
#include <string>
#include <stdlib.h> 

#include "libalgo/source/exceptions/BadDataException.h"

template <typename Point, typename List>
void File::load3DPoints(const char * file, List &l)
{
	//Load file and split to the words
	const TVector2D <std::string> file_content = loadFileToWords(file);

	//Convert to 3D points
	convertTo3DPoints<Point>(l, file_content);
}

template <typename Point, typename List>
void File::load2DPoints(const char * file, List &l)
{
	//Load file and split to the words
	const TVector2D <std::string> file_content = loadFileToWords(file);

	//Convert to 2D points
	convertTo2DPoints<Point>(l, file_content);
}

template <typename Point, typename List>
void File::load3DPointsD(const char * file, List &l)
{
	//Load file and split to the words
	const TVector2D <std::string> file_content = loadFileToWords(file);

	//Convert to 3D points
	convertTo3DPointsD<Point>(l, file_content);
}

template <typename Point, typename List>
void File::load2DPointsD(const char * file, List &l)
{
	//Load file and split to the words
	const TVector2D <std::string> file_content = loadFileToWords(file);

	//Convert to 2D points
	convertTo2DPointsD<Point>(l, file_content);
}

template <typename Point, typename List>
void File::convertTo3DPoints(List &l, const TVector2D <std::string> &file_content)
{
	//Convert loaded points to the list of 3D points
	for (unsigned int i = 0; i < file_content.size(); i++)
	{
		// 3D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
		if (file_content[i].size() == 8)
		{
			l.push_back(Point(file_content[i][0].c_str(), atof(file_content[i][1].c_str()) + atof(file_content[i][2].c_str()) / 60.0 + atof(file_content[i][3].c_str()) / 3600.0,
				atof(file_content[i][4].c_str()) + atof(file_content[i][5].c_str()) / 60.0 + atof(file_content[i][6].c_str()) / 3600.0, atof(file_content[i][7].c_str())));
		}

		//3D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
		else if (file_content[i].size() == 7)
		{
			l.push_back(Point(atof(file_content[i][0].c_str()) + atof(file_content[i][1].c_str()) / 60.0 + atof(file_content[i][2].c_str()) / 3600.0,
				atof(file_content[i][3].c_str()) + atof(file_content[i][4].c_str()) / 60.0 + atof(file_content[i][5].c_str()) / 3600.0, atof(file_content[i][6].c_str())));
		}

		//Common 3D point with label
		else if (file_content[i].size() == 4)
		{
			l.push_back(Point(file_content[i][0].c_str(), atof(file_content[i][1].c_str()), atof(file_content[i][2].c_str()), atof(file_content[i][3].c_str())));
		}

		//Common 3D point without label
		else if (file_content[i].size() == 3)
		{
			l.push_back(Point(atof(file_content[i][0].c_str()), atof(file_content[i][1].c_str()), atof(file_content[i][2].c_str())));
		}

		//Throw exception
		else throw BadDataException("BadDataException: unknown data format. ", "Can't read the file.");
	}
}	


template <typename Point, typename List>
void File::convertTo2DPoints(List &l, const TVector2D <std::string> &file_content)
{
	//Convert loaded points to the list of 2D points
	for (unsigned int i = 0; i < file_content.size(); i++)
	{
		// 2D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
		if (file_content[i].size() == 7)
		{
			l.push_back(Point(file_content[i][0].c_str(), atof(file_content[i][1].c_str()) + atof(file_content[i][2].c_str()) / 60.0 + atof(file_content[i][3].c_str()) / 3600.0,
				atof(file_content[i][4].c_str()) + atof(file_content[i][5].c_str()) / 60.0 + atof(file_content[i][6].c_str()) / 3600.0));
		}

		//2D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
		else if (file_content[i].size() == 6)
		{
			l.push_back(Point(atof(file_content[i][0].c_str()) + atof(file_content[i][1].c_str()) / 60.0 + atof(file_content[i][2].c_str()) / 3600.0,
				atof(file_content[i][3].c_str()) + atof(file_content[i][4].c_str()) / 60.0 + atof(file_content[i][5].c_str()) / 3600.0));
		}

		//Common 2D point with label
		else if (file_content[i].size() == 3)
		{
			l.push_back(Point(file_content[i][0].c_str(), atof(file_content[i][1].c_str()), atof(file_content[i][2].c_str())));
		}

		//Common 2D point without label
		else if (file_content[i].size() == 2)
		{
			l.push_back(Point(atof(file_content[i][0].c_str()), atof(file_content[i][1].c_str())));
		}

		//Throw exception
		else throw BadDataException("BadDataException: unknown data format. ", "Can't read the file.");
	}
}

template <typename Point, typename List>
void File::convertTo3DPointsD(List &l, const TVector2D <std::string> &file_content)
{
	//Convert loaded points to the list of 3D points (dynamically allocated)
	for ( unsigned int i = 0; i < file_content.size(); i++ )
        {
                // 3D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                if ( file_content[i].size() == 8 )
                {
                        l.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0, atof ( file_content[i][7].c_str() ) ) );
                }

                //3D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                else if ( file_content[i].size() == 7 )
                {
                        l.push_back ( new Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0, atof ( file_content[i][6].c_str() ) ) );
                }

                //Common 3D point with label
                else if ( file_content[i].size() == 4 )
                {
                        l.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ), atof ( file_content[i][3].c_str() ) ) );
                }

                //Common 3D point without label
                else if ( file_content[i].size() == 3 )
                {
                        l.push_back ( new Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                }

                        //Throw exception
                        else throw BadDataException ( "BadDataException: unknown data format. ", "Can't read the file." );
        }
}


template <typename Point, typename List>
void File::convertTo2DPointsD(List &l, const TVector2D <std::string> &file_content)
{
	//Convert loaded points to the list of 2D points (dynamically allocated)
        for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        // 2D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        if ( file_content[i].size() == 7 )
                        {
                                l.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0 ) );
                        }

                        //2D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        else if ( file_content[i].size() == 6 )
                        {
                                l.push_back ( new Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0 ) );
                        }

                        //Common 2D point with label
                        else if ( file_content[i].size() == 3 )
                        {
                                l.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                        }

                        //Common 2D point without label
                        else if ( file_content[i].size() == 2 )
                        {
                                l.push_back ( new Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ) ) );
                        }

                        //Throw exception
                        else throw BadDataException ( "BadDataException: unknown data format. ", "Can't read the file." );
                }
}

#endif