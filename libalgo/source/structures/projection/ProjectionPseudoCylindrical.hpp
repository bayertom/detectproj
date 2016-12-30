// Description: Pseudocylindrical projection, derived from Projection
// Supported equations in the non-closed form
//For evaluation, stored in the postfix notation

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

#ifndef ProjectionPseudoCylindrical_HPP
#define ProjectionPseudoCylindrical_HPP


template <typename T>
ProjectionPseudoCylindrical <T> ::ProjectionPseudoCylindrical(const T R_, const T lat0_, const T latp_, const T lonp_, const TTransformedLongtitudeDirection lon_dir_, const char * ftheta_equat_,
	const char * theta0_equat_, const TPostfixNotationDel * ftheta_equat_postfix_, const TPostfixNotationDel * theta0_equat_postfix_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, 
	const TPostfixNotationDel &x_equat_postfix_, const TPostfixNotationDel &y_equat_postfix_, const char * projection_family_, const char * projection_name_ ) : Projection <T>(R_, lon0_, dx_, dy_, c_, x_equat_, y_equat_, x_equat_postfix_, 
	y_equat_postfix_, projection_family_, projection_name_), lat0(lat0_), cart_pole(latp_, lonp_), lon_dir(lon_dir_), ftheta_equat_postfix(ftheta_equat_postfix_), theta0_equat_postfix(theta0_equat_postfix_)
{

	if (ftheta_equat != NULL)
	{
		ftheta_equat = new char[strlen(ftheta_equat_) + 1];
		strcpy(ftheta_equat, ftheta_equat_);
	}

	else
	{
		ftheta_equat = NULL;
	}

	if (theta0_equat_ != NULL)
	{
		theta0_equat_ = new char[strlen(theta0_equat_) + 1];
		strcpy(theta0_equat, theta0_equat_);
	}

	else
	{

		theta0_equat = NULL;
	}

}



template <typename T>
ProjectionPseudoCylindrical <T> ::ProjectionPseudoCylindrical(const ProjectionPseudoCylindrical <T> &proj) : Projection <T>(proj), lat0(proj.lat0), cart_pole(proj.cart_pole), lon_dir(proj.lon_dir), ftheta_equat_postfix(proj.ftheta_equat_postfix), theta0_equat_postfix(proj.theta0_equat_postfix)
{
	if (proj.ftheta_equat != NULL)
	{
		ftheta_equat = new char[strlen(proj.ftheta_equat) + 1];
		strcpy(ftheta_equat, proj.ftheta_equat);
	}

	else
	{
		ftheta_equat = NULL;
	}

	if (proj.theta0_equat != NULL)
	{
		theta0_equat = new char[strlen(proj.theta0_equat) + 1];
		strcpy(theta0_equat, proj.theta0_equat);
	}

	else
	{
		theta0_equat = NULL;
	}

}


template <typename T>
ProjectionPseudoCylindrical <T>:: ~ProjectionPseudoCylindrical()
{
	if (ftheta_equat != NULL)
	{
		delete[] ftheta_equat;
		ftheta_equat = NULL;
	}


	if (theta0_equat != NULL)
	{
		delete[] theta0_equat;
		theta0_equat = NULL;
	}
}


template <typename T>
void ProjectionPseudoCylindrical <T> ::setFThetaEquat(const char * ftheta_equat_)
{
	if (ftheta_equat_ != NULL)
	{
		ftheta_equat = new char[strlen(ftheta_equat_) + 1];
		strcpy(ftheta_equat, ftheta_equat_);
	}

	else
	{
		ftheta_equat = NULL;
	}
}


template <typename T>
void ProjectionPseudoCylindrical <T> ::setTheta0Equat(const char * theta0_equat_)
{
	if (theta0_equat_ != NULL)
	{
		theta0_equat = new char[strlen(theta0_equat_) + 1];
		strcpy(theta0_equat, theta0_equat_);
	}

	else
	{
		theta0_equat = NULL;
	}
}


template <typename T>
void ProjectionPseudoCylindrical <T> ::FThetaEquatToPostfix()
{
	//Convert the infix notation to the postfix notation
	char ftheta_equat_postfix_[MAX_TEXT_LENGTH];

	try
	{
		//Parse ftheta equation
		if (ftheta_equat != NULL)
		{
			//Infix to postfix
			char ftheta_equat_postfix_[MAX_TEXT_LENGTH];
			ArithmeticParser::infixToPostfix(ftheta_equat, ftheta_equat_postfix_);

			//Delimit postfix notation
			ftheta_equat_postfix = ArithmeticParser::delimitPostfixNotation(ftheta_equat_postfix_);
		}
	}

	//Throw exception
	catch (MathException <T> & error)
	{
		//Throw exception
		throw;
	}

	//Throw exception
	catch (Exception & error)
	{
		//Throw exception
		throw error;
	}
}


template <typename T>
void ProjectionPseudoCylindrical <T> ::Theta0EquatToPostfix()
{
	//Convert the infix notation to the postfix notation
	char theta0_equat_postfix_[MAX_TEXT_LENGTH];

	try
	{

		if (theta0_equat != NULL)
		{
			//Infix to postfix
			char theta0_equat_postfix_[MAX_TEXT_LENGTH];
			ArithmeticParser::infixToPostfix(theta0_equat, theta0_equat_postfix_);

			//Delimit postfix notation
			theta0_equat_postfix = ArithmeticParser::delimitPostfixNotation(theta0_equat_postfix_);
		}
	}

	//Throw exception
	catch (MathException <T> & error)
	{
		//Throw exception
		throw;
	}

	//Throw exception
	catch (Exception & error)
	{
		//Throw exception
		throw error;
	}
}

#endif