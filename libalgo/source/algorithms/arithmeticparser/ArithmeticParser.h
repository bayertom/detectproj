// Description: Arithmetic parserm converting the infix notation to the postfix notation, based on modified shunting-yard algorithm
//Unary minus represented by _

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


#ifndef ArithmeticParser_H
#define ArithmeticParser_H


#include <ostream>
#include <iostream>
#include <vector>
#include <map>
#include <cstring>


//Extern declarations
extern const char * vars[];
extern const char * consts[];
extern const char * functs[];
extern const char * opers[];


//Type of the +- operator (Binary or unary)
typedef enum
{
        UnaryOperator = 0,
        BinaryOperator,
} TSignOperatorType;


//List of types of + - operators in postfix notation
typedef std::vector <TSignOperatorType> TSignOperatorTypes;


//Comparator for variables, constants and functions represented by strings
struct compVarConstFunctMap
{
        bool operator() ( const std::string & a, const std::string & b ) const
        {
               return std::strcmp ( a.c_str(), b.c_str() ) < 0;
        }
};


//Map of variables, functions, constants (for faster searching than array )
typedef std::map <std::string, unsigned int, compVarConstFunctMap> TVarConsFunctMap;


//Variables  in arithmetic expression
enum variables
{
        v_x = 0,
        v_y,
        v_R,
        v_a,
        v_b,
        v_c,
        v_u,
        v_v,
        v_lat,
        v_lon,
        v_phi,
        v_lam,
        v_u0,
        v_lat0,
        v_phi0,
        v_u1,
        v_lat1,
        v_phi1,
        v_u2,
        v_lat2,
        v_phi2,
        v_v0,
        v_lon0,
        v_lam0,
	v_theta
};


//Constants in arithmetic expression
enum constants
{
        c_pi = 0,
        c_Pi,
        c_PI,
        c_e,
        c_E,
        c_Ro,
        c_RO
};


//Functions in arithmetic expression
enum functions
{
        f_sin = 0,
        f_cos,
        f_tan,
        f_tg,
        f_cot,
        f_cotg,
        f_asin,
        f_acos,
        f_atan,
        f_ln,
        f_log,
        f_exp,
        f_sqr,
        f_sqrt,
        f_abs,
        f_sign,
};


//Operators in arithmetic expression
enum operators
{
	o_plus = 0,
	o_minus,
	o_multiply,
	o_divide,
	o_power,
	o_unary_minus
};


//Delimited postfix notation
typedef std::vector <std::string> TPostfixNotationDel;


//Arithmetic parser converting equation from infix to postfix notation based on modified Shunting-yard algorithm
class ArithmeticParser
{
        public:
                template <typename T>
		static T parseEquation(const char * equation, char ** equation_postfix, const T x, const bool print_exception = true, std::ostream * output = &std::cout);

                template <typename T>
		static T parseEquation(const char * equation, char ** equation_postfix, const T x, const T y, const bool print_exception = true, std::ostream * output = &std::cout);

                template <typename T>
		static T parseEquation(const char * equation, char ** equation_postfix, const T lat, const T lon, const T R, const T a, const T b, const T c, const T lat0, const T lat1, const T lat2, const T theta, const bool print_exception = true, std::ostream * output = &std::cout);

		template <typename T>
		static T parseEquation(const TPostfixNotationDel * equation_postfix, const T lat, const T lon, const T R, const T a, const T b, const T c, const T lat0, const T lat1, const T lat2, const T theta, const bool print_exception = true, std::ostream * output = &std::cout);

		static void infixToPostfix(const char * infix, char * postfix);

		template <typename T>
		static T evaluatePostfixEquation(const TPostfixNotationDel * equation_postfix, const T x, const T y, const T lat, const T lon, const T R, const T a, const T b, const T c, const T lat0, const T lat1, const T lat2, const T theta);

		static TPostfixNotationDel delimitPostfixNotation(char *equation_postfix);

        private:

		static void init(TVarConsFunctMap & vars_list, TVarConsFunctMap & consts_list, TVarConsFunctMap & functs_list, TVarConsFunctMap & opers_list );

                static void findToken ( const char ** equation, char * operator_text );
		
};


#include "ArithmeticParser.hpp"


#endif
