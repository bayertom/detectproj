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


#include "ArithmeticParser.h"


//List of variables
const char * vars[] =
{
        "x",
        "y",
        "R",
        "a",
        "b",
        "c",
        "u",
        "v",
        "lat",
        "lon",
        "phi",
        "lam",
        "u0",
        "lat0",
        "phi0",
        "u1",
        "lat1",
        "phi1",
        "u2",
        "lat2",
        "phi2",
        "v0",
        "lon0",
        "lam0",
	"theta"
};


//List of constants
const char * consts [] =
{
        "pi",
        "Pi",
        "PI",
        "e",
        "E",
        "Ro",
        "RO"
};


//List of functions
const char * functs [] =
{
        "sin",
        "cos",
        "tan",
        "tg",
        "cot",
        "cotg",
        "asin",
        "acos",
        "atan",
        "ln",
        "log",
        "exp",
        "sqr",
        "sqrt",
        "abs",
        "sign"
};


//List of operators
const char * opers[] =
{
	"+",
	"-",
	"*",
	"/",
	"^",
	"_"
};


void ArithmeticParser::init(TVarConsFunctMap & vars_list, TVarConsFunctMap & consts_list, TVarConsFunctMap & functs_list, TVarConsFunctMap & opers_list)
{
        //Add variables, constants, functions to the map to enable the search

        //Variables
        vars_list[vars[v_x]] = v_x;
        vars_list[vars[v_y]] = v_y;
        vars_list[vars[v_R]] = v_R;
        vars_list[vars[v_a]] = v_a;
        vars_list[vars[v_b]] = v_b;
        vars_list[vars[v_c]] = v_c;
        vars_list[vars[v_u]] = v_u;
        vars_list[vars[v_v]] = v_v;
        vars_list[vars[v_lat]] = v_lat;
        vars_list[vars[v_lon]] = v_lon;
        vars_list[vars[v_phi]] = v_phi;
        vars_list[vars[v_lam]] = v_lam;
        vars_list[vars[v_u0]] = v_u0;
        vars_list[vars[v_lat0]] = v_lat0;
        vars_list[vars[v_phi0]] = v_phi0;
        vars_list[vars[v_u1]] = v_u1;
        vars_list[vars[v_lat1]] = v_lat1;
        vars_list[vars[v_phi1]] = v_phi1;
        vars_list[vars[v_u2]] = v_u2;
        vars_list[vars[v_lat2]] = v_lat2;
        vars_list[vars[v_phi2]] = v_phi2;
        vars_list[vars[v_v0]] = v_v0;
        vars_list[vars[v_lon0]] = v_lon0;
        vars_list[vars[v_lam0]] = v_lam0;
	vars_list[vars[v_theta]] = v_theta;

        //Constants
        consts_list[consts[c_pi]] = c_pi;
        consts_list[consts[c_Pi]] = c_Pi;
        consts_list[consts[c_PI]] = c_PI;
        consts_list[consts[c_e]] = c_e;
        consts_list[consts[c_pi]] = c_E;
        consts_list[consts[c_Ro]] = c_Ro;
        consts_list[consts[c_RO]] = c_RO;

        //Functions
        functs_list[functs[f_sin]] = f_sin;
        functs_list[functs[f_cos]] = f_cos;
        functs_list[functs[f_tan]] = f_tan;
        functs_list[functs[f_tg]] = f_tg;
        functs_list[functs[f_cot]] = f_cot;
        functs_list[functs[f_cotg]] = f_cotg;
        functs_list[functs[f_asin]] = f_asin;
        functs_list[functs[f_acos]] = f_acos;
        functs_list[functs[f_atan]] = f_atan;
        functs_list[functs[f_ln]] = f_ln;
        functs_list[functs[f_log]] = f_log;
        functs_list[functs[f_exp]] = f_exp;
        functs_list[functs[f_sqr]] = f_sqr;
        functs_list[functs[f_sqrt]] = f_sqrt;
        functs_list[functs[f_abs]] = f_abs;
        functs_list[functs[f_sign]] = f_sign;

	//Operators
	opers_list[opers[o_plus]] = o_plus;
	opers_list[opers[o_minus]] = o_minus;
	opers_list[opers[o_multiply]] = o_multiply;
	opers_list[opers[o_divide]] = o_divide;
	opers_list[opers[o_power]] = o_power;
	opers_list[opers[o_unary_minus]] = o_unary_minus;
}


void ArithmeticParser::infixToPostfix ( const char * infix, char * postfix)
{
        //Convert an equation from thr infix notation to the postfix notation
	//Unary - is represented by !
        TSignOperatorType plus_minus_oper_next = UnaryOperator;
	TSignOperatorTypes sign_types_postfix;

	//Stack with operators, implemented as vector
        std::stack < std::string, std::vector<std::string> > operators;

        //List of variables, constants, functions: initialize
        TVarConsFunctMap vars_list, consts_list, functs_list, opers_list;
        init ( vars_list, consts_list, functs_list, opers_list );

        //Stack of binary or unary +- operators, implemented as vector
	std::stack <TSignOperatorType, std::vector <TSignOperatorType> > sign_types_temp;

	//Empty equation return
	if (infix == NULL)
	{
		return;
	}

        //Process all characters of the infix equation
        while ( *infix != '\0' )
        {
                //Space, tab, new line: jump over
		if ((*infix == '\t') || (*infix == ' ') || (*infix == '\n'))
                {
			infix++;
			continue;
                }

		//Get valid sequence of chars: number, variable, operator, bracket, function
		char token[32];
		findToken(&infix, token);

                //Number
                if ( isdigit ( ( unsigned char ) *token ) )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, token );

                        //Move pointer about length of the token
                        postfix += strlen ( token );

                        //Set next operator +- as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Variable or const
                else if ( vars_list.find ( std::string ( token ) ) != vars_list.end() || consts_list.find ( std::string ( token ) ) != consts_list.end() )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, token );

                        //Move pointer about length of the token
                        postfix += strlen ( token );

                        //Set next operator +- as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Left bracket
                else if ( *token == '(' )
                {
                        //Add to the stack
                        operators.push ( std::string ( token ) );

                        //Set next operator +- as unary
                        plus_minus_oper_next = UnaryOperator;
                }

                //Right bracket
                else if ( *token == ')' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                char c[32];

                                //Run until operator with lowest priority found
                                do
                                {
                                        //Error, stop parsing
                                        if ( operators.empty() )
                                        {
                                                //Miss "(" bracket;
                                                throw ParseException ( "ParseException: can not parse equation,", "missing ( ." );
                                        }

                                        //New operator on the top of the stack
                                        strcpy ( c, operators.top().c_str() );

                                        //Type of the plus minus operator on the top of the stack
                                        TSignOperatorType plus_minus;
                                        plus_minus = ( !sign_types_temp.empty() ? sign_types_temp.top() : UnaryOperator );

                                        //Add new operator to the postfix notation
                                        if ( *c != '(' ) //Do not add bracket
                                        {
                                                //Add space before function string
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add + - type operator to the output and remove from the stack
						if ((*c == '+') || (*c == '-') || (*c == '_'))
                                                {
                                                        //Add + - operator type to the output
                                                        sign_types_postfix.push_back ( plus_minus );

                                                        //Remove + - operator type on the top of the stack
                                                        sign_types_temp.pop();
                                                }

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );
                                        }

                                        //Remove operator on the top of the stack
                                        operators.pop();

                                }
                                while ( *c != '(' ) ;
                        }

                        //Error
                        else
                        {
                                //Miss first "(" bracket;
                                throw ParseException ( "ParseException: can not parse equation, ", " missing ( ." );
                        }

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }


                //***************************************************************************

                //Functions (Priority Level 4)
                else if ( functs_list.find ( std::string ( token ) ) != functs_list.end() )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();
					
                                        //Found operator with lower priority, stop
					if ((*c == '^') || (*c == '*') || (*c == '/') || (*c == '+') || (*c == '-') || (*c == '(') || (*c == '_'))
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //Add token to the stack
                        operators.push ( token );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Power (Priority level 3)
                else if ( *token == '^' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Found operator with lower priority, stop
					if ((*c == '*') || (*c == '/') || (*c == '+') || (*c == '-') || (*c == '(') || (*c == '_'))
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or empty stack
                        operators.push ( std::string ( token ) );

                        //Set next operator + - as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Multiple, divide (Priority level 2)
                else if ( ( *token == '*' ) || ( *token == '/' ) )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Found operator with lower priority, stop
					if ((*c == '+') || (*c == '-') || (*c == '(') || (*c == '_') )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or empty stack
                        operators.push ( std::string ( token ) );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Add, subtract (Priority level 1)
		else if ((*token == '+') || (*token == '-') || (*token == '_'))
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Type of the plus minus operator on the top of the stack
                                        TSignOperatorType plus_minus;
                                        plus_minus = ( !sign_types_temp.empty() ? sign_types_temp.top() : UnaryOperator );

                                        //Found operator with lower priority
                                        if ( *c == '(' )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notationnotation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Add plus minus type to the output and remove from the stack
						if ((*c == '+') || (*c == '-') || (*c == '_'))
                                                {
                                                        //Remove plus minus type of the operator from the top of the stack
                                                        sign_types_temp.pop();

                                                        //Add plus minus the type to the output
                                                        sign_types_postfix.push_back ( plus_minus );
                                                }

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }
			
                        //No operator with lowest priority found or stack empty
			//Unary -
			if (plus_minus_oper_next == UnaryOperator )
				operators.push("_");

			//Binary -
			else
				operators.push ( std::string ( token ) );

                        //Set the next + - operator as unary or binary
                        sign_types_temp.push ( plus_minus_oper_next );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Throw exception
                else
                {
                        //Unknown variable or function
                        throw ParseException ( "ParseException: can not parse equation, unknown function or variable. ", token );
                }
        }

        //Add rest of the stack to the postfix
        while ( ! operators.empty() )
        {
                //Get operator on the top
                char c[32];

                //New operator on the top of the stack
                strcpy ( c, operators.top().c_str() );

                //Type of the plus minus operator on the top of the stack
                TSignOperatorType plus_minus;
                plus_minus = ( !sign_types_temp.empty() ? sign_types_temp.top() : UnaryOperator );

                //Add new operator to the postfix notation
                if ( *c != '(' )
                {
                        //Add space before function string
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add new operator to the postfix notation
                        strcat ( postfix, c );
                        postfix += strlen ( c );

                        //Add plus minus type to the output and remove from the stack
			if ((*c == '+') || (*c == '-') || (*c == '_'))
                        {
                                //Remove plus minus type of the operator from the top of the stack
                                sign_types_postfix.push_back ( plus_minus );

                                //Remove type of the +- operator on the top of the stack
                                sign_types_temp.pop();
                        }
                }

                //Error, stop
                else
                {
                        //Miss ")" bracket
                        throw ParseException ( "ParseException: can not parse equation, ", "misssing ) ." );
                        break;
                }

                //Remove operator on the top of the stack
                operators.pop();

        }

        //End postfix notation
        *postfix = '\0';
}



void ArithmeticParser::findToken ( const char ** equation, char * token )
{
        //Find first possible valid sequence of characters in the infix notation
        //	0 = space .................. NA (non-algebraic)
        //	1 = text/number.............. A (algebraic)
        //	2 = decimal separator ...... NA
        //	3 = arithmetic operator .... NA
	//      4 = bracket .................NA
        do
        {
                //Actual char and increment equation char
                const int c = * ( * equation ) ++;

                //Compare actual char end next char
                //Alphanumeric char followed by non-alphanumeric char or vice versa (digit / letter and something)
                //Detect changes: [1->0], [1->2], [1->3], [1->4], [2->1], [3->1], [4->1]
                if ( ( bool ) isalnum ( ( unsigned char ) ** equation ) != ( bool ) isalnum ( ( unsigned char ) c ) )
                {
			*token++ = c;

			//Detect changes: [1->2], [2->1], if nor c neither c++ is not the decimal separator, stop (the only accepted combination)
			if ((c != '.') && (c != ',') && (**equation != '.') && (**equation != ','))
			{
				break;
			}
                }

                //Two alphanumeric chars or two non-alphanumeric chars
                //Detect changes: [1->1], [3->0], [4-> 0], [3->4], [4->3], [3->3], [4->4]
                else
                {
                        // Detect changes [3->0], [4->0], [3->4], [4->3], [4->4], [3->3]
                        if ( ( c == '+' ) || ( c == '-' ) || ( c == '*' ) || ( c == '/' ) || ( c == '^' ) || ( c == '(' ) || ( c == ')' ) )
                        {
                                *token++ = c;
                                break;
                        }

                        //Detect changes [1->1]
                        else 
                        {
                                *token++ = c;
                        }
                }

        } while ( **equation != '\0' );

        //Correctly end char with \0
        *token = '\0';
}


TPostfixNotationDel ArithmeticParser::delimitPostfixNotation(char *equation_postfix)
{
	//Delimit postfix notation from char
	char *equation_delimit;
	TPostfixNotationDel equation_postfix_del;

	//Delimiters: space, tab, new line
	equation_delimit = strtok(equation_postfix, " \t\n");

	while (equation_delimit)
	{
		//Convert to string
		std::string postfix_del_string(equation_delimit);
			
		//Add to the list
		equation_postfix_del.push_back(postfix_del_string);
		
		equation_delimit = strtok(NULL, " \t\n");
	}

	return equation_postfix_del;
}
