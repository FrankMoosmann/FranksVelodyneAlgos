/*!
    \file   ParamStructCreation.h
    \brief  Provides a preprocessor-macro for easy creation of a parameter-struct
    \author Frank Moosmann (<frank.moosmann@kit.edu>)
    \date   15.11.2011

    The big advantage is that each variable has to be specified only once
    including a description and its default value
    From this information a parameter struct is created containing all specified variables
    A default constructor is also created which initializes all variables as specified
    Additionally, an options_description (boost::program_options) is generated,
    which can be directly used when parsing arguments

    Below is a small example code on how to use this file

    Copyright: Frank Moosmann
*/

/*****************************************
            USAGE-EXAMPLE
******************************************
#include "ParamStructCreation.h"
#define PARAM_ARR (3, ( \     // THIS size value must be adjusted accordingly
       (double,      param1, 2.0, "description for first parameter", 0) \
      ,(int,         param2, 536, "description 2", 0) \
      ,(vector<int>, pvec,      , "vector parameter", 1) \
      ) )
CREATE_PARAM_STRUCT(teststruct, "options for test", PARAM_ARR)
CREATE_PARAM_STRUCT_STREAMOP(teststruct, PARAM_ARR, "\n", 1) \

teststruct ts;
boost::program_options::options_description progOpt;
progOpt.add_options()
  ("help", "produce help message")
  ("verbose", "specify for verbose output")
;
progOpt.add(ts.parameterDescription);
// missing: parse argv with the help of boost::program_options
std::cout << std::endl << "parameter: " << std::endl << ts << std::endl;
******************************************/

#ifndef PARAMSTRUCTCREATION_H_
#define PARAMSTRUCTCREATION_H_

#include <boost/program_options.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/array/size.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/control/if.hpp>

/****************************************
         "private" HELPER-MACROS
         don't use them directly
*****************************************/

// the following helper-macro extracts a specified field from the above array:
#define ARR_VAL(z, row, field, array) \
  BOOST_PP_TUPLE_ELEM( 5, field, BOOST_PP_ARRAY_ELEM( row, array ) )
// the following helper-macro generates the variable definition for the specified nth row in the array
// of the form         double myVariable;
#define GEN_DEFINITION(z, n, array) \
  ARR_VAL( z, n, 0, array ) ARR_VAL( z, n, 1, array ) ;
// the following helper-macro generates the variable initialization for the specified nth row in the array
// of the form         ,myVariable(-1)
#define GEN_INITIALIZER(z, n, array) \
  BOOST_PP_COMMA_IF(n) ARR_VAL( z, n, 1, array ) ( ARR_VAL( z, n, 2, array ) )
// the following helper-macro generates the program-options entry for the specified nth row in the array
// of the form         ("myVariable", po::value(&myVariable)->default_value(myVariable), "")
// or                  ("myVariable", po::value(&myVariable)->multitoken(), "")
#define GEN_OPTDESC(z, n, array) \
  (BOOST_PP_STRINGIZE(ARR_VAL( z, n, 1, array )), \
      po::value(&ARR_VAL( z, n, 1, array ))BOOST_PP_IF(ARR_VAL( z, n, 4, array ),->multitoken(),->default_value(ARR_VAL( z, n, 2, array ))), \
      ARR_VAL( z, n, 3, array ))
// the following helper-macro generates an outstream-command for the specified nth row in the array
// of the form         [<< sep] [<< "myVariable: "] << param.myVariable
#define GEN_OSTREAM(z, n, arr_sep_nflag) \
  BOOST_PP_IF(n, << BOOST_PP_TUPLE_ELEM(3,1,arr_sep_nflag), ) \
  BOOST_PP_IF(BOOST_PP_TUPLE_ELEM(3,2,arr_sep_nflag), << BOOST_PP_STRINGIZE(ARR_VAL( z, n, 1, BOOST_PP_TUPLE_ELEM(3,0,arr_sep_nflag) )) << ": " , ) \
  << paramvec.ARR_VAL( z, n, 1, BOOST_PP_TUPLE_ELEM(3,0,arr_sep_nflag) )



/****************************************
           "public" MACROS
    the only functions to be used
        from outside this file
*****************************************/

// this macro defines the parameter struct including a constructor and program-option-description
// name + description can be arbitrarily specified
// the array must contain rows of (type, value_name, init_value, description_string, multiparam_flag)
#define CREATE_PARAM_STRUCT(structname, sturctdescr, array) \
  struct structname { \
    boost::program_options::options_description parameterDescription; \
    BOOST_PP_REPEAT( BOOST_PP_ARRAY_SIZE( array ), GEN_DEFINITION, array ) \
    structname() : \
      parameterDescription(sturctdescr) \
      ,BOOST_PP_REPEAT( BOOST_PP_ARRAY_SIZE( array ), GEN_INITIALIZER, array ) \
    { \
      namespace po = boost::program_options; \
      parameterDescription.add_options() \
        BOOST_PP_REPEAT( BOOST_PP_ARRAY_SIZE( array ), GEN_OPTDESC, array ) \
        ; \
    } \
  };

// this macro optionally defines an operator<< for the parameter struct
// assumption: the struct "structname" was already defined by the above function on the same array
// separator: separator between variables
// nameflag: if not 0 the variable name is printed before each value
#define CREATE_PARAM_STRUCT_STREAMOP(structname, array, separator, nameflag) \
  std::ostream &operator<<(std::ostream &ostr, const structname &paramvec) { \
    ostr BOOST_PP_REPEAT( BOOST_PP_ARRAY_SIZE( array ), GEN_OSTREAM, (array, separator, nameflag) ) ; \
    return ostr; \
  };


#endif /* PARAMSTRUCTCREATION_H_ */
